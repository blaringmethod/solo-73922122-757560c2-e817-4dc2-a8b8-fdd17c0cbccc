"""FastAPI application: bounded uploads, status, cancellation, downloads."""

from __future__ import annotations

import re
import sqlite3
import tempfile
import uuid
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, UploadFile, status
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field, field_validator

from . import config, db, worker
from .fastq import FastqFormatError
from .indexcheck import (
    INDEX_LEN,
    MAX_SAMPLES,
    MIN_SAMPLES,
    IndexCheckError,
    SampleSpec,
    run_index_check,
)
from .worker import UploadTooLarge, start_worker, task_dir

_ACGT_RE = re.compile(r"^[ACGT]{%d}$" % INDEX_LEN)


@asynccontextmanager
async def lifespan(_: FastAPI):
    db.init_db()
    worker.recover_on_startup()
    yield


app = FastAPI(
    title="Paired FASTQ Service",
    version="1.1.0",
    description=(
        "Upload shuffled paired-end FASTQ files (R1/R2, uncompressed, "
        "max 256 MiB each). The service re-pairs them by read id, delivers "
        "both mates in R1 order and preserves the original bytes. "
        "Succeeded tasks can be audited with the index-check endpoint: "
        "streaming dual 8-base ACGT index assignment with <=1 mismatch per "
        "end, ambiguity / suspected crossover / unassigned classification."
    ),
    lifespan=lifespan,
)


class TaskOut(BaseModel):
    id: str
    status: str
    stage: str
    progress: float
    message: str
    r1_name: str
    r2_name: str
    r1_bytes: int
    r2_bytes: int
    created_at: str
    updated_at: str
    download_r1: str | None = None
    download_r2: str | None = None


class IndexSampleIn(BaseModel):
    name: str = Field(
        ..., min_length=1, max_length=128, description="Sample name (unique in the sheet)"
    )
    r1_index: str = Field(
        ...,
        min_length=INDEX_LEN,
        max_length=INDEX_LEN,
        description=f"R1 index: exactly {INDEX_LEN} ACGT characters",
    )
    r2_index: str = Field(
        ...,
        min_length=INDEX_LEN,
        max_length=INDEX_LEN,
        description=f"R2 index: exactly {INDEX_LEN} ACGT characters",
    )

    @field_validator("r1_index", "r2_index")
    @classmethod
    def _acgt_only(cls, value: str) -> str:
        if not _ACGT_RE.match(value):
            raise ValueError("index must contain exactly 8 uppercase ACGT characters")
        return value

    @field_validator("name")
    @classmethod
    def _name_non_blank(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("sample name must not be blank")
        return value


class IndexCheckIn(BaseModel):
    samples: list[IndexSampleIn] = Field(
        ...,
        description=f"Sample sheet: {MIN_SAMPLES}..{MAX_SAMPLES} uniquely named samples",
    )

    @field_validator("samples")
    @classmethod
    def _count_and_unique_names(cls, value: list[IndexSampleIn]) -> list[IndexSampleIn]:
        if not (MIN_SAMPLES <= len(value) <= MAX_SAMPLES):
            raise ValueError(
                f"expected between {MIN_SAMPLES} and {MAX_SAMPLES} samples, "
                f"got {len(value)}"
            )
        names = [item.name for item in value]
        if len(set(names)) != len(names):
            raise ValueError("sample names must be unique")
        return value


class IndexSampleCount(BaseModel):
    name: str
    count: int


class IndexCheckOut(BaseModel):
    task_id: str
    total_pairs: int
    samples: list[IndexSampleCount]
    ambiguous: int
    suspected_cross_sample: int
    unassigned: int


def _serialize(row: sqlite3.Row) -> TaskOut:
    download_r1 = download_r2 = None
    if row["status"] == "succeeded" and row["result_r1"]:
        download_r1 = f"/tasks/{row['id']}/download/r1"
        download_r2 = f"/tasks/{row['id']}/download/r2"
    return TaskOut(
        id=row["id"],
        status=row["status"],
        stage=row["stage"],
        progress=round(row["progress"], 4),
        message=row["message"],
        r1_name=row["r1_name"],
        r2_name=row["r2_name"],
        r1_bytes=row["r1_bytes"],
        r2_bytes=row["r2_bytes"],
        created_at=row["created_at"],
        updated_at=row["updated_at"],
        download_r1=download_r1,
        download_r2=download_r2,
    )


async def _spool_bounded(upload: UploadFile, side: str, dest: Path) -> int:
    """Stream an upload to ``dest`` without holding it in memory.

    SpooledTemporaryFile rolls over to a real file after
    ``UPLOAD_SPOOL_THRESHOLD`` and the copy proceeds in fixed chunks; the
    per-file 256 MiB limit is enforced while streaming.
    """
    spool = tempfile.SpooledTemporaryFile(max_size=config.UPLOAD_SPOOL_THRESHOLD)
    total = 0
    try:
        while True:
            chunk = await upload.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > config.MAX_FILE_BYTES:
                raise UploadTooLarge(side, config.MAX_FILE_BYTES)
            spool.write(chunk)
        spool.rollover()
        spool.seek(0)
        with dest.open("wb") as out:
            while True:
                chunk = spool.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    finally:
        spool.close()
    return total


@app.post("/tasks", response_model=TaskOut, status_code=status.HTTP_202_ACCEPTED)
async def create_task(
    r1: UploadFile = File(..., alias="r1", description="R1 FASTQ, uncompressed, <= 256 MiB"),
    r2: UploadFile = File(..., alias="r2", description="R2 FASTQ, uncompressed, <= 256 MiB"),
) -> TaskOut:
    task_id = uuid.uuid4().hex
    base = task_dir(task_id)
    base.mkdir(parents=True, exist_ok=False)
    r1_dest = base / "input_R1.fastq"
    r2_dest = base / "input_R2.fastq"
    try:
        size1 = await _spool_bounded(r1, "R1", r1_dest)
        size2 = await _spool_bounded(r2, "R2", r2_dest)
    except UploadTooLarge as exc:
        worker.cleanup_task_dir(task_id)
        raise HTTPException(
            status_code=getattr(
                status, "HTTP_413_CONTENT_TOO_LARGE",
                status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            ),
            detail=str(exc),
        )
    except Exception:
        worker.cleanup_task_dir(task_id)
        raise
    finally:
        await r1.close()
        await r2.close()

    conn = db.connect()
    try:
        db.create_task(
            conn,
            task_id,
            r1_name=Path(r1.filename or "r1.fastq").name,
            r2_name=Path(r2.filename or "r2.fastq").name,
            r1_bytes=size1,
            r2_bytes=size2,
        )
        row = db.get_task(conn, task_id)
    finally:
        conn.close()

    worker.register_cancel_event(task_id)
    start_worker(task_id)
    return _serialize(row)


@app.get("/tasks/{task_id}", response_model=TaskOut)
def get_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        return _serialize(row)
    finally:
        conn.close()


@app.post("/tasks/{task_id}/cancel", response_model=TaskOut)
def cancel_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        previous = db.try_set_cancelling(conn, task_id)
        if previous is None:
            raise HTTPException(status_code=404, detail="task not found")
        row = db.get_task(conn, task_id)
        if row["status"] in ("succeeded", "failed", "cancelled"):
            # The task reached a terminal state first; tell the caller plainly
            # that it cannot be cancelled.
            raise HTTPException(
                status_code=409,
                detail=f"task is already {row['status']} and cannot be cancelled",
            )
        # Flip the cooperative flag and return without waiting for the worker.
        worker.request_cancel(task_id)
        return _serialize(row)
    finally:
        conn.close()


def _result_path(task_id: str, side: int) -> Path:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        if row["status"] != "succeeded":
            raise HTTPException(
                status_code=409,
                detail=f"task is {row['status']}; results are not published",
            )
        rel = row["result_r1"] if side == 1 else row["result_r2"]
        if not rel:
            raise HTTPException(status_code=404, detail="result file missing")
    finally:
        conn.close()

    base = task_dir(task_id).resolve()
    path = (base / rel).resolve()
    # Guard against path traversal even though the rel value is server-generated.
    if base not in path.parents or not path.is_file():
        raise HTTPException(status_code=404, detail="result file missing")
    return path


@app.get("/tasks/{task_id}/download/r1")
def download_r1(task_id: str) -> FileResponse:
    path = _result_path(task_id, 1)
    return FileResponse(
        path, media_type="text/plain", filename="R1.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


@app.get("/tasks/{task_id}/download/r2")
def download_r2(task_id: str) -> FileResponse:
    path = _result_path(task_id, 2)
    return FileResponse(
        path, media_type="text/plain", filename="R2.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


@app.post(
    "/tasks/{task_id}/index-check",
    response_model=IndexCheckOut,
    summary="索引核查：按双端 8 位 ACGT 索引归属读段对",
    description=(
        "读取成功任务配对结果中每条读段序列的前 8 位索引（R1/R2 各一个）。"
        "每端至多允许一次错配；按两端总错配最少的合法样本组合归属，"
        "最优组合并列计为歧义。无合法组合时，两端各自唯一指向不同样本计为"
        "疑似串样，其余计为未归属。逐条复核双端读段号、四行结构与序列/质量"
        "长度；结果缺失或读取损坏时整次报错、不返回任何部分统计，也不修改"
        "原文件。长读段按固定块流式核查，遵守 256 MiB 内存上限。"
    ),
)
def index_check(task_id: str, payload: IndexCheckIn) -> IndexCheckOut:
    samples = [
        SampleSpec(item.name, item.r1_index, item.r2_index) for item in payload.samples
    ]
    try:
        result = run_index_check(task_id, samples)
    except IndexCheckError as exc:
        raise HTTPException(status_code=exc.status_code, detail=exc.detail) from exc
    except FastqFormatError as exc:
        raise HTTPException(status_code=422, detail=f"index check failed: {exc}") from exc
    return IndexCheckOut(**result)


@app.get("/health")
def health() -> JSONResponse:
    return JSONResponse({"status": "ok"})
