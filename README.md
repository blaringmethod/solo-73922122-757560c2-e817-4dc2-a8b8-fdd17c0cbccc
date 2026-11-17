# Paired FASTQ Re-pairing Service

测序湖入湖工具：接收顺序被打乱的双端未压缩 FASTQ（R1/R2），按读段号重新
配对，按 R1 顺序交付两份文件，原文逐字节保留。FastAPI + SQLite 实现，单文件
上限 256 MiB，进程内存上限 256 MiB，可处理总量超过内存的输入。

成功任务另提供**索引核查（index check）**：读段对按双端 8 位 ACGT 索引前缀
归属到样本表（2～32 个样本），输出各样本计数与异常分类，供质控员判断哪些
样本需要隔离。核查全程流式读取，长读段也不会撑破内存上限。

## 配对规则

- 读段号取头行首个空白前的字段，去掉开头 `@`，再去掉末尾的 `/1` 或 `/2`
  （例如 `@read-7/1 sample=x` → `read-7`）。
- 输出按 **R1 的原始顺序**排列；R2 按相同读段号对齐到同一位置。
- 四行记录原文（含头行描述、换行风格）原样写入结果，不做任何改写。

## 失败即整次失败，不开放部分结果

任一情况发生，任务整体标记为 `failed`/`cancelled`，已写入的暂存全部删除，
下载接口返回 409：

- 同一文件内读段号重复（含 `/1`、`/2` 归一化后重复）
- 孤端（某读段号只出现在一侧）
- 格式错误：非四行式、头行不以 `@` 开头、第三行不以 `+` 开头、空文件
- 序列与质量字符串长度不等

## 内存与大输入

- 上传经 `SpooledTemporaryFile` 落盘（超过 1 MiB 即写临时文件），分块复制，
  单文件流式校验 256 MiB 上限，超限返回 413。
- 解析按 512 条/批写入 SQLite（WAL、小 page cache、临时表走文件
  `temp_store=FILE`、关闭 mmap），任意时刻内存里只有一个批次。
- 配对发布使用 `JOIN ... ORDER BY a.ord` 流式游标，分批 `fetchmany` 写出。
- 因此输入总量可以远超 256 MiB；自测以约 393 MiB 总量（两份各 196.5 MiB、
  各 200 万条）验证，峰值 RSS 约 19 MiB（见下文自测结果）。

## 任务与取消语义

- `POST /tasks` 上传后立即返回 `202`，后台守护线程处理。
- `GET /tasks/{id}` 查询状态、阶段（upload/parse/validate/publish/done）与
  0–1 进度。
- `POST /tasks/{id}/cancel`：置取消标记后**立即返回**，不等待配对结束；
  工作线程在下一个批次边界中止。取消与完成之间用条件 UPDATE 仲裁：
  - 取消先到 → 任务最终 `cancelled`，结果永不发布，暂存目录删除；
  - 完成先到 → 取消请求得到 `409`，明确告知任务已完成、不可取消，结果保留。
- 取消后可立刻创建新任务；任务目录以 UUID 隔离，清理只删除本任务目录，
  不会误删新任务文件。
- 服务重启（启动钩子）：所有非终态任务标记 `failed`（原因
  `interrupted by service restart`）并清理暂存与记录；`succeeded` 任务的
  结果文件保留，仍可下载。

## 索引核查（index check）规则

`POST /tasks/{id}/index-check` 的请求体是一份样本表：2～32 个唯一样本名，
每个样本给 R1、R2 各一个恰好 8 位的大写 ACGT 索引。服务读取成功任务两份
结果文件中每条序列的**前 8 个碱基**作为观测索引，逐条按下列规则归属：

- 每端对每个样本计算 8 位索引的错配数，**每端至多允许 1 个错配**，超过即
  该端不指向该样本。
- 只考虑两端都指向**同一个样本**的合法组合，取两端总错配最少的组合归属。
- 最优总错配存在并列（多个样本同为最小）→ 计为**歧义（ambiguous）**，不
  归入任何样本。
- 不存在合法组合时：
  - 两端各自**唯一**指向一个样本、且两个样本不同 → **疑似串样**
    （suspected_cross_sample）；
  - 其余情况（一端无指向、一端多指向、两端指向同一但不合法等）→
    **未归属（unassigned）**。
- 返回读段对总数、每个样本的计数以及三类异常计数；所有计数之和等于读段
  对总数，每条读段对恰好出现一次，便于决定隔离哪些样本。

核查同时**逐条按原配对规则复核**：双端读段号必须一致、每条记录四行结构、
头行 `@`、第三行 `+`、序列与质量非空且等长、两端记录数相同。任一结果文件
缺失或读取中发现损坏，整次请求报错（422/500），不返回任何部分统计，也不
修改原结果文件（只读打开）。

长读段按 64 KiB 定长块流式扫描，序列/质量行只计量长度并保留前 8 位，头部/
`+` 行设 1 MiB 结构性上限；内存占用与行长度无关。自测以单条 100 MiB 行长、
两份各 800 MiB 的结果文件验证，核查子进程峰值 RSS 约 13.5 MiB。

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/tasks` | multipart 字段 `r1`、`r2`，返回任务对象（202） |
| GET  | `/tasks/{id}` | 状态/进度/消息/下载链接 |
| POST | `/tasks/{id}/cancel` | 请求取消（立即返回） |
| GET  | `/tasks/{id}/download/r1` | 下载配对后的 R1 |
| GET  | `/tasks/{id}/download/r2` | 下载配对后的 R2 |
| POST | `/tasks/{id}/index-check` | 提交样本表，返回双端索引归属与异常统计 |
| GET  | `/health` | 健康检查 |

索引核查请求与响应示例：

```jsonc
// POST /tasks/{id}/index-check
{
  "samples": [
    {"name": "sample-1", "r1_index": "AAAAAAAA", "r2_index": "CCCCCCCC"},
    {"name": "sample-2", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"}
  ]
}
// 200
{
  "task_id": "…",
  "total_pairs": 4,
  "samples": [{"name": "sample-1", "count": 2}, {"name": "sample-2", "count": 0}],
  "ambiguous": 0,
  "suspected_cross_sample": 1,
  "unassigned": 1
}
```

样本表非法（样本数不在 2～32、索引不是 8 位 ACGT、名称重复/空白）返回 422；
任务不存在返回 404；任务未成功返回 409。

启动后打开交互式文档：<http://localhost:8000/docs>，可直接在页面上演示
上传、取消、下载全流程。

## 运行（Docker Compose，内存上限 256 MiB）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 即交付要求的内存上限；数据保存在
命名卷 `fastq-data`（容器内 `/data`，含 SQLite 与各任务目录），重启后成功
结果仍可下载。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
FASTQ_DATA_DIR=./data .venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                 # 需要一个运行在 localhost:8000 的服务
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本会生成一对乱序 FASTQ，演示上传 → 轮询进度 → 下载校验，并额外演示一次
取消、一次错误输入，以及成功任务的索引核查（精确归属、单错配、疑似串样、
未归属与 422/409 守卫）。

## 自测

```bash
# 短测试（功能/取消/重启语义 + 取消与完成的仲裁 + 索引核查 + 流式扫描器）
bash scripts/run_tests.sh
# 或分别运行
.venv/bin/python -I tests/test_api.py          # 30 项：配对、整次失败、取消、重启恢复
.venv/bin/python -I tests/test_cancel_race.py  # 9 项：取消/完成条件 UPDATE 仲裁
.venv/bin/python -I tests/test_index_check.py  # 35 项：索引归属/歧义/串样/未归属/整次报错
.venv/bin/python -I tests/test_scanner.py      # 20 项：分块扫描、CRLF、截断与行长上限

# 内存压力测试 1（配对）：两份各 196.5 MiB（合计约 393 MiB，超过 256 MiB 内存
# 上限），子进程实测峰值 RSS 约 19 MiB，并逐条校验配对顺序与原文
.venv/bin/python -I tests/test_memory.py 2000000

# 内存压力测试 2（索引核查）：单条序列/质量行 100 MiB、两份结果各约 800 MiB
# （合计约 1.6 GiB），核查子进程峰值 RSS 约 13.5 MiB，统计结果逐条核对
.venv/bin/python -I tests/test_index_check_memory.py 4
```

## 目录结构

```
app/
  config.py      # 上限、路径、批大小等配置（环境变量可覆盖）
  db.py          # SQLite schema 与任务/记录访问
  fastq.py       # 流式解析、校验、JOIN 配对发布
  indexcheck.py  # 索引核查：分块流式扫描、双端索引归属与异常分类
  worker.py      # 后台线程、取消仲裁、启动恢复、任务目录清理
  main.py        # FastAPI 路由、有界流式上传与索引核查 API
tests/           # 功能自测、扫描器单测与两套内存压力测试
scripts/demo.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 遗留项 / 边界说明

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 服务为单实例 + 单 SQLite 文件设计；多实例横向扩展需要换用共享存储/队列。
- 没有内置任务 TTL 清理，成功结果会一直保留在数据卷中，需外部运维定期清理。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile 仅做
  了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证（含
  `/docs`、上传、取消、下载、重启恢复）。在装有 Docker 的环境执行
  `docker compose up --build` 即可，`mem_limit: 256m` 即内存上限。
