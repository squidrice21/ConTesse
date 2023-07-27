# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

import asyncio
import logging
from functools import partial

from ContessUtils.Log import logger
from ContessUtils.AtomicInteger import AtomicInteger

RUNNER_TAG = AtomicInteger()

# capture the output of a process


async def log_stream(level, stream, tag):
    lines = []
    while not stream.at_eof():
        data = await stream.readline()
        line = data.decode(encoding="utf-8", errors="ignore").rstrip()
        if line:
            logger.log(level, "[%s] %s", tag, line)
            lines.append(line)
    return lines

# Run a process in background. Does not block.


async def run_subprocess_in_bg(case, cmd, **kwargs):
    logger.info("Running %s `%s`", case, " ".join([str(c) for c in cmd]))
    runner_tag = RUNNER_TAG.increment_and_get()
    process = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        **kwargs,
    )

    task_code = asyncio.ensure_future(process.wait())
    task_out = asyncio.ensure_future(log_stream(
        logging.DEBUG, process.stdout, runner_tag))
    task_err = asyncio.ensure_future(log_stream(
        logging.DEBUG, process.stderr, runner_tag))

    code = await task_code
    out = await task_out
    err = await task_err

    return code, out, err

gather = asyncio.gather


async def run_func_in_bg(func, *args, **kwargs):
    return await asyncio.get_event_loop().run_in_executor(executor=None, func=partial(func, *args, **kwargs))


def run_async_entry_point(coroutine):
    return asyncio.get_event_loop().run_until_complete(coroutine)
