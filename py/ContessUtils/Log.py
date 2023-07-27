# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

import logging
from pathlib import Path


class Logger(logging.Logger):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLevel('DEBUG')
        self.console_handler : logging.Handler = logging.StreamHandler()
        self.console_handler.setLevel('INFO')
        self.addHandler(self.console_handler)

    def add_log_file(self, path):
        path = Path(path)
        path.parent.resolve().mkdir(exist_ok=True, parents=True)
        hdlr = logging.FileHandler(path, 'w')
        hdlr.setLevel("DEBUG")
        self.addHandler(hdlr)

logging.setLoggerClass(Logger)
logger : Logger = logging.getLogger('ContessUtils')
