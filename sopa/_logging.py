from __future__ import annotations

import logging

log = logging.getLogger(__name__)


class ColorFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    blue = "\x1b[36;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    prefix = "[%(levelname)s] (%(name)s)"
    suffix = "%(message)s"

    FORMATS = {
        logging.DEBUG: f"{grey}{prefix}{reset} {suffix}",
        logging.INFO: f"{blue}{prefix}{reset} {suffix}",
        logging.WARNING: f"{yellow}{prefix}{reset} {suffix}",
        logging.ERROR: f"{red}{prefix}{reset} {suffix}",
        logging.CRITICAL: f"{bold_red}{prefix}{reset} {suffix}",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def configure_logger(log: logging.Logger):
    log.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(ColorFormatter())

    log.addHandler(consoleHandler)
    log.propagate = False
