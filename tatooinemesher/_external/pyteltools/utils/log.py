"""
Minimal stub for the original ``pyteltools.utils.log`` module.

The upstream module wires ``coloredlogs`` and reads logging settings from
``pyteltools.conf``. Here we expose a plain ``new_logger(name)`` factory based
on the standard library so the vendored Serafin module gets a usable logger
without pulling in extra dependencies.
"""

import logging


def new_logger(name):
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
