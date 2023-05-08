import logging
from numpy import set_printoptions as numpy_set_printoptions

defaultConfig = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "standard": {"format": "[%(levelname).5s] %(message)s"},
    },
    "handlers": {
        "default": {
            "level": "INFO",
            "formatter": "standard",
            "class": "logging.StreamHandler",
        },
        "core": {
            "level": "DEBUG",
            "filename": "LOG/core.log",
            "class": "logging.FileHandler",
            "mode": "w",
            "formatter": "standard",
        },
    },
    "loggers": {
        "": {"handlers": ["default"], "level": "ERROR", "propagate": True},
        "compMesh": {
            "handlers": ["core"],
            "level": "DEBUG",
            "propagate": False,
        },
        "varForm": {
            "handlers": ["core"],
            "level": "DEBUG",
            "propagate": False,
        },
    },
}