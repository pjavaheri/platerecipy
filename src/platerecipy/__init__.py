
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ for package-wide logging ~~~~~~~~~~~~~~~~~~~~~~~~~~
import logging
import logging.config

logging.config.dictConfig({
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "default",
            "level": "DEBUG",
        },
    },
    "loggers": {
        "platerecipy": {  
            "level": "DEBUG", #"INFO", # "DEBUG"
            "handlers": ["console"],
            "propagate": False, # to prevent it from propagating through root
        },
    },
    "root": {
        "level": "WARNING",  # Default level for other packages
        "handlers": ["console"],
    },
})
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~ for consistent numerical types ~~~~~~~~~~~~~~~~~~~~~~~
from numpy import int32, float64

_INT     = int32
_FLOAT   = float64
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
_IS_WINDOWS = (sys.platform == 'win32')