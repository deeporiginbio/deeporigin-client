import sys
import enum
import datetime
import threading


class LogLevel(enum.Enum):
    DEBUG = 1
    INFO = 2
    WARNING = 3
    ERROR = 4


class Logger:
    def __init__(self, log_level, log_dest):
        self.level = self.str_to_loglevel(log_level)
        self.file = open(log_dest, "a") if isinstance(log_dest, str) else log_dest if log_dest else sys.stdout
        self.lock = threading.Lock()
        self.depth = 0

    def log(self, level, message, depth=None):
        if self.level.value <= level.value:
            if depth is None:
                depth = self.depth

            indent = "  " * depth  # Using two spaces per depth level for indentation
            log_entry = f"{indent}: {message}"
            with self.lock:
                print(log_entry, file=self.file)
                self.file.flush()

    def log_info(self, message, depth=None):
        self.log(LogLevel.INFO, message, depth)

    def log_warning(self, message, depth=None):
        self.log(LogLevel.WARNING, message, depth)

    def log_error(self, message, depth=None):
        self.log(LogLevel.ERROR, message, depth)

    def log_debug(self, message, depth=None):
        self.log(LogLevel.DEBUG, message, depth)

    def str_to_loglevel(self, level_str):
        try:
            return LogLevel[level_str]
        except KeyError:
            raise ValueError("Invalid log level: {}".format(level_str))

    def close(self):
        self.file.close()

    def add_depth(self):
        self.depth += 1

    def sub_depth(self):
        if self.depth > 0:
            self.depth -= 1

    def get_current_date(self):
        cur_date = datetime.datetime.now()
        return cur_date.strftime("%Y-%m-%d %H:%M:%S")

    def get_state_info(self):
        d = self.get_current_date()
        return "date={}".format(d)


DEFAULT_LOGGER = Logger("WARNING", None)
