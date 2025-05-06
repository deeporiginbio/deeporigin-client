from contextlib import redirect_stderr, redirect_stdout
import io
import sys
import unittest
import unittest.mock
import warnings

import pytest

from deeporigin import __version__, cli
from deeporigin.exceptions import DeepOriginException
from deeporigin.warnings import DeepOriginWarning


class TestCase(unittest.TestCase):
    def test_raw_help(self):
        with unittest.mock.patch("sys.argv", ["", "--help"]):
            with self.assertRaises(SystemExit) as context:
                cli.main()
                self.assertRegex(context.Exception, "usage: deeporigin")

    def test_help(self):
        with cli.App(argv=["--help"]) as app:
            with self.assertRaises(SystemExit) as context:
                app.run()
                self.assertRegex(context.Exception, "usage: deeporigin")

        with cli.App(argv=[]) as app:
            with self.assertRaises(SystemExit) as context:
                app.run()
                self.assertRegex(context.Exception, "usage: deeporigin")

    def test_version(self):
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["--version"]) as app:
                app.run()

        stdout = stdout_capture.getvalue().strip()

        self.assertEqual(__version__, stdout)

    def test_data(self):
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["data"]) as app:
                app.run()

        stdout = stdout_capture.getvalue().strip()

        self.assertRegex(stdout, "List data in the data hub on Deep Origin")

    @pytest.mark.skipif(
        sys.version_info < (3, 11), reason="requires Python 3.10 or higher"
    )
    def test_except_hook(self):
        mock_code = MockCode("mock_filename.py", "mock_function")
        mock_frame = MockFrame(mock_code, {})
        mock_tb = MockTraceback([mock_frame], [1])

        cli.except_hook(sys.excepthook, DeepOriginException, "error", mock_tb)
        cli.except_hook(sys.excepthook, Exception, "error", mock_tb)

    def test_format_warning(self):
        cli.format_warning(
            warnings.formatwarning, "message", DeepOriginWarning, "filename", 1
        )
        cli.format_warning(
            warnings.formatwarning, "message", UserWarning, "filename", 1
        )


class MockTraceback:
    def __init__(self, frames, line_nums):
        if len(frames) != len(line_nums):
            raise ValueError("Ya messed up!")
        self._frames = frames
        self._line_nums = line_nums
        self.tb_frame = frames[0]
        self.tb_lineno = line_nums[0]

    @property
    def tb_next(self):
        if len(self._frames) > 1:
            return MockTraceback(self._frames[1:], self._line_nums[1:])


class MockFrame:
    def __init__(self, f_code, f_globals):
        self.f_code = f_code
        self.f_globals = f_globals


class MockCode(object):
    def __init__(self, co_filename, co_name):
        self.co_filename = co_filename
        self.co_name = co_name
