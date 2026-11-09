"""HTML reports: one self-contained page per run, built only from the run's tables. base.py is the frame of the
pages; inspect.py is the report of `inspect`."""
from .inspect import DOCS, inspect_report

__all__ = ["DOCS", "inspect_report"]
