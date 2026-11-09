"""HTML reports: one self-contained page per run, built only from the run's tables. base.py is the frame of the
pages; inspect.py and starts.py are the reports of the two commands."""
from .inspect import DOCS, inspect_report
from .starts import starts_report

__all__ = ["DOCS", "inspect_report", "starts_report"]
