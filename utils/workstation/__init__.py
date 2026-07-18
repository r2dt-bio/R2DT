"""Local curator workstation: private compare dashboard + job runner.

See ``.ai/local-curator-workstation-plan.md``. Editing / export phases are
not wired yet — this package covers dashboard + generate from uploads.
"""

from utils.workstation.server import run_server

__all__ = ["run_server"]
