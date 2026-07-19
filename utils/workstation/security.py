"""Loopback Host/Origin checks for the local workstation HTTP API."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional
from urllib.parse import urlparse

# Loopback names accepted in Host / Origin / Referer for mutating requests.
_LOOPBACK_HOSTS = frozenset({"localhost", "127.0.0.1", "::1"})
_MUTATING_METHODS = frozenset({"POST", "PUT", "DELETE", "PATCH"})


def path_is_within(path: Path, base: Path) -> bool:
    """Return True if ``path`` is ``base`` or a descendant of ``base``.

    Uses ``Path.relative_to`` so a sibling directory that merely shares a
    string prefix (e.g. ``{base}-evil``) is rejected. Both paths are
    resolved before comparison.
    """
    try:
        path.resolve().relative_to(base.resolve())
        return True
    except ValueError:
        return False


def _hostname_from_host_header(host: str) -> str:
    """Extract the hostname from an HTTP Host header value."""
    host = host.strip()
    if not host:
        return ""
    if host.startswith("["):
        end = host.find("]")
        if end == -1:
            return ""
        return host[1:end].rstrip(".").lower()
    # Strip :port for IPv4 / hostname (never unbracketed IPv6 in Host).
    if host.count(":") == 1:
        host = host.rsplit(":", 1)[0]
    return host.rstrip(".").lower()


def _hostname_from_url(url: str) -> str:
    """Extract hostname from an Origin or Referer URL."""
    parsed = urlparse(url.strip())
    if parsed.scheme not in ("http", "https"):
        return ""
    return (parsed.hostname or "").rstrip(".").lower()


def local_mutating_request_error(  # pylint: disable=too-many-return-statements
    method: str, headers: Mapping[str, str]
) -> Optional[str]:
    """Return a rejection reason for unsafe cross-origin / rebound requests.

    State-changing routes require a loopback ``Host``. When ``Origin`` (or
    ``Referer``) is present, it must also be loopback. Missing Origin/Referer
    is allowed so local ``curl`` and similar tools keep working. This blocks
    casual localhost CSRF and DNS rebinding; it is not authentication.
    """
    if method.upper() not in _MUTATING_METHODS:
        return None
    host = headers.get("Host") or headers.get("host") or ""
    if not host:
        return "Missing Host header"
    if _hostname_from_host_header(host) not in _LOOPBACK_HOSTS:
        return "Host must be localhost or 127.0.0.1"
    origin = headers.get("Origin") or headers.get("origin")
    if origin is not None and origin != "":
        if origin.strip().lower() == "null":
            return "Origin not allowed"
        if _hostname_from_url(origin) not in _LOOPBACK_HOSTS:
            return "Origin not allowed"
        return None
    referer = headers.get("Referer") or headers.get("referer")
    if referer:
        if _hostname_from_url(referer) not in _LOOPBACK_HOSTS:
            return "Referer not allowed"
    return None
