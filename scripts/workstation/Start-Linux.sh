#!/bin/bash
# Linux launcher for curator testing (no git checkout).
# Run: ./Start-Linux.sh   (or double-click if your file manager allows)
set -euo pipefail

# TODO(maintainer): rnacentral/r2dt:pr-219 is a temporary PR-preview tag --
# switch the default below to rnacentral/r2dt:latest once that PR merges and
# a release build includes the workstation command. Override without editing
# this file via R2DT_WORKSTATION_IMAGE.
IMAGE="${R2DT_WORKSTATION_IMAGE:-rnacentral/r2dt:pr-219}"
PORT="8765"
WS="${R2DT_WORKSPACE:-$HOME/.r2dt-workstation}"
URL="http://127.0.0.1:${PORT}/compare"
HOME_URL="http://127.0.0.1:${PORT}/"

pause_if_terminal() {
  if [[ -t 0 ]]; then
    echo
    read -r -p "Press Enter to close this window…"
  fi
}

open_browser() {
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$1" >/dev/null 2>&1 || true
  fi
}

echo "R2DT workstation (Compare) — Linux"
echo "Image: ${IMAGE}"
echo "Workspace: ${WS}"
echo

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker was not found."
  echo "Install Docker Engine or Docker Desktop for Linux, then run Start-Linux.sh again."
  echo "https://docs.docker.com/engine/install/"
  pause_if_terminal
  exit 1
fi

if ! docker info >/dev/null 2>&1; then
  echo "Docker is not running (or this user cannot talk to the daemon)."
  echo "Start Docker, and ensure your user can run 'docker info' without sudo."
  echo "Then run Start-Linux.sh again."
  pause_if_terminal
  exit 1
fi

if ! mkdir -p "${WS}"; then
  echo
  echo "Could not create the workspace folder: ${WS}"
  echo "Check that you have permission to write there, or set R2DT_WORKSPACE"
  echo "to a folder you can write to, then try again."
  pause_if_terminal
  exit 1
fi

if curl -sf -o /dev/null "${HOME_URL}"; then
  echo "Workstation is already running — opening Compare…"
  open_browser "${URL}"
  pause_if_terminal
  exit 0
fi

echo "Pulling image (first time can take several minutes)…"
if ! docker pull "${IMAGE}"; then
  echo
  echo "Could not download the R2DT image (${IMAGE})."
  echo "Check your network connection and try again. If this keeps happening,"
  echo "contact the person who sent you this folder — the image tag may need updating."
  pause_if_terminal
  exit 1
fi

echo
echo "Starting on ${URL}"
echo "Leave this window open. Press Ctrl+C to stop."
echo

(
  for _ in $(seq 1 60); do
    if curl -sf -o /dev/null "${HOME_URL}"; then
      open_browser "${URL}"
      exit 0
    fi
    sleep 0.5
  done
) >/dev/null 2>&1 &

# `|| true` so a failed/interrupted run (including a deliberate Ctrl+C stop)
# still reaches pause_if_terminal below instead of exiting the script early.
docker run --rm \
  -p "127.0.0.1:${PORT}:${PORT}" \
  -v "${WS}:/workspace" \
  -w /rna/r2dt \
  "${IMAGE}" \
  python3 r2dt.py workstation \
    --workspace /workspace \
    --port "${PORT}" \
    --bind 0.0.0.0 \
    --docker-image "${IMAGE}" || true

pause_if_terminal
