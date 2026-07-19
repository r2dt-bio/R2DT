@echo off
REM Windows launcher for curator testing (no git checkout).
REM Double-click Start-Windows.bat
setlocal EnableExtensions

REM TODO(maintainer): rnacentral/r2dt:pr-219 is a temporary PR-preview tag --
REM switch the default below to rnacentral/r2dt:latest once that PR merges and
REM a release build includes the workstation command. Override without
REM editing this file via the R2DT_WORKSTATION_IMAGE environment variable.
if defined R2DT_WORKSTATION_IMAGE (
  set "IMAGE=%R2DT_WORKSTATION_IMAGE%"
) else (
  set "IMAGE=rnacentral/r2dt:pr-219"
)
set "PORT=8765"
if defined R2DT_WORKSPACE (
  set "WS=%R2DT_WORKSPACE%"
) else (
  set "WS=%USERPROFILE%\.r2dt-workstation"
)
set "URL=http://127.0.0.1:%PORT%/compare"
set "HOME_URL=http://127.0.0.1:%PORT%/"

echo R2DT workstation (Compare) - Windows
echo Image: %IMAGE%
echo Workspace: %WS%
echo.

where docker >nul 2>&1
if errorlevel 1 (
  echo Docker was not found.
  echo Install Docker Desktop from https://www.docker.com/products/docker-desktop/
  echo then start it and double-click Start-Windows.bat again.
  echo.
  pause
  exit /b 1
)

docker info >nul 2>&1
if errorlevel 1 (
  echo Docker is not running.
  echo Open Docker Desktop, wait until it is ready, then double-click Start-Windows.bat again.
  echo.
  pause
  exit /b 1
)

if not exist "%WS%" mkdir "%WS%"

curl -sf -o NUL "%HOME_URL%" >nul 2>&1
if not errorlevel 1 (
  echo Workstation is already running - opening Compare...
  start "" "%URL%"
  echo.
  pause
  exit /b 0
)

echo Pulling image (first time can take several minutes)...
docker pull %IMAGE%
if errorlevel 1 (
  echo.
  echo Could not download the R2DT image ^(%IMAGE%^).
  echo Check your network connection and try again. If this keeps happening,
  echo contact the person who sent you this folder - the image tag may need updating.
  echo.
  pause
  exit /b 1
)

echo.
echo Starting on %URL%
echo Leave this window open. Press Ctrl+C to stop.
echo.

start "" /b powershell -NoProfile -Command "for ($i=0; $i -lt 60; $i++) { try { Invoke-WebRequest -UseBasicParsing '%HOME_URL%' -TimeoutSec 1 | Out-Null; Start-Process '%URL%'; exit 0 } catch {} Start-Sleep -Milliseconds 500 }"

docker run --rm -p 127.0.0.1:%PORT%:%PORT% -v "%WS%":/workspace -w /rna/r2dt %IMAGE% python3 r2dt.py workstation --workspace /workspace --port %PORT% --bind 0.0.0.0 --docker-image %IMAGE%

echo.
pause
endlocal
