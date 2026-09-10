@echo off
setlocal
cd /d "%~dp0"
set FRONTEND_DIR=%CD%
set REPO_ROOT=%CD%\..

echo == Symbology Studio ==

where python >nul 2>nul
if errorlevel 1 (
  echo Error: python was not found. Please install Python 3.10+ first. 1>&2
  exit /b 1
)

if not exist "%REPO_ROOT%\bootstrap.exe" (
  echo The bootstrap binary is missing; please build it first with: make
)

if not exist server\.venv (
  echo Setting up the Python environment (one-time)...
  python -m venv server\.venv
  server\.venv\Scripts\pip install --quiet --upgrade pip
  server\.venv\Scripts\pip install --quiet -r server\requirements.txt
)

if not exist web\dist (
  where npm >nul 2>nul
  if not errorlevel 1 (
    echo Building the web interface (one-time)...
    pushd web
    call npm install --no-audit --no-fund
    call npm run build
    popd
  ) else (
    echo Note: npm not found; skipping the prebuilt UI.
  )
)

set URL=http://127.0.0.1:8321
echo.
echo Starting the server at %URL%
echo Press Ctrl+C to stop.
start "" "%URL%"

cd server
.venv\Scripts\python run.py
