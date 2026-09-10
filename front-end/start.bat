@echo off
setlocal EnableExtensions DisableDelayedExpansion
cd /d "%~dp0"
if errorlevel 1 exit /b 1
set "PYTHONUTF8=1"

echo == Symbology Studio ==
echo For full calculations on Windows, use WSL as described in README.md.

where python >nul 2>nul
if errorlevel 1 (
  echo Error: Python was not found. Please install Python 3.10+ first. 1>&2
  exit /b 1
)

python -c "import sys; sys.exit(sys.version_info < (3, 10))"
if errorlevel 1 (
  echo Error: Python 3.10+ is required. 1>&2
  exit /b 1
)

if not exist "server\.venv\Scripts\python.exe" (
  echo Setting up the Python environment...
  python -m venv "server\.venv"
  if errorlevel 1 goto failed
)

rem Retry dependencies after an interrupted installation, including an existing venv.
"server\.venv\Scripts\python.exe" -m pip install --quiet -r "server\requirements.txt"
if errorlevel 1 goto failed

if not exist "web\dist\index.html" (
  call :build_web
  if errorlevel 1 goto failed
)

set "URL=http://127.0.0.1:8321"
echo.
echo Starting the server at %URL%
echo Press Ctrl+C to stop.
if not defined SYMBOLOGY_NO_BROWSER start "" "%URL%"

cd server
if errorlevel 1 goto failed
".venv\Scripts\python.exe" run.py
exit /b %errorlevel%

:build_web
where npm >nul 2>nul
if errorlevel 1 (
  echo Error: Node.js and npm are required to build the web interface. 1>&2
  exit /b 1
)
echo Building the web interface...
pushd web
if errorlevel 1 exit /b 1
call npm ci --no-audit --no-fund
if errorlevel 1 goto web_failed
call npm run build
if errorlevel 1 goto web_failed
if not exist "dist\index.html" goto web_failed
popd
exit /b 0

:web_failed
popd
exit /b 1

:failed
echo Error: Setup failed. Resolve the error above and run start.bat again. 1>&2
exit /b 1
