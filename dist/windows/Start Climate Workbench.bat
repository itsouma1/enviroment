@echo off
rem Starts the local workbench and opens it in your default browser.
rem Everything runs on this computer; nothing is uploaded or shared.
setlocal
cd /d "%~dp0"

if not exist "%~dp0climate.exe" (
  echo climate.exe was not found next to this file.
  echo Keep this .bat in the same folder as climate.exe and the data folder.
  pause
  exit /b 1
)

if not exist "%~dp0data" (
  echo The data folder is missing, so the workbench would have nothing to show.
  echo Extract the whole zip, not only climate.exe.
  pause
  exit /b 1
)

echo Starting the climate workbench on http://127.0.0.1:8080 ...
echo Leave this window open while you use it; closing it stops the server.
echo.

rem Open the browser a few seconds after the server starts, in the background.
start "" /b cmd /c "timeout /t 3 /nobreak >nul & start "" http://127.0.0.1:8080"

rem Run the server in this window so its messages stay readable.
"%~dp0climate.exe" ui --port 8080

echo.
echo The workbench server stopped.
pause
