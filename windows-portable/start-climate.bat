@echo off
setlocal
cd /d "%~dp0"
start "Climate Workbench server" /D "%~dp0" "%~dp0climate.exe" ui --host 127.0.0.1 --port 8080
timeout /t 2 /nobreak >nul
start "" "http://127.0.0.1:8080"
echo Climate Workbench is opening in your browser.
echo Keep the server window open while you use it.
