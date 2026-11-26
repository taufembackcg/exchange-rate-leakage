@echo off
rem Runs run_pipeline.R with the R installation found on this PC. Options are passed on,
rem for example:  run_pipeline.bat --cores=8 --designs=main
setlocal
cd /d "%~dp0"
set "RSCRIPT="
for /f "delims=" %%i in ('where Rscript.exe 2^>nul') do if not defined RSCRIPT set "RSCRIPT=%%i"
if not defined RSCRIPT for /d %%d in ("%ProgramFiles%\R\R-*") do if exist "%%d\bin\Rscript.exe" set "RSCRIPT=%%d\bin\Rscript.exe"
if not defined RSCRIPT (echo Rscript.exe not found. & pause & exit /b 1)
"%RSCRIPT%" --no-save --no-restore run_pipeline.R %*
pause
