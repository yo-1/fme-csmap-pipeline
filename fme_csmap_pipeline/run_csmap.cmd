@echo off
rem Edit the following path to match your Conda installation.
call "C:\Users\YOUR_NAME\miniforge3\condabin\conda.bat" run --no-capture-output -n csmap python "%~dp0csmap_pipeline.py" --config "%~dp0config.json" %*
exit /b %errorlevel%
