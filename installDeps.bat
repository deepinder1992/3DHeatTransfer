@echo off

echo =====================================
echo Installing dependencies for Heat3D
echo Windows (MSYS2)
echo =====================================

where pacman >nul 2>nul
if %errorlevel% neq 0 (
    echo ERROR: MSYS2 not found.
    echo Please install MSYS2 from https://www.msys2.org
    pause
    exit /b 1
)

echo Updating MSYS2 packages...
pacman -Syu --noconfirm

echo Installing required packages...
pacman -S --noconfirm ^
    mingw-w64-x86_64-gcc ^
    mingw-w64-x86_64-cmake ^
    mingw-w64-x86_64-openmp ^
    git

echo.
echo NOTE:
echo - CUDA must be installed separately from NVIDIA installer
echo - Make sure nvcc is available in PATH
echo.

echo Done.
pause
