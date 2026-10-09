param([string]$Compiler = 'g++')
$ErrorActionPreference = 'Stop'
& $Compiler -std=c++17 -O2 -Wall -Wextra -Wpedantic -shared -static -static-libgcc -static-libstdc++ (Join-Path $PSScriptRoot 'kinetics.cpp') -o (Join-Path $PSScriptRoot 'mskinetics.dll')
if ($LASTEXITCODE -ne 0) { throw 'Native kinetics build failed' }
