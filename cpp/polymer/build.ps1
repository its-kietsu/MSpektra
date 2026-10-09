param([string]$Compiler = 'g++')
$ErrorActionPreference='Stop'
& $Compiler -std=c++17 -O2 -Wall -Wextra -Wpedantic -shared -static -static-libgcc -static-libstdc++ (Join-Path $PSScriptRoot 'polymer.cpp') -o (Join-Path $PSScriptRoot 'mspolymer.dll')
if($LASTEXITCODE -ne 0){throw 'Native polymer build failed'}
