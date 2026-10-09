param([string]$Compiler = 'g++')
$ErrorActionPreference = 'Stop'
$sourceRoot = $PSScriptRoot
if (-not (Get-Command $Compiler -ErrorAction SilentlyContinue)) { throw "C++ compiler missing: $Compiler" }
$output = Join-Path $sourceRoot 'msmaxent.dll'
& $Compiler -std=c++20 -O2 -Wall -Wextra -Wpedantic -shared -static -static-libgcc -static-libstdc++ '-I' (Join-Path $sourceRoot 'include') (Join-Path $sourceRoot 'src\core.cpp') (Join-Path $sourceRoot 'src\maxent_native.cpp') '-o' $output
if ($LASTEXITCODE -ne 0) { throw 'Native MaxEnt pipeline build failed' }
Write-Output "Built $output"
