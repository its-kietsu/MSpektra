param([string]$Compiler = 'g++')
$ErrorActionPreference='Stop'
foreach($item in @(@('polymer.cpp','mspolymer.dll'),@('singly_charged.cpp','mspolymer_links.dll'),@('charge_states.cpp','mspolymer_charge5.dll'))){
    & $Compiler -std=c++17 -O2 -Wall -Wextra -Wpedantic -shared -static -static-libgcc -static-libstdc++ (Join-Path $PSScriptRoot $item[0]) -o (Join-Path $PSScriptRoot $item[1])
    if($LASTEXITCODE -ne 0){throw "Native polymer build failed: $($item[0])"}
}
# Copy the outputs to _portable/msengine and run release_tools/pe_details.py
# before packaging; numerical code is unchanged by resource stamping.
