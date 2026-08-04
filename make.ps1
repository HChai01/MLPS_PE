$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = if ($env:PYTHON) { $env:PYTHON } else { "python" }
$Venv = Join-Path $Root ".mlps_pe_runtime"
$Requirements = Join-Path $Root "python\toolkit\requirements-toolkit.txt"

& $Python -m venv $Venv
$VenvPython = Join-Path $Venv "Scripts\python.exe"
& $VenvPython -m pip install --upgrade pip
& $VenvPython -m pip install -r $Requirements
Write-Host "Built MLPS_PE for Windows."
Write-Host "Run with: $VenvPython $(Join-Path $Root 'python\toolkit\mlps_pe.py') --help"
