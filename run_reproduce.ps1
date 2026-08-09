$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $repoRoot
python src/paper_grade_analysis.py
python src/generate_paper_figures.py
python src/build_updated_3d.py
