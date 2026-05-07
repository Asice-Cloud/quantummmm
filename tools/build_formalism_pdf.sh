#!/usr/bin/env bash
set -euo pipefail
mkdir -p results
if command -v pdflatex >/dev/null 2>&1; then
  pdflatex -interaction=nonstopmode -output-directory=results docs/nonabelian_formalism.tex
  echo "Built results/nonabelian_formalism.pdf using pdflatex"
elif command -v pandoc >/dev/null 2>&1; then
  pandoc docs/nonabelian_formalism.tex -o results/nonabelian_formalism.pdf
  echo "Built results/nonabelian_formalism.pdf using pandoc"
else
  echo "No pdflatex or pandoc found. To build locally run: pdflatex -output-directory=results docs/nonabelian_formalism.tex"
  exit 1
fi
