#!/bin/bash

echo \\begin{table}[htbp]          > tab5.tex
echo \\begin{center}               >>tab5.tex
echo \\begin{tabular}{ccc}         >>tab5.tex
echo dim \& hgm \& exact-hgm\\\\ >>tab5.tex
echo \\hline                       >>tab5.tex

./a.out 11  3 1.0 >> tab5.tex
./a.out 11  4 1.0 >> tab5.tex
./a.out 11  5 1.0 >> tab5.tex
./a.out 11  6 1.0 >> tab5.tex
./a.out 11  7 1.0 >> tab5.tex
./a.out 11  8 1.0 >> tab5.tex
./a.out 11  9 1.0 >> tab5.tex
./a.out 11 10 1.0 >> tab5.tex

echo \\hline >> tab5.tex
echo \\end{tabular} >> tab5.tex
echo \\end{center} >> tab5.tex
echo \\caption{Comparison at specific parameters} >> tab5.tex
echo \\label{tab:feller} >> tab5.tex
echo \\end{table} >> tab5.tex
