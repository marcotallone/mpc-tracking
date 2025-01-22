#!/bin/bash

# Specify the extensions of the files you want to delete
extensions=(".aux" ".log" ".out" ".toc" ".lof" ".lot" ".fls" ".blg" ".bbl" ".bcf" ".fdb_latexmk" ".synctex.gz"".listing""-blx.bib" ".run.xml")

# Walk through the current directory and its subdirectories
for ext in "${extensions[@]}"; do
    # Find files with the current extension and delete them
    find . -type f -name "*.$ext" -exec rm -v {} \;
done