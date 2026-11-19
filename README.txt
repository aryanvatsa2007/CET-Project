CL 202 PROJECT, TOPIC 2
A computational engine for the calculus of thermodynamics
=========================================================

Group and authors: [TO FILL: group number and all five members]

Folder contents
---------------
manuscript.pdf            the manuscript (Industrial & Engineering Chemistry Research format,
                          achemso class, two-column layout)
manuscript.tex            editable LaTeX source
references.bib            bibliography (BibTeX)
generated_numbers.tex     every number quoted in the manuscript (written by Codes/run_all.py)
generated_tables.tex      every table body of the manuscript (written by Codes/run_all.py)
fig2_... fig8_....pdf     figures (written by Codes/run_all.py); Figure 1 is drawn in the LaTeX source
toc_graphic.pdf           table-of-contents graphic required by the journal
manuscript.bbl            compiled bibliography (lets the source compile without BibTeX)
AUDIT_REPORT.txt          result of the consistency audit (Codes/tools/audit.py)
SUBMISSION_CHECKLIST.txt  what is done and what the group still has to do
README.txt                this file

Literature/               reference list with verification status, download checklist;
                          PDF copies of the references go here (to be added by the group)
Codes/                    the code package, its README.md, tests, tools and the website
Data/                     all input and output files (.txt, .csv) and DATA_INDEX.txt

To rebuild everything:    python Codes/run_all.py
To compile the paper:     pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript
Interactive mode:         cd Codes && python -m thermocalculus ask   (the engine prompts for missing data)
To run the tests:         cd Codes && python -m pytest -q
To run the audit:         python Codes/tools/audit.py
Website:                  open Codes/website/index.html, or publish that folder with GitHub Pages
