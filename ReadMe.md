Paper: Identification and correction of sample mix-ups in expression genetic data
---------------------------------------------------------------------------------

The full manuscript (with supplementary tables and figures) is
[here](http://www.biostat.wisc.edu/~kbroman/publications/samplemixups_wsupp.pdf).

The paper will appear at [G3](http://g3journal.org/content/early/2015/08/19/g3.115.019778.abstract),
and is also available at [arXiv](http://arxiv.org/abs/1402.2633).

The primary manuscript files are `samplemixups_nolegends.Rnw`
and `samplemixups_supp_nolegends.tex`.

The Perl script `add_legends.pl`
adds all of the legends, and then the `.Rnw` file is run through
[knitr](http://yihui.name/knitr/) to create a
LaTeX file, and the two LaTeX files are sent through `pdflatex` and
`xelatex`, respectively, to create PDFs.

The `Makefile` tells the full story.

The `Analysis/R` subdirectory has an
[asciidoc](http://www.methods.co.nz/asciidoc/) file for the analyses
in the work. That directory has its own `Makefile`.
