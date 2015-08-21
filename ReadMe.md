Paper: Identification and correction of sample mix-ups in expression genetic data
---------------------------------------------------------------------------------

The full manuscript (with supplementary tables and figures) is
[here](http://www.biostat.wisc.edu/~kbroman/publications/samplemixups_wsupp.pdf).

The paper will appear at [G3](http://g3journal.org/content/early/2015/08/19/g3.115.019778.abstract),
and is also available at [arXiv](http://arxiv.org/abs/1402.2633).

**I'm still working out how to include the data, or links to the data
  (because they're _really_ big).
  At the moment, this isn't fully reproducible.**

The data are available at the
[Mouse Phenome Database](http://phenome.jax.org/db/q?rtn=projects/projdet&reqprojid=532),
though not in exactly the form used in this repository.

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

---

## To Do

- R script to grab data from MPD and convert to form needed for rest
  of the scripts

- Seem to rely on old annotations in some bits; fix this

- Create zip file with all of the intermediate files (put on figshare?)

- Do clean tests, with and without the intermediate files

- Rename some stuff to make it cleaner?
