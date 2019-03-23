Paper: Identification and correction of sample mix-ups in expression genetic data
---------------------------------------------------------------------------------

[![doi badge](https://zenodo.org/badge/DOI/10.5281/zenodo.2603993.svg)](https://doi.org/10.5281/zenodo.2603993)

The full manuscript (with supplementary tables and figures) is
[here](http://www.biostat.wisc.edu/~kbroman/publications/samplemixups_wsupp.pdf).

The paper is available at [arXiv](http://arxiv.org/abs/1402.2633) and
as a formal journal article at
[G3](http://g3journal.org/content/early/2015/08/19/g3.115.019778.abstract):

> Broman KW, Keller MP, Broman AT, Kendziorski C, Yandell BS, Sen
> &#346;, Attie AD (2015) Identification and correction of sample
> mix-ups in expression genetic data: A case study.
> [G3](http://g3journal.org) 5:2177-2186
> [![PubMed](https://kbroman.org/pages/icons16/pubmed-icon.png)](https://www.ncbi.nlm.nih.gov/pubmed/26290572)
> [![pdf](https://kbroman.org/pages/icons16/pdf-icon.png)](http://www.g3journal.org/content/ggg/5/10/2177.full.pdf)
> [![data](https://kbroman.org/pages/icons16/data-icon.png)](http://bit.ly/B6BTBR)
> [![R/lineup software](https://kbroman.org/pages/icons16/R-icon.png)](https://github.com/kbroman/lineup)
> [![doi](https://kbroman.org/pages/icons16/doi-icon.png)](https://doi.org/10.1534/g3.115.019778)


**I'm still working out how to include the data, or links to the data
  (because they're about 2 GB compressed).
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

Cached intermediate results are available at figshare:
[`samplemixups_rcache.zip`](http://files.figshare.com/2219404/samplemixups_rcache.zip)
(This contains a bunch of `.RData` files that go in `Analysis/R/Rcache`.)

---

## Necessary tools

- [R](https://www.r-project.org)
- Perl
- Python 2.7
- [GNU make](https://www.gnu.org/software/make/)
- [Asciidoc](http://www.methods.co.nz/asciidoc/)
- R packages: knitr, qtl, broman, lineup, ascii, data.table, igraph,
  beeswarm, RColorBrewer

---

## To Do

- Do clean tests, with and without the intermediate files
