### lists of figs, supp figs and supp table files
figs = Figs/fig1.pdf Figs/fig2.pdf Figs/fig3.eps Figs/fig4.eps Figs/fig5.eps Figs/fig6.eps Figs/fig7.eps Figs/fig8.eps Figs/fig9.eps

suppfigs = SuppFigs/figS1.eps SuppFigs/figS2.eps SuppFigs/figS3.jpg SuppFigs/figS4.eps SuppFigs/figS5.eps SuppFigs/figS6.eps SuppFigs/figS7.eps SuppFigs/figS8.eps SuppFigs/figS9.eps SuppFigs/figS10.eps SuppFigs/figS11.eps SuppFigs/figS12.eps SuppFigs/figS13.eps SuppFigs/figS14.eps SuppFigs/figS15.eps SuppFigs/figS16.eps SuppFigs/figS17.eps SuppFigs/figS18.eps SuppFigs/figS19.eps SuppFigs/figS20.eps

supptabs = SuppTables/tableS1.tex SuppTables/tableS2.tex SuppTables/tableS3.tex SuppTables/tableS4.tex SuppTables/tableS5.tex

all: pdf Figs/cover.tiff

### pdf files
pdf: samplemixups.pdf samplemixups_supp.pdf samplemixups_wsupp.pdf

samplemixups_wsupp.pdf: samplemixups.pdf samplemixups_supp.pdf
	pdftk $^ cat output $@

### main manuscript
# add legends
ProcessedFiles/samplemixups.Rnw: samplemixups_nolegends.Rnw legends.txt
	add_legends.pl samplemixups_nolegends.Rnw ProcessedFiles/samplemixups.Rnw

# Sweave to .tex
ProcessedFiles/samplemixups.tex: ProcessedFiles/samplemixups.Rnw
	cd ProcessedFiles;Rscript -e 'library(knitr); knit("samplemixups.Rnw")'

# pdflatex to .pdf
samplemixups.pdf: ProcessedFiles/samplemixups.tex samplemixups.bib genetics.bst $(figs) $(suppfigs) $(supptabs)
	cd ProcessedFiles;pdflatex samplemixups
	cd ProcessedFiles;bibtex samplemixups
	cd ProcessedFiles;pdflatex samplemixups
	cd ProcessedFiles;pdflatex samplemixups
	cd ProcessedFiles;pdflatex samplemixups
	\mv ProcessedFiles/samplemixups.pdf .


# add legends to supplement
ProcessedFiles/samplemixups_supp.tex: samplemixups_supp_nolegends.tex legends.txt
	add_legends.pl samplemixups_supp_nolegends.tex ProcessedFiles/samplemixups_supp.tex

# supplement to pdf
samplemixups_supp.pdf: ProcessedFiles/samplemixups_supp.tex $(suppfigs) $(supptabs)
	cd ProcessedFiles;xelatex samplemixups_supp
	\mv ProcessedFiles/samplemixups_supp.pdf .

### supplemental tables
SuppTables/tableS1.tex: SuppTables/tableS1.Rnw Analysis/R/Rcache/dnadup.RData
	cd SuppTables;Rscript -e 'library(knitr); knit("tableS1.Rnw")'

SuppTables/tableS2.tex: SuppTables/tableS2.Rnw Analysis/R/Rcache/expr_tab.RData
	cd SuppTables;Rscript -e 'library(knitr); knit("tableS2.Rnw")'

SuppTables/tableS3.tex: SuppTables/tableS3.Rnw Analysis/R/Rcache/expr_corr_counts.RData
	cd SuppTables;Rscript -e 'library(knitr); knit("tableS3.Rnw")'

SuppTables/tableS4.tex: SuppTables/tableS4.Rnw Analysis/R/Rcache/agouti_tab.RData
	cd SuppTables;Rscript -e 'library(knitr); knit("tableS4.Rnw")'

SuppTables/tableS5.tex: SuppTables/tableS5.Rnw Analysis/R/Rcache/tufted_tab.RData
	cd SuppTables;Rscript -e 'library(knitr); knit("tableS5.Rnw")'

### figures
Figs/fig1.pdf: R/eve_scheme.R Analysis/OrigData/F2.mlratio.kidney.RData Analysis/OrigData/F2.mlratio.liver.RData Analysis/R/Rcache/expr_corr.RData Analysis/R/Rcache/expr_corr_betw_tissues.RData
	cd R;R CMD BATCH --no-save eve_scheme.R

Figs/fig2.pdf: R/gve_scheme.R Analysis/R/Rcache/dgve.RData
	cd R;R CMD BATCH --no-save gve_scheme.R

Figs/fig3.eps: R/eve_similarity.R Analysis/R/Rcache/expr_mixup_summaries.RData
	cd R;R CMD BATCH --no-save eve_similarity.R

Figs/fig4.eps: R/expr_swaps.R
	cd R;R CMD BATCH --no-save expr_swaps.R

Figs/fig5.eps: R/gve.R Analysis/R/Rcache/dgve.RData Analysis/R/Rcache/pmark.RData
	cd R;R CMD BATCH --no-save gve.R

Figs/fig6.eps: R/gve_similarity.R Analysis/R/Rcache/dgve_min_and_self.RData
	cd R;R CMD BATCH --no-save gve_similarity.R

Figs/fig7.eps: R/genotype_plates.R R/genotype_plates_arrows.csv Analysis/R/calls.RData
	cd R;R CMD BATCH --no-save genotype_plates.R

Figs/fig8.eps: R/insulin_lod.R R/my_plot_scanone.R Analysis/R/Rcache/insulin_scan.RData
	cd R;R CMD BATCH --no-save insulin_lod.R

Figs/fig9.eps: R/eqtl_counts.R Analysis/R/Rcache/neqtl_new.RData Analysis/R/Rcache/neqtl_old.RData
	cd R;R CMD BATCH --no-save eqtl_counts.R

### supplemental figures
SuppFigs/figS1.eps: R/xchr_fig.R
	cd R;R CMD BATCH --no-save xchr_fig.R

SuppFigs/figS2.eps: R/sex_vs_X.R Analysis/OrigData/rawg.RData Analysis/OrigData/clean_cross.RData Analysis/FinalData/necropsy.csv
	cd R;R CMD BATCH --no-save sex_vs_X.R

SuppFigs/figS3.jpg: R/expr_scatterplot_allprobes.R Analysis/R/Rcache/dgve.RData
	cd R;R CMD BATCH --no-save expr_scatterplot_allprobes.R

SuppFigs/figS4.eps: R/betw_tissue_corr.R Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save betw_tissue_corr.R

SuppFigs/figS5.eps: R/eve_hist.R Analysis/R/Rcache/expr_corr_betw_tissues.RData
	cd R;R CMD BATCH --no-save eve_hist.R

SuppFigs/figS6.eps: R/eve_similarity_supp.R Analysis/R/Rcache/expr_mixup_summaries.RData
	cd R;R CMD BATCH --no-save eve_similarity_supp.R

SuppFigs/figS7.eps: R/expr_scatterplots_swap.R Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save expr_scatterplots_swap.R

SuppFigs/figS8.eps: R/expr_scatterplots_dup.R Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save expr_scatterplots_dup.R

SuppFigs/figS9.eps: R/expr_corr_dup.R Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save expr_corr_dup.R

SuppFigs/figS10.eps: R/expr_scatterplots_mix.R Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save expr_scatterplots_mix.R

SuppFigs/figS11.eps: R/expr_corr_mix.R Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/expr_corr.RData
	cd R;R CMD BATCH --no-save expr_corr_mix.R

SuppFigs/figS12.eps: R/xist_and_y.R Analysis/OrigData/annot.final.RData Analysis/R/Rcache/mlratios.RData Analysis/R/Rcache/mlratios_revised.RData Analysis/FinalData/necropsy.csv
	cd R;R CMD BATCH --no-save xist_and_y.R

SuppFigs/figS13.eps: R/gve_supp.R Analysis/R/Rcache/dgve.RData Analysis/R/Rcache/pmark.RData
	cd R;R CMD BATCH --no-save gve_supp.R

SuppFigs/figS14.eps: R/local_eqtl_locations.R R/my_plot_map.R Analysis/FinalData/aligned_geno_with_pmap.RData Analysis/R/Rcache/dgve.RData
	cd R;R CMD BATCH --no-save local_eqtl_locations.R

SuppFigs/figS15.eps: R/gve_hist.R Analysis/R/Rcache/dgve.RData
	cd R;R CMD BATCH --no-save gve_hist.R

SuppFigs/figS16.eps: R/gve_similarity_alltissues.R Analysis/R/Rcache/dgve_min_and_self.RData
	cd R;R CMD BATCH --no-save gve_similarity_alltissues.R

SuppFigs/figS17.eps: R/gve_similarity_2ndbest.R Analysis/R/Rcache/dgve_min_and_self.RData
	cd R;R CMD BATCH --no-save gve_similarity_2ndbest.R

SuppFigs/figS18.eps: R/gve_new.R
	cd R;R CMD BATCH --no-save gve_new.R

SuppFigs/figS19.eps: R/coatcolor_lod.R R/my_plot_scanone.R Analysis/R/Rcache/agouti_scan.RData Analysis/R/Rcache/tufted_scan.RData
	cd R;R CMD BATCH --no-save coatcolor_lod.R

SuppFigs/figS20.eps: R/eqtl_counts_10.R Analysis/R/Rcache/neqtl_new_10.RData Analysis/R/Rcache/neqtl_old_10.RData
	cd R;R CMD BATCH --no-save eqtl_counts_10.R

### clean directories of intermediate files
clean:
	\rm -f *.aux *.bbl *.blg *.log *.bak *~ *.Rout */*~ */*.Rout */*.aux */*.log *.dvi
	cd ProcessedFiles;\rm -f *.tex *.Rnw *.aux *.log *.dvi *.bbl *.blg

### remove everything
cleanall: clean
	\rm -f *.pdf Figs/*.eps Figs/*.pdf SuppFigs/*.pdf

### proposed figure for G3 cover
Figs/cover.tiff: R/cover_fig.R
	cd R;R CMD BATCH --no-save cover_fig.R
