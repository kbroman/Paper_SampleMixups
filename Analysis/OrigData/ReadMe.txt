Original data (before correction of sample mix-ups)

rawg.RData: even more raw genotypes, plus with info on sex-swapped mice
  (created from genotypes4rqtl.csv with ../R/raw_xchr_genotypes.R)


clean_cross.RData: clean genotype data; contains sample duplicates; no corrections of sample mix-ups
                   within cross: Cox et al map
                   newmap = map estimated from data
   f2g: f2 genotype data
        phenotype Sex has corrections from genotypes
        SexID is the sex from the identifier

F2.mlratio.[tissue].RData  Original expression data (mlratios)
  each contains [tissue].mlratio, which is (# ind'ls) x (# genes)
     with # genes = 40572

