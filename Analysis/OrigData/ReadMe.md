## Original data (before correction of sample mix-ups)

To be populated by the script [`grab_data.R`](../R/grab_data.R)

`F2.mlratio.[tissue].RData`  Original expression data (mlratios)
  each contains `[tissue].mlratio`, which is (# genes) x (# ind'ls)
     with # genes = 40572

clean_cross.RData: clean genotype data; contains sample duplicates;
no corrections of sample mix-ups within cross: Cox et al map
(`f2g` is the F2 genotype data; phenotype `Sex` has corrections from
genotypes; `SexID` is the sex from the individual identifiers)

`rawg.RData`: has even more raw genotypes, including info on sex-swapped mice
