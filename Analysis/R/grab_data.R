# grab the data from the Mouse Phenome Database
# and put it in the form needed for this repository

library(data.table)

mpd_url <- "http://phenome.jax.org/grpdoc_qtla/attie_2015"
raw_zipfile <- "Attie_2015_eqtl_raw.zip"
clean_zipfile <- "Attie_2015_eqtl_clean.zip"

# temporary directory
dir <- tempdir()
raw_path <- file.path(dir, raw_zipfile)
clean_path <- file.path(dir, clean_zipfile)

# download raw data
cat(" -Downloading raw zip file (nearly 1 GB)\n")
utils::download.file(paste0(mpd_url, "/", raw_zipfile), raw_path)

# unzip
cat(" -Unzipping raw data\n")
unzipped_files <- utils::unzip(raw_path, exdir=dir)

# function to read csv file
myfread <-
    function(file)
{
    x <- data.table::fread(file, data.table=FALSE, header=TRUE)
    rownames(x) <- x[,1]
    x[,-1]
}

# reorg microarray data
cat(" -Reading and saving mlratios\n")
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
for(tissue in tissues) {
    cat(" ---Reading and saving", tissue, "\n")
    obj <- paste0(tissue, ".mlratio")
    assign(obj,
           t(myfread(file.path(dir, "Raw", paste0(tissue, "_mlratio_raw.csv")))))
    save(list=obj, file=file.path("..", "OrigData", paste0("F2.mlratio.", tissue, ".RData")))
}

# read very raw genotypes and make some slight changes
suppressWarnings(rawg <- read.cross("csv", file.path(dir, "Raw"), "genotypes_veryraw.csv",
                                    genotypes=0:2, convertXdata=FALSE))
# reorder chromosomes
rawg$geno <- rawg$geno[c(1:19, "X", "un")]
# drop uninformative markers from X
suppressWarnings(genocount <- apply(geno.table(rawg, chr="X")[,c("AA", "AB")], 1, min))
rawg <- drop.markers(rawg, names(genocount)[genocount < 5])
# drop some individuals
rawg <- rawg[,-grep("Kit", as.character(rawg$pheno$longid))]
rawg <- rawg[,nmissing(rawg) < 1000]
save(rawg, file=file.path("..", "OrigData", "rawg.RData"))

# read raw data; reorder chr numbers
suppressWarnings(f2g <- read.cross("csv", file.path(dir, "Raw"), "genotypes_raw.csv",
                                   genotypes=c("BB", "BR", "RR"), alleles=c("B", "R")))
# reorder chromosomes
f2g$geno <- f2g$geno[c(1:19, "X", "un")]
save(f2g, file=file.path("..", "OrigData", "clean_cross.RData"))

cat(" -Downloading clean zip file (nearly 1 GB)\n")
utils::download.file(paste0(mpd_url, "/", clean_zipfile), clean_path)

# unzip
cat(" -Unzipping clean data\n")
unzipped_files <- utils::unzip(clean_path, exdir=dir)

# copy annotation file
annot <- fread(file.path(dir, "microarray_annot.csv"), data.table=FALSE)
save(annot, file=file.path("..", "OrigData", "annot.amit_rev.RData"))
