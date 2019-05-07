# grab the data from the Mouse Phenome Database
# and put it in the form needed for this repository

library(data.table)
library(qtl)

mpd_url <- "https://phenomedoc.jax.org/QTL_Archive/attie_2015/"
raw_zipfile <- "Attie_2015_eqtl_raw.zip"
clean_zipfile <- "Attie_2015_eqtl_clean.zip"
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")

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
for(tissue in tissues) {
    cat(" ---Reading and saving", tissue, "\n")
    obj <- paste0(tissue, ".mlratio")
    assign(obj,
           t(myfread(file.path(dir, "Raw", paste0(tissue, "_mlratio_raw.csv")))))
    save(list=obj, file=file.path("..", "OrigData", paste0("F2.mlratio.", tissue, ".RData")))
}

# drop some arrays and save in mlratios.RData
omit <- read.csv("omitted_samples.csv", as.is=TRUE)
omit <- omit[omit$Tissue != "DNA",]
omit <- omit[-grep("duplicates", omit$Note),]
for(tissue in tissues) {
    obj <- paste0(tissue, ".mlratio")
    assign(obj, get(obj)[is.na(match(rownames(get(obj)), omit$MouseNum[omit$Tissue==tissue])),])
}
save(list=paste0(tissues, ".mlratio"), file.path("Rcache", "mlratios.RData"))

# read very raw genotypes and make some slight changes
suppressWarnings(rawg <- read.cross("csv", file.path(dir, "Raw"), "genotypes_veryraw.csv",
                                    genotypes=0:2, convertXdata=FALSE))
# save plate info
longid <- as.character(rawg$pheno$longid)
platepos <- as.character(rawg$pheno$PlatePos)
mousenum <- sapply(strsplit(longid, "_"), function(a) a[[2]])
plate <- as.character(as.numeric(sapply(strsplit(platepos, "[_:]"), function(a) a[2])))
well <- sapply(strsplit(platepos, "[_:]"), function(a) a[3])
strain <- sapply(strsplit(longid, "[_:]"), function(a) a[3])
strain[longid == "MouseKitControl_755356"] <- "CNTL"
plateinfo <- data.frame(longid=longid,
                        platepos=platepos,
                        mousenum=mousenum,
                        plate=plate,
                        well=well,
                        strain=strain,
                        stringsAsFactors=FALSE)
write.table(plateinfo, file=file.path("..", "OrigData", "plateinfo.csv"),
          row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)

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

# reorg microarray data
cat(" -Reading and saving mlratios\n")
for(tissue in tissues) {
    cat(" ---Reading and saving", tissue, "\n")
    obj <- paste0(tissue, ".mlratio")
    assign(obj,
           myfread(file.path(dir, "Clean", paste0(tissue, "_mlratio_clean.csv"))))
    save(list=obj, file=file.path("..", "FinalData", paste0(tissue, "_mlratio_final.RData")))
}

# copy annotation file
annot <- fread(file.path(dir, "Clean", "microarray_annot.csv"), data.table=FALSE)
save(annot, file=file.path("..", "OrigData", "annot.amit_rev.RData"))

# clean genotype data
f2g <- read.cross("csv", file.path(dir, "Clean"), "genotypes_clean.csv",
                  genotypes=c("BB", "BR", "RR"), alleles=c("B", "R"))
f2g$geno <- x$geno[c(1:19, "X")]

# genetic and physical maps
gmap <- myfread(file.path(dir, "Clean", "markers_genetic_map.csv"))
pmap <- myfread(file.path(dir, "Clean", "markers_physical_map.csv"))
mnames <- split(rownames(pmap), factor(pmap$chr, levels=c(1:19,"X")))
gmap <- split(gmap$pos, factor(gmap$chr, levels=c(1:19,"X")))
pmap <- split(pmap$pos, factor(pmap$chr, levels=c(1:19,"X")))
for(i in seq(along=gmap)) {
    names(gmap[[i]]) <- mnames[[i]]
    names(pmap[[i]]) <- mnames[[i]]
}

# save genotypes + physical map
save(f2g, pmap, file=file.path("..", "FinalData", "aligned_geno_with_pmap.RData"))
