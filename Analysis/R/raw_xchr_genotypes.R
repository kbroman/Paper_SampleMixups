##############################
# load raw genotype data
#
# identify the mice with sex swaps
##############################
rawg <- read.cross("csv", "../OrigData", "genotypes4rqtl.csv", genotypes=0:2,
                   convertXdata=FALSE, alleles=c("B","R"))

# mouse with missing strain
rawg <- subset(rawg, ind=(rawg$pheno$Strain!=""))

# clean up sex and strain
rawg$pheno$Sex <- factor(as.character(rawg$pheno$Sex), levels=c("Female", "Male"))
rawg$pheno$Strain <- factor(as.character(rawg$pheno$Strain), levels=unique(as.character(rawg$pheno$Strain)))

# no genotype data; drop them
rawg$pheno[ntyped(rawg)==0,]
rawg <- subset(rawg, ind=(ntyped(rawg)>0))

# pull out X chr data
gx <- pull.geno(rawg, chr="X")
dim(gx)
strain <- rawg$pheno$Strain
sex <- rawg$pheno$Sex

# markers not segregating
ming <- apply(apply(gx[strain=="F2",], 2, function(a) table(factor(a, levels=1:3))), 2, min)
table(ming)
todrop <- names(ming[ming<120])
rawg <- drop.markers(rawg, todrop)
gx <- pull.geno(rawg, chr="X")
dim(gx)

# no. markers left on X: 20
nmar(rawg)["X"]

# B6: X chr
all(apply(gx[strain=="B6",], 2, function(a) length(unique(a[!is.na(a)]))) == 1)
b6xg <- apply(gx[strain=="B6",], 2, function(a) unique(a[!is.na(a)]))

# BTBR: X chr
all(apply(gx[strain=="BTBR",], 2, function(a) length(unique(a[!is.na(a)])) == 1))
btbrxg <- apply(gx[strain=="BTBR",], 2, function(a) unique(a[!is.na(a)]))

# B6 and BTBR really homozygous for different alleles?
all((b6xg == 1 & btbrxg==3) | (b6xg == 3 & btbrxg==1))

# swap genotype codes
toswap <- names(b6xg)[b6xg==3]
for(mar in toswap)
  rawg$geno[["X"]]$data[,mar] <- 4 - rawg$geno[["X"]]$data[,mar]
gx <- pull.geno(rawg, chr="X")

# genotypes in F1
f1male.gx <- gx[strain=="F1" & sex=="Male",]
f1female.gx <- gx[strain=="F1" & sex=="Female",]
all(is.na(f1male.gx) | f1male.gx==3)
all(is.na(f1female.gx) | f1female.gx==2)

# F2 males
f2male.gx <- gx[strain=="F2" & sex=="Male",]
sum(!is.na(f2male.gx) & f2male.gx==2)

# at least one het
wh <- apply(f2male.gx, 1, function(a) any(!is.na(a) & a==2))
sum(wh) # 18 
tab <- apply(f2male.gx[wh,], 1, function(a) table(factor(a[!is.na(a)], levels=1:3)))
tab[,order(tab[2,])] # 4 look like errors (one het); 14 look like females

# not male
malef2id <- as.character(rawg$pheno$MouseNum)[strain=="F2" & sex=="Male"]
notmale <- malef2id[apply(f2male.gx, 1, function(a) sum(!is.na(a) & a==2)>1)]

# not sure male
sum(apply(f2male.gx, 1, function(a) all(is.na(a) | a==3))) # 50

# F2 females
f2female.gx <- gx[strain=="F2" & sex=="Female",]
sum(!is.na(f2female.gx) & f2female.gx==1)

# look wrong
wh <- apply(f2female.gx, 1, function(a) any(!is.na(a) & a==1))
sum(wh)
tab <- apply(f2female.gx[wh,], 1, function(a) table(factor(a[!is.na(a)], levels=1:3)))
tab[,order(tab[1,])] # all 21 of these look like males

# not female
femalef2id <- as.character(rawg$pheno$MouseNum)[strain=="F2" & sex=="Female"]
notfemale <- femalef2id[apply(f2female.gx, 1, function(a) sum(!is.na(a) & a==1)>1)]

# not sure female
sum(apply(f2female.gx, 1, function(a) all(is.na(a) | a==3))) # 53

# save outcome
save(rawg, notmale, notfemale, file="../OrigData/rawg.RData")

