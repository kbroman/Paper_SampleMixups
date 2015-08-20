library(B6BTBR07a)
data(topstuff)

ts.mousenum <- sapply(strsplit(rownames(topstuff), "_"), function(a) a[[2]])
ts.plate <- as.character(as.numeric(sapply(strsplit(as.character(topstuff[,3]), "[_:]"), function(a) a[2])))
ts.well <- sapply(strsplit(as.character(topstuff[,3]), "[_:]"), function(a) a[3])
ts.strain <- sapply(strsplit(rownames(topstuff), "[_:]"), function(a) a[3])
ts.strain[rownames(topstuff) == "MouseKitControl_755356"] <- "CNTL"

names(ts.plate) <- ts.mousenum
names(ts.well) <- ts.mousenum

load("../OrigData/rawg.RData")
tmp <- strsplit(as.character(rawg$pheno$PlatePos), "[_:]")
rawg.plate <- as.numeric(sapply(tmp, function(a) a[[2]]))
rawg.well <- sapply(tmp, function(a) a[[3]])
names(rawg.plate) <- names(rawg.well) <- rawg$pheno$MouseNum

oplate <- ts.plate[names(rawg.plate)]
all(oplate == rawg.plate) # TRUE

owell <- ts.well[names(rawg.well)]
all(owell == rawg.well) # TRUE

# missing
ts.mousenum[(!(ts.mousenum %in% rawg$pheno$MouseNum))]
# "Mouse3377" "755356"    "Mouse3265" "Mouse3338"
