# study sexes of swaps in expression arrays
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
load("Rcache/xist_and_y.RData")
necropsy <- read.csv("../Analysis/FinalData/necropsy.csv", as.is=TRUE)
sex <- sub("F", "Female", sub("M", "Male",necropsy$Sex))
names(sex) <- necropsy$MouseNum

arraySwaps <- vector("list", length(tissues))
names(arraySwaps) <- tissues

arraySwaps$adipose <- list(c("Mouse3583", "Mouse3584"),
                           c("Mouse3187", "Mouse3188", "Mouse3200"))
                               
arraySwaps$gastroc <- list(c("Mouse3655", "Mouse3659"))  # Female/Male

arraySwaps$hypo <- list(c("Mouse3179", "Mouse3188"),
                        c("Mouse3208", "Mouse3210"), # Male/Female
                        c("Mouse3347", "Mouse3348"),
                        c("Mouse3367", "Mouse3369"),
                        c("Mouse3381", "Mouse3382"),
                        c("Mouse3449", "Mouse3451"), # Female/Male
                        c("Mouse3452", "Mouse3454"),
                        c("Mouse3589", "Mouse3590"),
                        c("Mouse3592", "Mouse3594"))

arraySwaps$islet <- list(c("Mouse3598", "Mouse3599"),
                         c("Mouse3295", "Mouse3296"))

arraySwaps$kidney <- list(c("Mouse3523", "Mouse3510"))

arraySwaps$liver <- list(c("Mouse3142", "Mouse3143"),
                         c("Mouse3136", "Mouse3141"))

lapply(arraySwaps, lapply, function(a,b) b[a], sex)
