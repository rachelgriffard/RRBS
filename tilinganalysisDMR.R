##Rachel Griffard
##12112022

###############################Find DMR - methylKit###################################
###Install packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("methylKit")
library(methylKit)

### Create methylRaw list object
# replace "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/" with WD


## BBN v BBN+5ppm
file.list5=list( 
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1092_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1093_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1094_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1095_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1096_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1098_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1099_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1101_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1102_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1103_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1104_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1105_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1106_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1077_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1078_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1079_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1080_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1081_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1083_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1085_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1086_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1087_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1088_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1089_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1090_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1091_1_val_1_bismark_bt2_pe.bismark.cov"
)

myobj5ppm=methRead(file.list5,
               sample.id=list("test1","test2","test3","test4","test5","test6","test7","test8",
                              "test9","test10","test11","test12","test13","ctrl1","ctrl2",
                              "ctrl3","ctrl4","ctrl5","ctrl6","ctrl7","ctrl8","ctrl9",
                              "ctrl10","ctrl11","ctrl12","ctrl13"),
               assembly="Mm10",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
               mincov=10,
               header=FALSE,
               pipeline="bismarkCoverage"
)

## BBN v BBN+10ppm
file.list10=list( 
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1107_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1108_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1109_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1110_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1111_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1112_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1113_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1114_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1115_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1116_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1117_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1118_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1119_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1120_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1121_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1077_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1078_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1079_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1080_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1081_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1083_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1085_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1086_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1087_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1088_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1089_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1090_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1091_1_val_1_bismark_bt2_pe.bismark.cov"
)

myobj10ppm=methRead(file.list10,
               sample.id=list("test1","test2","test3","test4","test5","test6","test7","test8",
                              "test9","test10","test11","test12","test13","test14","test15","ctrl1","ctrl2",
                              "ctrl3","ctrl4","ctrl5","ctrl6","ctrl7","ctrl8","ctrl9",
                              "ctrl10","ctrl11","ctrl12","ctrl13"),
               assembly="Mm10",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
               mincov=10,
               header=FALSE,
               pipeline="bismarkCoverage"
)

## BBN v BBN+50ppm
file.list50=list( 
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1122_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1123_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1124_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1125_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1128_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1129_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1130_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1131_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1132_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1133_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1134_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1135_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1136_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1077_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1078_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1079_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1080_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1081_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1083_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1085_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1086_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1087_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1088_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1089_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1090_1_val_1_bismark_bt2_pe.bismark.cov",
  "/panfs/pfs.local/home/r816g589/work/RRBS/processed_20221212/D1091_1_val_1_bismark_bt2_pe.bismark.cov"
)

myobj50ppm=methRead(file.list50,
               sample.id=list("test1","test2","test3","test4","test5","test6","test7","test8",
                              "test9","test10","test11","test12","test13","ctrl1","ctrl2",
                              "ctrl3","ctrl4","ctrl5","ctrl6","ctrl7","ctrl8","ctrl9",
                              "ctrl10","ctrl11","ctrl12","ctrl13"),
               assembly="Mm10",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
               mincov=10,
               pipeline="bismarkCoverage",
               header=FALSE
               )

#####Tiling analysis - 1000bp (Identification of DMR)
## BBN v BBN+5ppm
tiles5 <- tileMethylCounts(myobj5ppm,win.size=1000,step.size=1000)
head(tiles5[[1]],3)
## BBN v BBN+10ppm
tiles10 <- tileMethylCounts(myobj10ppm,win.size=1000,step.size=1000)
head(tiles10[[1]],3)
## BBN v BBN+50ppm
tiles50 <- tileMethylCounts(myobj50ppm,win.size=1000,step.size=1000)
head(tiles50[[1]],3)

### Unite
## BBN v BBN+5ppm
meth5 <- unite(tiles5, destrand=F)
clusterSamples(meth5, dist="correlation", method="ward", plot=TRUE)
meth10 <- unite(tiles10, destrand=F)
clusterSamples(meth10, dist="correlation", method="ward", plot=TRUE)
meth50 <- unite(tiles50, destrand=F)
clusterSamples(meth50, dist="correlation", method="ward", plot=TRUE)

#Find differentially methylated regions (DMRs)
diff5 <- calculateDiffMeth(meth5)
diff_all5 <- getMethylDiff(diff5, difference=25, qvalue=0.01)
write.csv(diff_all5,"/panfs/pfs.local/home/r816g589/work/RRBS/DiffAll_5ppm.csv")

diff10 <- calculateDiffMeth(meth10)
diff_all10 <- getMethylDiff(diff10, difference=25, qvalue=0.01)
write.csv(diff_all10,"/panfs/pfs.local/home/r816g589/work/RRBS/DiffAll_10ppm.csv")

diff50 <- calculateDiffMeth(meth50)
diff_all50 <- getMethylDiff(diff50, difference=25, qvalue=0.01)
write.csv(diff_all50,"/panfs/pfs.local/home/r816g589/work/RRBS/DiffAll_50ppm.csv")