require(analyzeMC)
require(backShift)
require(reshape2)
require(Hmisc)
###############
### ADJUST ####
###############
resultPath <- "C:/Users/Daniel Tong/Desktop/school stuff/DSC/DSC 180A/BM/_test/"
thresConv <- 75
saveResultsPath <- "C:/Users/Daniel Tong/Desktop/school stuff/DSC/DSC 180A/BM/"
takeout <- "cd8+_res.rds"
folderName <- "effects"
woInhibitors <- c("")

allInhibitors <-  c("Akti", "BTKi", "Crassin", "Dasatinib", "GDC-0941",
                    "Go69", "H89","IKKi","Imatinib","Jak1i","Jak2i","Jak3i",
                    "Lcki","Lestaurtinib","PP2","Rapamycin","Ruxolitinib",
                    "SB202","Sorafenib","SP6","Staurosporine","Streptonigrin",
                    "Sunitinib","Syki","Tofacitinib","U0126","VX680")
allActs <-  c("BCR","GCSF","GMCSF","IFNa","IFNg","IL12","IL2","IL3","LPS",
              "PMA-IONO","pVO4","Ref")

files <- aggregateResults(resultPath, thresConv, takeout, woInhibitors,
                          allInhibitors, allActs, saveResultsPath, folderName)
