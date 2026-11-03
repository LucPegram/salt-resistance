# Pre-processing of microarray data from salt-treated Col-0 A. thaliana.
# 2026 Nov. 03

library(limma)

nm_f1 <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/Hydroponics array/Text files/",
                "B1 50N root - US83103526_256398510023_S01_GE1_1105_Oct12_1_3.txt")

df1 <- read.delim(nm_f1, sep = '\t', header = F)

