# The Rscript to convert the hg38 to hg37 for bim file
rm(list = ls())
study <- commandArgs(T)[1]

# Read bim file
bim.file.dir <- "/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/XXXX_batch/_mapt_region/bed_files/MAPT.topmed.bim"
bim.file.dir <- gsub("XXXX", study, bim.file.dir)
bim.file <- read.table(bim.file.dir, as.is = T, header = F)

# Read legend file
rsid <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/ref.paper.impute.legend',
                   header = T, as.is = T)

# Read locid hg38 file
locid.38 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/candidate.locid.hg38.pos.txt', 
                       header = F, as.is = T)

rsid.38 <- rsid
rsid.38$position <- locid.38[,1]

bim.file$V2 <- rsid.38$rsID[match(bim.file$V4, rsid.38$position)]
bim.file$V4 <- rsid$position[match(bim.file$V2, rsid$rsID)]

write.table(bim.file, bim.file.dir, col.names = F, row.names = F,
            quote = F, sep = "\t")
