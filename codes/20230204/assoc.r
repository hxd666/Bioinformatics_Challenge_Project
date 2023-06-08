# association test for 1) FAQ, 2) CDS, 3) BRIEF BRI
setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230204')

load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/imputatoin.RData')

new_pheno <- haven::read_sav('/restricted/projectnb/cte/Challenge_Project_2022/January 2023 AT8 Update/APOEAnalysisscales.sav')

# 1) FAQ
pheno$faqtot <- new_pheno$faqtot[match(pheno$subjid, new_pheno$subjid)]
sum(is.na(pheno$faqtot)) # 59 is missing

# 2) CDS
pheno$CDStot <- new_pheno$CDStot[match(pheno$subjid, new_pheno$subjid)]
sum(is.na(pheno$CDStot)) # 114 is missing

# 3) BRIEF BRI
pheno$tbri   <- new_pheno$tbri[match(pheno$subjid, new_pheno$subjid)]
sum(is.na(pheno$tbri)) # 126 is missing

model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)

# Only European
subset = pheno$race %in% c(1) & !is.na(pheno$H1B1G1)
for (.pheno in c('faqtot', 'CDStot', 'tbri')){
  form1 <- as.formula(paste(.pheno, '~ H1B1G1 + agedeath +', paste(sapply(1:10, function(x){paste0('PC',x)}), collapse = " + ")))
  fit   <- glm(form1, data = pheno, na.action = na.omit, family = gaussian(), subset = subset)
  result <- summary(fit)$coefficients
  assoc_subject <- !is.na(pheno[,.pheno]) & subset
  n <- sum(assoc_subject)
  effect <- round(result['H1B1G1', 1], 3)
  se <- round(result['H1B1G1', 2], 3)
  p  <- result['H1B1G1', 4]
  freq <- round(sum(pheno$H1B1G1[assoc_subject])/(2*n), 3)
  output.row <- c('Euro', '1', .pheno, 'H1B1G1', freq, n, effect, se, p)
  output.row <- unname(output.row)
  model.all.data <- rbind(model.all.data, output.row)
}

model.all.data <- model.all.data[2:nrow(model.all.data), ]
outputname <- 'results/Model1.Euro.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
