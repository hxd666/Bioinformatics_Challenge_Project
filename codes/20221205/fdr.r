# recalculate the FDR
predictors <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221204/predictors_list',
                         as.is = T)
all.models.fn <- sapply(predictors$V1, function(x){paste0('results/Model.All.',x,'.assoc.txt')})
cases.models.fn <- sapply(predictors$V1, function(x){paste0('results/Model.Cases.',x,'.assoc.txt')})
age.models.fn <- sapply(predictors$V1, function(x){paste0('results/Model.Age_Stratified.',x,'.assoc.txt')})

all.models.fn.filter <- c()
for (.models in all.models.fn){
  if (file.exists(.models)){
    all.models.fn.filter <- c(all.models.fn.filter, .models)
  }
}

cases.models.fn.filter <- c()
for (.models in cases.models.fn){
  if (file.exists(.models)){
    cases.models.fn.filter <- c(cases.models.fn.filter, .models)
  }
}

age.models.fn.filter <- c()
for (.models in age.models.fn){
  if (file.exists(.models)){
    age.models.fn.filter <- c(age.models.fn.filter, .models)
  }
}

