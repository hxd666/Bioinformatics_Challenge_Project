rm( list=ls() )
require(MASS)
require(mice)
library(AER)
options( width=180 )

#########################################################################################
####	The things to change
#########################################################################################
#------------------------------------------------------
#---	input parameters
#------------------------------------------------------
t.args <- commandArgs(T)
risk   <- t.args[1] # AT8PositiveTissue.log
study  <- t.args[2] # batch1 or batch2
chr    <- t.args[3] # number of chromosome

risk   <- 'AT8PositiveTissue.log'
#study  <- 'batch2'
#chr    <- '2'
dir.name <- chr
#------------------------------------------------------
#---	Set up a temporary working directory
#------------------------------------------------------
tmpdir<-    try( Sys.getenv("TMPDIR"), silent=TRUE )
if ("try-error" %in% class(tmpdir) | tmpdir == "")      {
  cat("...working not in queue", fill=T)
  tmpdir<-	getwd()
}       else    {
  if (nchar(tmpdir) < 5) stop("TMPDIR variable is not set")
}
cat( "...temporary directory: ", tmpdir, fill=T )

#------------------------------------------------------
#---	set up the output directory
#------------------------------------------------------
dir.out <-		"/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/linear_model/"
dir.out <- paste(dir.out, study, "/", sep = "")
dir.out <- paste(dir.out, risk, "/", sep = "")

#------------------------------------------------------
#---	set up phenotype and covariate files
#------------------------------------------------------
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/imputatoin.RData')
data.pheno <-	lapply(1:10, function(x){complete(imp,x)})

# filter by race
data.pheno <- lapply(data.pheno, function(x){x[x$race %in% c(1),]})
if (study == "batch1"){
  data.pheno <- lapply(data.pheno, function(x){x[x$Batch == 'I',]})
} else if (study == "batch2"){
  data.pheno <- lapply(data.pheno, function(x){x[x$Batch == 'II',]})
}
#------------------------------------------------------
#--- set up the covariates 
#------------------------------------------------------
age.col.name<-	"agedeath"
ph.col.name<-			risk
pcs.col.name<-	unname(sapply(1:10, function(x){paste0('PC',x)}))
data.pheno <- lapply(data.pheno, function(x){x[, c("subjid", "pin", age.col.name, ph.col.name, pcs.col.name)]})


#------------------------------------------------------
#---	set up the output file name
#------------------------------------------------------
ph.col.name.4outfile <- risk
study.outname<-	paste( "UNITE","_", study, "_", ph.col.name.4outfile, sep="" )
fn.out<-  paste( tmpdir, "AAA_chrXXX.assoc", sep="/" )

fn.out<-  gsub( "AAA", study.outname, fn.out )
fn.out<-	gsub( "XXX", chr, fn.out )

fn.out.gz<-	paste( dir.out, "chrXXX/AAA_chrXXX.assoc.gz", sep="" )

fn.out.gz<-	gsub( "AAA", study.outname, fn.out.gz )
fn.out.gz<-	gsub( "XXX", dir.name, fn.out.gz )

#------------------------------------------------------
#--- set up vcf sample
#------------------------------------------------------
if (study == 'batch1'){
  sample.dir <- "/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_TOPMed_imputation/chrXXXX.dose.sample"
  sample.dir <- gsub("XXXX",chr,sample.dir)
  data.sample <- read.table(sample.dir, header = F, as.is = T)
  data.sample$V1 <- sapply(data.sample$V1, function(x){strsplit(x, split = '_')[[1]][1]})
  data.sample <- data.sample$V1[match(data.pheno[[1]]$pin, data.sample$V1)]
  data.sample <- sapply(data.sample, function(x){paste(x,x,sep='_')})
  data.sample <- data.frame('pin'  = data.sample, stringsAsFactors = F)
  data.sample.dir <- paste0(dir.out, 'chr', chr, '/sample.txt')
  if (!file.exists(data.sample.dir)){
    write.table(data.sample, data.sample.dir, col.names = F, row.names = F, quote = F)
  }
} else if (study == 'batch2'){
  sample.dir <- "/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/_TOPMed_imputation/chrXXXX.dose.sample"
  sample.dir <- gsub("XXXX",chr,sample.dir)
  data.sample <- read.table(sample.dir, header = F, as.is = T)
  data.sample$V1[data.sample$V1=='k-0487'] <- "K-0487"
  data.sample$V1[data.sample$V1=='k-0554'] <- "K-0554"
  data.sample <- data.sample$V1[match(data.pheno[[1]]$subjid, data.sample$V1)]
  data.sample <- data.frame('subjid'  = data.sample, stringsAsFactors = F)
  data.sample$subjid[data.sample$subjid=='K-0487'] <- "k-0487"
  data.sample$subjid[data.sample$subjid=='K-0554'] <- "k-0554"
  data.sample.dir <- paste0(dir.out, 'chr', chr, '/sample.txt')
  if (!file.exists(data.sample.dir)){
    write.table(data.sample, data.sample.dir, col.names = F, row.names = F, quote = F)
  }
}

#------------------------------------------------------
#---	Set up dosage files
#------------------------------------------------------
if (study == 'batch1'){
  fn.dosage1.org <-	"/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_TOPMed_imputation/chrXXXX.dose.vcf.gz"
} else if (study == 'batch2'){
  fn.dosage1.org <-	"/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/_TOPMed_imputation/chrXXXX.dose.vcf.gz"
}
fn.dosage1 <-	gsub( "XXXX", chr, fn.dosage1.org )
if( !file.exists( fn.dosage1 ) ){
    cat("[Error] No dosage file\n")
    cat( fn.dosage1, fill=T )
    }


#-------------------------------------------------------------------
#---	function for one association test between a SNP and an outcome
#-------------------------------------------------------------------

comp.func.glm_dx.int<-  function(df,i.freq,i.marker,i.var,i.data.dos1,i.num.samp,outFile)  {
  #print(i.marker)
  i.chr <- gsub('chr','',df['CHR'])
  i.bp <- df['BP']
  i.locID <- paste(i.chr, i.bp, sep = '-')
  #print(i.data.dos1[,i.marker])
  if (i.freq<=0.01 | i.freq>=0.99 | is.na(i.freq) | i.var < 0.4)  {
    i.glm.out<-	c( i.freq, i.var, i.num.samp, rep(NA, 4) )
  }	else    {
    tryCatch({
      #print(i.freq)
      #------------------------------------------------------
      #---	study model
      #------------------------------------------------------
      all.covar.names.4glm<-	c( age.col.name, pcs.col.name )
      i.form <- paste(risk, i.marker, sep = " ~ ")
      i.form <- paste(c(i.form, all.covar.names.4glm), collapse = ' + ')
      i.form <- as.formula( i.form )
      
      lme1 <- lapply(i.data.dos1, function(x){glm(i.form, x, family = gaussian(), na.action = na.omit)})
      i.glm.out.table <- summary(pool(lme1))
      i.glm.out<-	c(i.freq, i.var, i.num.samp)
      
      
      #---	snp term
      if ( sum( i.glm.out.table$term %in% i.marker ) == 1 )    {
        i.glm.out.snp<-	unlist( i.glm.out.table[i.glm.out.table$term %in% i.marker, c(2,3,4,6)])
      }	else	{
        i.glm.out.snp<-	rep(NA, 4)
      }

      #---  the covariance between interaction term and snp main effect term
      # i.glm.out.cov.table <- vcov_gee(lme1)
      
      #---	combine those three summary stat
      i.glm.out <-	c( i.glm.out, i.glm.out.snp)
      #print(i.marker)
      
    }, error = function(e){
      i.glm.out <<- c( i.freq, i.var, i.num.samp, rep(NA, 4) )
    })
  } 
  #---	interaction term
  i.outline<-	as.vector( unlist( df[1:5]) )
  i.outline[1] <- i.chr
  names(i.outline)<-	NULL
  names(i.glm.out)<-	NULL
  #i.outline <- paste(i.outline,collapse = '\t')
  #i.glm.out <- paste(i.glm.out,collapse = '\t')
  i.output <- paste(c(i.outline,i.glm.out,i.locID),collapse = '\t')
  i.output <- paste(i.output,'\n',sep = '')
  #cat(i.outline, file=outFile)
  #cat("\t", file=outFile)
  #cat(i.glm.out, file=outFile, fill=256)
  return(c(i.output))
  # cat( i.outline, i.glm.out, fill=256 )
}




#-------------------------------------------------------------------
#---	open the output object to write down
#-------------------------------------------------------------------
files.conn<-  file(fn.out,"w")

cat("CHR\tBP\tSNP\tA1\tA2\tFREQ\tINFO\tN\tBETA.snp\tSE.snp\tstat.snp\tP.snp\tLocID\n", file=files.conn)


#-------------------------------------------------------------------
#---	prepare phenotype and covaritates
#-------------------------------------------------------------------
cat("_______________________________________________________________________", fill=T)
cat( study, fill=T) 
cat( ph.col.name, fill=T) 
cat( chr, fill=T)
cat( pcs.col.name, fill=T)
cat( fn.out, fill=T)
cat("_______________________________________________________________________", fill=T)


buffer.nlines<- 1000     ### the number of lines to read dosage lines at a time
cmd.tabix<-	paste("bcftools query -f '%CHROM\\t%POS\\t%ID\\t%ALT{0}\\t%REF[\\t%DS]\\n' -S", data.sample.dir, fn.dosage1, sep = " ")
#cmd.tabix<- gsub("XXXX", fn.dosage1, cmd.tabix)
#cmd.tabix<- gsub("YYYY", chr, cmd.tabix)
cat( cmd.tabix, fill=T )
#cmd.tabix <- 'less /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/linear_model/pipe.txt'

#-------------------------------------------------------------------
#--- Run the association tests by 100 snps at once
#-------------------------------------------------------------------
pipe.tabix<-    pipe( cmd.tabix, "r" )
cat("...computing regression", fill=T)
k<-	0
while( length( stream<- readLines( pipe.tabix, n=buffer.nlines) ) > 0 )    {
  #---	show how much the job progress into the screen
  k<-     k+1;	cat('.'); cat(k)
  
  #---	grab the dosage data from the memory
  stream.df<-     data.frame(do.call(rbind,strsplit(stream,"\t")),stringsAsFactors=F)
  
  #---	extract the snp info from the stream.df
  i.data.info<-	stream.df[, 1:5]; colnames( i.data.info )<-     c('CHR','BP','SNP','A1','A2')
  
  #---	make unique marker name 
  i.data.info$order<-	1:nrow( i.data.info )
  i.data.info$marker<-	paste( "X", i.data.info$CHR, i.data.info$BP, i.data.info$order, sep="_" )
  
  #---	extract the dosage info from the stream.df
  #---	transpose the dosage data (Now, it will have subjects in rows and snps in columns.)
  i.data.tdos1<-	t( stream.df[, 6:ncol(stream.df)] )
  
  #i.data.tdos1<- i.data.tdos1[match(data.pheno$vcfID, data.sample$vcfID), ]
  
  if (is.null(dim(i.data.tdos1))){
    i.data.tdos1 <- matrix(i.data.tdos1, ncol = 1)
  }
  #---	convert the string into the numberic values
  i.data.tdos1<-	apply( i.data.tdos1, 2, as.numeric )
  
  #---	convert the data into dataframe 
  i.data.tdos1<-	as.data.frame( i.data.tdos1, as.is=T, stringsAsFactors=F )
  
  #---	the number of subjects in the sample file should match the number of subjects in the dosage data.
  if( nrow( data.pheno[[1]] ) != nrow( i.data.tdos1 ) )       {
    cat("the number of samples written in the sample file differs from the number of sample in the dosage file", fill=T)
    quit("yes")
  }
  
  all.covar.names<-	c( age.col.name, pcs.col.name )
  
  colnames( i.data.tdos1 )<-	i.data.info$marker
  if (study == 'batch1'){
    i.data.tdos1$ID<-			data.pheno[[1]]$pin
  } else if (study == 'batch2'){
    i.data.tdos1$ID<-     data.pheno[[1]]$subjid
  }
  
  #---	merge the dosage data and the phenotype/covariate data
  i.data.tdos1.0<-	i.data.tdos1
  i.data.tdos1<-	lapply(1:10, function(x){cbind(data.pheno[[x]], i.data.tdos1)})
  
  #---	conduct association tests
  # =====================================================
  #------------------ PARALLEL --------------------------
  # =====================================================
  for (.i in 1:nrow(i.data.info)){
    df <- i.data.info[.i,]
    #df
    #print('one')
    df <- unlist(df)
    i.marker<-	df['marker']
    i.bp<-		df['BP']
    i.chr <- gsub('chr','',df['CHR'])
    i.data.dos1<-	lapply(i.data.tdos1, function(x){x[, c('ID',ph.col.name, all.covar.names, i.marker) ]})
    i.data.dos1<-	lapply(i.data.dos1, function(x){na.omit(x)})
    i.num.samp<-	nrow( i.data.dos1[[1]] )
    
    i.freq<-		mean( i.data.dos1[[1]][, c(i.marker)], na.rm=T ) / 2
    
    if (i.freq==1 | i.freq==0 | is.na(i.freq))	{
      i.var<-	NA
    }	else	{
      i.var<-	var ( i.data.dos1[[1]][, c(i.marker)], na.rm=T ) / (2*i.freq*(1-i.freq))
    }
    df <- comp.func.glm_dx.int(df = df,i.freq = i.freq,i.marker = i.marker,i.var = i.var,i.data.dos1 = i.data.dos1,i.num.samp = i.num.samp,outFile = files.conn)
    cat(df, file=files.conn)
  }  
}
cat("_______________________________________________________________________", fill=T)

#-------------------------------------------------------------------
#---	close the access of dosage1 file and the output object
#-------------------------------------------------------------------
close( pipe.tabix )
close(files.conn)


#-------------------------------------------------------------------
#---	archive output files and place output to the project space
#-------------------------------------------------------------------
cmd.gzip<-    paste("gzip -c9f ", fn.out)
cmd.gzip<-    paste(cmd.gzip, " > ", fn.out.gz, sep="")
cat(cmd.gzip, fill=T)
system(cmd.gzip)


#-------------------------------------------------------------------
#---	make archive files executable ???
#-------------------------------------------------------------------
cmd.chmod<-   paste("chmod 775 ", fn.out.gz, sep="" )
cat(cmd.chmod, fill=T)
system(cmd.chmod)


#-------------------------------------------------------------------
#---	delete remaining files from scratch
#-------------------------------------------------------------------
cmd.rm <-     paste("rm -rf ", fn.out)
cat(cmd.rm, fill=T)
system(cmd.rm)

#cmd.rm <-     paste("rm -rf", sink.tmp)
cat(cmd.rm, fill=T)
system(cmd.rm)

cat ("END_without_Error", fill=T)

