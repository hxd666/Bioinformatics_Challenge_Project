#!/bin/bash -l

module load plink/1.90b6.21
module load perl

# Step 1: generate bed files includes the Samples of interest
sed -i 's/K-0487/k-0487/g' UNITE_batch2.ID.tbx
sed -i 's/K-0554/k-0554/g' UNITE_batch2.ID.tbx

awk '{print $1"\t"$1}' UNITE_batch1.ID.tbx | sed '1d' > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/UNITE1.ID
awk '{print 0"\t"$1}' UNITE_batch2.ID.tbx | sed '1d' > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/UNITE2.ID

#sed -i 's/K-0487/k-0487/g' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/UNITE2.ID
#sed -i 's/K-0554/k-0554/g' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/UNITE2.ID


plink --bfile /restricted/projectnb/cte/_gwas.cte_lgen_org/_variants_without_GC.score/CRARY-313HTSsamples-1-27-17_FinalReport_maf1_geno5_hwe1e6_mind5 --keep /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/UNITE1.ID --make-bed --out /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1

awk '{print 0,$2,$3,$4,$5,$6}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam.tmp
mv /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam.tmp  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam

plink --bfile /restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC/IJME_Neurology_11-8-21.v2.QC --keep /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/UNITE2.ID --make-bed --out /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/unite2


# step 2: Convert extracted bed files to the required format for PCA
mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input
mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input

echo "genotypename:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.bed
snpname:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.bim
indivname:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam
outputformat:	EIGENSTRAT
genotypeoutname:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.eigenstratgeno
snpoutname:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.snp
indivoutname:	/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.ind" > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/convert.par

/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/convertf -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/convert.par

sed -i 's/0://' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.ind
sed -i 's/???/U/g' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.ind

echo "genotypename:     /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/unite2.bed
snpname:        /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/unite2.bim
indivname:      /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/unite2.fam
outputformat:   EIGENSTRAT
genotypeoutname:        /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.eigenstratgeno
snpoutname:     /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.snp
indivoutname:   /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.ind" > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/convert.par

/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/convertf -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/convert.par

sed -i 's/0://' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.ind
sed -i 's/???/U/g' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.ind

# Step 3: Compute the PCA
mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA
mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA

echo '#!/usr/bin/perl' > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/pca.perl

echo "\$ENV{'PATH'} = '/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/';

\$command = 'smartpca.perl ';
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.ind ';
\$command .= ' -k 20 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.plot ';
\$command .= ' -e /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.eval ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.log ';
\$command .= ' -m 0 ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'smarteigenstrat.perl '; 
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/input/unite1.snp ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca ';
\$command .= ' -k 10 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.chisq ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.chisq.log ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'gc.perl /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.chisq /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.chisq.GC';
print(\"\$command\n\");
system(\"\$command\");" >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/pca.perl

chmod a+x /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/pca.perl
/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/pca.perl


awk '{print $2}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/unite1.fam > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt.1
sed '1,21d' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt.2
paste -d ' ' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt.1 /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt.2 > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt
rm -rf /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt.*



echo '#!/usr/bin/perl' > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/pca.perl

echo "\$ENV{'PATH'} = '/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/';

\$command = 'smartpca.perl ';
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.ind ';
\$command .= ' -k 20 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.plot ';
\$command .= ' -e /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.eval ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.log ';
\$command .= ' -m 0 ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'smarteigenstrat.perl '; 
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/input/unite2.snp ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2bd_batch/PCA/unite2.pca ';
\$command .= ' -k 10 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.chisq ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.chisq.log ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'gc.perl /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.chisq /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.chisq.GC';
print(\"\$command\n\");
system(\"\$command\");" >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/pca.perl

chmod a+x /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/pca.perl
/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/pca.perl


awk '{print $2}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/unite2.fam > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite1.pca.txt.1
sed '1,21d' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt.2
paste -d ' ' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt.1 /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt.2 > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt
rm -rf /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt.*
