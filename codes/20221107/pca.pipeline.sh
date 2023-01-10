# ADD the reference panel while calculating PC.
study=${1}

mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/

echo "geno1: /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/input1.bed
snp1:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/input1.pedsnp
ind1:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/input1.pedind
geno2: /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input2.bed
snp2:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input2.pedsnp
ind2:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input2.pedind
genooutfilename:   /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.bed
snpoutfilename:    /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.pedsnp
indoutfilename:    /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.pedind
outputformat: PACKEDPED
allowdups:  YES " > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.par

/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/mergeit -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.par

echo "geno1: /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.bed
snp1:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.pedsnp
ind1:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge1.pedind
geno2: /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input3.bed
snp2:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input3.pedsnp
ind2:  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input3.pedind
genooutfilename:   /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.eigenstratgeno
snpoutfilename:    /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.snp
indoutfilename:    /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.ind
outputformat: EIGENSTRAT
allowdups:  YES " > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.par

/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/mergeit -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.par


mkdir /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/

module load perl

echo '#!/usr/bin/perl' > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/pca.perl

echo "\$ENV{'PATH'} = '/restricted/projectnb/adgc/_CANDIDATE.nb/2020.03.27_CTE/xdhan/new_gwas_data/QC.v2/PCA/EIG-7.2.1/bin/';

\$command = 'smartpca.perl ';
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.ind ';
\$command .= ' -k 20 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.plot ';
\$command .= ' -e /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.eval ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.log ';
\$command .= ' -m 0 ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'smarteigenstrat.perl '; 
\$command .= ' -i /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.eigenstratgeno ';
\$command .= ' -a /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.snp ';
\$command .= ' -b /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/merge/merge2.ind ';
\$command .= ' -p /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca ';
\$command .= ' -k 10 ';
\$command .= ' -o /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.chisq ';
\$command .= ' -l /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.chisq.log ';
print(\"\$command\n\");
system(\"\$command\");

\$command = 'gc.perl /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.chisq /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.chisq.GC';
print(\"\$command\n\");
system(\"\$command\");" >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/pca.perl

chmod a+x  /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/pca.perl
/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/pca.perl


awk '{print $2}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/input1.pedind > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.1

awk '{print $2}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input2.pedind >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.1

awk '{print $2}' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input3.pedind >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.1

sed '1,21d' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.2

paste -d ' ' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.1 /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.2 > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt

rm -rf /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/${study}/PCA/${study}.pca.txt.*
