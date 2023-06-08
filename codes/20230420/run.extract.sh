#!/bin/bash -l

zcat /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/metal/CA1/UNITE_CA1.assoc.1.fltdf.txt.tbx.gz | head -n 1 > rs6910507.assoc.1.fltdf.txt.tbx

while read -r endopheno; do
	echo "${endopheno}"
	zgrep '6-29995149' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/metal/${endopheno}/UNITE_${endopheno}.assoc.1.fltdf.txt.tbx.gz >> rs6910507.assoc.1.fltdf.txt.tbx
done < endopheno.list

echo "CTE"
zgrep '6-29995149' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/simple_model/metal/CTE/UNITE_CTE.assoc.1.fltdf.txt.tbx.gz >> rs6910507.assoc.1.fltdf.txt.tbx

echo "CTEstage"
zgrep '6-29995149' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/simple_model/metal/CTEstage/UNITE_CTEstage.assoc.1.fltdf.txt.tbx.gz >> rs6910507.assoc.1.fltdf.txt.tbx

echo "AT8PositiveTissue.log"
zgrep '6-29995149' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/linear_model/metal/AT8PositiveTissue.log/UNITE_AT8PositiveTissue.log.assoc.1.fltdf.txt.tbx.gz >> rs6910507.assoc.1.fltdf.txt.tbx

paste -d "\t" pheno.list rs6910507.assoc.1.fltdf.txt.tbx > rs6910507.assoc.1.fltdf.txt.final.tbx
