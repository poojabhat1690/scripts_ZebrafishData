ml samtools

touch  /clustertmp/pooja/STARmapping/allMapped/factors_normalize.txt
touch /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/scripts_ZebrafishData/ucscBrowser/factors_normalize.txt

cd /clustertmp/pooja/STARmapping/allMapped/
var=$(find *.bam)

for f in $var
do

mapped=`samtools flagstat /clustertmp/pooja/STARmapping/allMapped/"$f" | cut -f1 -d + - | sed -n '5p' -`
million=1000000
fac=`echo $million / $mapped | bc -l` 



echo $f $fac >>/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/scripts_ZebrafishData/ucscBrowser/factors_normalize.txt

done

