cd /sc/orga/scratch/webste01/TRF_analysis/data 
for f in *fa;
do
trf $f 2 7 7 80 10 50 20 -f -d -h -ngs > ${f}_out
done
mv *_out ../results/trf_results

