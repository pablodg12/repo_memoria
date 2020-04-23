N=10
M=1000
ee=25438
T=300

for i in {1..300}
do
	qsub -F "${i}_${M}_${T}K ${ee} ${M}" run.sh
done
