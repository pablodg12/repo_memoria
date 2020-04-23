M=1000
ee=25438
T=300
for N in {20,20}
do
	for i in {1..1}
	do
        	sbatch run.sh "${i}_${N}_${M}_${T}K_25 ${ee} ${M} ${N}"
		echo "${i}_${N}_${M}_${T}K_25_e"
	done
done
