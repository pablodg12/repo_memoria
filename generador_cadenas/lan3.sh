M=10000
ee=254963
T=300
for N in {100..100..10}
do
        for i in {1..10}
        do
                sbatch run2.sh "${i}_${N}_${M}_${T}K_25 ${ee} ${M} ${N} /home/pibarra/respaldo2/restart_${i}_${N}_${M}_${T}K_25.npy"
                echo "${i}_${N}_${M}_${T}K_25_e"
        done
done
