qsub -F "1_10000_500K 187221 10000 restart_1.npy" run2.sh
qsub -F "2_10000_500K 187221 10000 restart_2.npy" run2.sh
qsub -F "3_10000_500K 187221 10000 restart_4.npy" run2.sh
qsub -F "1_10000_300K 254963 10000 restart_3.npy" run2.sh
qsub -F "2_10000_300K 254963 10000 restart_5.npy" run2.sh
qsub -F "3_10000_300K 254963 10000 restart_hola.npy" run2.sh
