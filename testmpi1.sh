#!/bin/bash
#
#SBATCH --job-name=dijkstra_mpi
#SBATCH --nodes=3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=84
#SBATCH --mem-per-cpu=16G
#SBATCH --output=dijkstra_mpi-%j.out 
#SBATCH --time=23:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sjhou@scu.edu #

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra_mpi-$current_time.log
file_path=/WAVE/users/unix/sjhou/Dijkstra-s-Algorithm-MPI-OpenMP- 
module load OpenMPI/3.1.4-GCC-8.3.0

for size in 100 500 1000 5000 10000; do
    for thread_num in 1 2 4 8 16 28; do
    # for k in seq 1 20; do   
        for k in seq 1 5; do
            source=$RANDOM
            let "source %= $size"
            mpirun -n $thread_num --map-by node $file_path/dijkstra_mpi "$file_path/data/${size}.graph" $source "$file_path/data/${size}_serial.out" >> "$file_path/logs/$log_file_name"
        done 
    done
done
