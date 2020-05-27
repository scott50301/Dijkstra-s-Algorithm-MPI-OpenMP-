#!/bin/bash
##SBATCH --job-name=pl2ap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=3
#SBATCH --output=pl2ap-%j.out
#SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sjhou@scu.edu#
./dijkstra_omp ./data/100.graph 1 test.out 1
./dijkstra_omp ./data/100.graph 1 test.out 2
./dijkstra_omp ./data/100.graph 1 test.out 4
./dijkstra_omp ./data/100.graph 1 test.out 8
./dijkstra_omp ./data/100.graph 1 test.out 16
./dijkstra_omp ./data/100.graph 1 test.out 28
./dijkstra_omp ./data/500.graph 1 test.out 1
./dijkstra_omp ./data/500.graph 1 test.out 2
./dijkstra_omp ./data/500.graph 1 test.out 4
./dijkstra_omp ./data/500.graph 1 test.out 8
./dijkstra_omp ./data/500.graph 1 test.out 16
./dijkstra_omp ./data/500.graph 1 test.out 28
./dijkstra_omp ./data/1000.graph 1 test.out 1
./dijkstra_omp ./data/1000.graph 1 test.out 2
./dijkstra_omp ./data/1000.graph 1 test.out 4
./dijkstra_omp ./data/1000.graph 1 test.out 8
./dijkstra_omp ./data/1000.graph 1 test.out 16
./dijkstra_omp ./data/1000.graph 1 test.out 28
./dijkstra_omp ./data/5000.graph 1 test.out 1
./dijkstra_omp ./data/5000.graph 1 test.out 2
./dijkstra_omp ./data/5000.graph 1 test.out 4
./dijkstra_omp ./data/5000.graph 1 test.out 8
./dijkstra_omp ./data/5000.graph 1 test.out 16
./dijkstra_omp ./data/5000.graph 1 test.out 28
./dijkstra_omp ./data/10000.graph 1 test.out 1
./dijkstra_omp ./data/10000.graph 1 test.out 2
./dijkstra_omp ./data/10000.graph 1 test.out 4
./dijkstra_omp ./data/10000.graph 1 test.out 8
./dijkstra_omp ./data/10000.graph 1 test.out 16
./dijkstra_omp ./data/10000.graph 1 test.out 28
./dijkstra_omp ./data/100.graph 2 test.out 1
./dijkstra_omp ./data/100.graph 2 test.out 2
./dijkstra_omp ./data/100.graph 2 test.out 4
./dijkstra_omp ./data/100.graph 2 test.out 8
./dijkstra_omp ./data/100.graph 2 test.out 16
./dijkstra_omp ./data/100.graph 2 test.out 28
./dijkstra_omp ./data/500.graph 2 test.out 1
./dijkstra_omp ./data/500.graph 2 test.out 2
./dijkstra_omp ./data/500.graph 2 test.out 4
./dijkstra_omp ./data/500.graph 2 test.out 8
./dijkstra_omp ./data/500.graph 2 test.out 16
./dijkstra_omp ./data/500.graph 2 test.out 28
./dijkstra_omp ./data/1000.graph 2 test.out 1
./dijkstra_omp ./data/1000.graph 2 test.out 2
./dijkstra_omp ./data/1000.graph 2 test.out 4
./dijkstra_omp ./data/1000.graph 2 test.out 8
./dijkstra_omp ./data/1000.graph 2 test.out 16
./dijkstra_omp ./data/1000.graph 2 test.out 28
./dijkstra_omp ./data/5000.graph 2 test.out 1
./dijkstra_omp ./data/5000.graph 2 test.out 2
./dijkstra_omp ./data/5000.graph 2 test.out 4
./dijkstra_omp ./data/5000.graph 2 test.out 8
./dijkstra_omp ./data/5000.graph 2 test.out 16
./dijkstra_omp ./data/5000.graph 2 test.out 28
./dijkstra_omp ./data/10000.graph 2 test.out 1
./dijkstra_omp ./data/10000.graph 2 test.out 2
./dijkstra_omp ./data/10000.graph 2 test.out 4
./dijkstra_omp ./data/10000.graph 2 test.out 8
./dijkstra_omp ./data/10000.graph 2 test.out 16
./dijkstra_omp ./data/10000.graph 2 test.out 28
./dijkstra_omp ./data/100.graph 3 test.out 1
./dijkstra_omp ./data/100.graph 3 test.out 2
./dijkstra_omp ./data/100.graph 3 test.out 4
./dijkstra_omp ./data/100.graph 3 test.out 8
./dijkstra_omp ./data/100.graph 3 test.out 16
./dijkstra_omp ./data/100.graph 3 test.out 28
./dijkstra_omp ./data/500.graph 3 test.out 1
./dijkstra_omp ./data/500.graph 3 test.out 2
./dijkstra_omp ./data/500.graph 3 test.out 4
./dijkstra_omp ./data/500.graph 3 test.out 8
./dijkstra_omp ./data/500.graph 3 test.out 16
./dijkstra_omp ./data/500.graph 3 test.out 28
./dijkstra_omp ./data/1000.graph 3 test.out 1
./dijkstra_omp ./data/1000.graph 3 test.out 2
./dijkstra_omp ./data/1000.graph 3 test.out 4
./dijkstra_omp ./data/1000.graph 3 test.out 8
./dijkstra_omp ./data/1000.graph 3 test.out 16
./dijkstra_omp ./data/1000.graph 3 test.out 28
./dijkstra_omp ./data/5000.graph 3 test.out 1
./dijkstra_omp ./data/5000.graph 3 test.out 2
./dijkstra_omp ./data/5000.graph 3 test.out 4
./dijkstra_omp ./data/5000.graph 3 test.out 8
./dijkstra_omp ./data/5000.graph 3 test.out 16
./dijkstra_omp ./data/5000.graph 3 test.out 28
./dijkstra_omp ./data/10000.graph 3 test.out 1
./dijkstra_omp ./data/10000.graph 3 test.out 2
./dijkstra_omp ./data/10000.graph 3 test.out 4
./dijkstra_omp ./data/10000.graph 3 test.out 8
./dijkstra_omp ./data/10000.graph 3 test.out 16
./dijkstra_omp ./data/10000.graph 3 test.out 28