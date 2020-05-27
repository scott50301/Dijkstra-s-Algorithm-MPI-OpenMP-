#!/bin/bash
echo 1
mpirun -n 1 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 1 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 1 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 1 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 1 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.out
echo 2
mpirun -n 2 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 2 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 2 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 2 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 2 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.out
echo 4
mpirun -n 4 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 4 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 4 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 4 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 4 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.ouecho
echo 8
mpirun -n 8 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 8 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 8 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 8 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 8 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.out
echo 16
mpirun -n 16 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 16 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 16 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 16 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 16 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.out
echo 28
mpirun -n 28 ./dijkstra_mpi ./data/100.graph 1 test_mpi.out
mpirun -n 28 ./dijkstra_mpi ./data/500.graph 1 test_mpi.out
mpirun -n 28 ./dijkstra_mpi ./data/1000.graph 1 test_mpi.out
mpirun -n 28 ./dijkstra_mpi ./data/5000.graph 1 test_mpi.out
mpirun -n 28 ./dijkstra_mpi ./data/10000.graph 1 test_mpi.out