#!/bin/bash -l
#SBATCH --job-name="MPL"
hostname date
mkdir /HIV

date
hostname
# ----------- Start Main Computation ---------------#
g++ src/main.cpp src/inf.cpp src/io.cpp -O3 -march=native -lgslcblas -lgsl -o bin/mpl
./bin/mpl -d HIV -i 703010505-3-poly-seq2state.dat -o 703010505-3-poly-seq2state-MPL.dat -g 1e5 -N 1e4 -m Zanini-extended.dat -sc covariance-703010505-3-poly-seq2state.dat -sn numerator-703010505-3-poly-seq2state.dat
./bin/mpl -d HIV -i 703010505-3-poly-seq2state.dat -o 703010505-3-poly-seq2state-SL.dat  -g 1e5 -N 1e4 -m Zanini-extended.dat -nc
./bin/mpl -d HIV -i 703010848-3-poly-seq2state.dat -o 703010848-3-poly-seq2state-MPL.dat -g 1e5 -N 1e4 -m Zanini-extended.dat -sc covariance-703010848-3-poly-seq2state.dat -sn numerator-703010848-3-poly-seq2state.dat
./bin/mpl -d HIV -i 703010848-3-poly-seq2state.dat -o 703010848-3-poly-seq2state-SL.dat  -g 1e5 -N 1e4 -m Zanini-extended.dat -nc
date
hostname
exit 0
