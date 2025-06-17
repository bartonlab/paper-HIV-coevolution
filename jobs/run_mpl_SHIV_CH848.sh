#!/bin/bash -l
#SBATCH --job-name="MPL_CH848"
hostname date
mkdir -p /HIV/CH848

date
hostname
# ----------- Start Main Computation ---------------#
./bin/mpl -d HIV/CH848 -i RM6163-poly.num -o RM6163-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6163-poly-AA.dat -sn numerator-RM6163-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6163-poly.num -o RM6163-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i RM6167-poly.num -o RM6167-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6167-poly-AA.dat -sn numerator-RM6167-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6167-poly.num -o RM6167-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i RM6700-poly.num -o RM6700-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6700-poly-AA.dat -sn numerator-RM6700-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6700-poly.num -o RM6700-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i RM6713-poly.num -o RM6713-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6713-poly-AA.dat -sn numerator-RM6713-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6713-poly.num -o RM6713-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i RM6714-poly.num -o RM6714-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6714-poly-AA.dat -sn numerator-RM6714-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6714-poly.num -o RM6714-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i RM6720-poly.num -o RM6720-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6720-poly-AA.dat -sn numerator-RM6720-poly-AA.dat
./bin/mpl -d HIV/CH848 -i RM6720-poly.num -o RM6720-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH848 -i 703010848-poly.num -o 703010848-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-703010848-poly-AA.dat -sn numerator-703010848-poly-AA.dat
./bin/mpl -d HIV/CH848 -i 703010848-poly.num -o 703010848-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
date
hostname
exit 0
