#!/bin/bash -l
#SBATCH --job-name="MPL_CH505"
hostname date
mkdir -p /HIV/CH505

date
hostname
# ----------- Start Main Computation ---------------#
./bin/mpl -d HIV/CH505 -i 703010505-poly.num -o 703010505-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-703010505-poly-AA.dat -sn numerator-703010505-poly-AA.dat
./bin/mpl -d HIV/CH505 -i 703010505-poly.num -o 703010505-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM5695-poly.num -o RM5695-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM5695-poly-AA.dat -sn numerator-RM5695-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM5695-poly.num -o RM5695-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6070-poly.num -o RM6070-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6070-poly-AA.dat -sn numerator-RM6070-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6070-poly.num -o RM6070-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6072-poly.num -o RM6072-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6072-poly-AA.dat -sn numerator-RM6072-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6072-poly.num -o RM6072-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6697-poly.num -o RM6697-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6697-poly-AA.dat -sn numerator-RM6697-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6697-poly.num -o RM6697-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6699-poly.num -o RM6699-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6699-poly-AA.dat -sn numerator-RM6699-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6699-poly.num -o RM6699-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6701-poly.num -o RM6701-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6701-poly-AA.dat -sn numerator-RM6701-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6701-poly.num -o RM6701-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
./bin/mpl -d HIV/CH505 -i RM6703-poly.num -o RM6703-poly-AA-MPL.dat -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -sc covariance-RM6703-poly-AA.dat -sn numerator-RM6703-poly-AA.dat
./bin/mpl -d HIV/CH505 -i RM6703-poly.num -o RM6703-poly-AA-SL.dat  -g 1e5 -N 1e4 -q 22 -m mutation_rate.txt -nc
date
hostname
exit 0
