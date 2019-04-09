#!/bin/bash
# PRL Project 3
# @author Michal Ormos
# @email xormos00@stud.fit.vutbr.cz
# @date April 2019

if [ "$#" -gt 1 ]; then
    echo "Illegal number of parameters"
    exit 1
elif [ "$#" -eq 0 ]; then
    echo "Illegal number of parameters"
    exit 1
# elif [ $1 -lt 1 ]; then
#     echo "invalid amount of numbers"
#     exit 1 
elif [ "$#" -eq 1 ]; then
    numbers=$1;
else
    numbers="ABCDE"
fi

# compile
mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp

input=$1;
input_length=${#input};
let processes=2*$input_length-2;
if [ $processes == 0 ]; then # in case of a single item, run only one processor
    processes=1;
fi
#run

mpirun --prefix /usr/local/share/OpenMPI -np $processes vuv $numbers

rm -f vuv