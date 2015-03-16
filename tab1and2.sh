#!/bin/bash

echo hirotsu_1_1: > tab1and2.txt
for i in $(seq 10 20)
do 
    echo $i
    /usr/bin/time -o tmp -v ./a.out 5 $i 40.0 > tmp2 
    cat tmp  >> tab1and2.txt
    cat tmp2 >> tab1and2.txt
done

echo hirotsu_1_2: >> tab1and2.txt
for i in $(seq 10 20)
do 
    echo $i
    /usr/bin/time -o tmp -v ./a.out 6 $i 40.0 > tmp2 
    cat tmp  >> tab1and2.txt
    cat tmp2 >> tab1and2.txt
done

echo hirotsu_2_2: >> tab1and2.txt
for i in $(seq 10 20)
do 
    echo $i
    /usr/bin/time -o tmp -v ./a.out 7 $i 40.0 > tmp2 
    cat tmp  >> tab1and2.txt
    cat tmp2 >> tab1and2.txt
done

echo hirotsu_2_2: >> tab1and2.txt
for i in $(seq 10 20)
do 
    echo $i
    /usr/bin/time -o tmp -v ./a.out 8 $i 40.0 > tmp2 
    cat tmp  >> tab1and2.txt
    cat tmp2 >> tab1and2.txt
done

cat tab1and2.txt | grep Probability > tmp
cat tab1and2.txt | grep User > tmp2
R -f tab1and2.R
