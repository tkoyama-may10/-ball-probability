#!/bin/bash

echo "" > tab3.txt
for i in $(seq 5 5 100)
do 
    echo $i
    /usr/bin/time -o tmp -v ./a.out 9 $i 20.0 > tmp2 
    cat tmp  >> tab3.txt
    cat tmp2 >> tab3.txt
done


cat tab3.txt | grep Probability > tmp
cat tab3.txt | grep User > tmp2
R -f tab3.R