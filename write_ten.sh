#!/bin/bash
# Loop through specific files and concatenate to create other files
for i in 0 5 10 12 15 17 20 22 25 30; do
  awk -v file_num=$i '{print file_num, $2}' w${i}.dat >> benzene_10_colvar.dat
done

for i in 0 10 15 20 25 30; do
  awk -v file_num=$i '{print file_num, $2}' w${i}.dat >> benzene_6_colvar.dat
done

for i in 0 5 10 15 20 25 30; do
  awk -v file_num=$i '{print file_num, $2}' w${i}.dat >> benzene_7_colvar.dat
done

