#!/bin/bash
cd /home/leonid/PycharmProjects/python/minor_2022/project/supergroup
for filename in *_z.bed; do
    sort -k1,1 -k2,2n $filename > sorted_$filename
    bedtools merge -i sorted_$filename  -c 5 -o max > merged_$filename
    awk '{print $1, $2, $3, $1"_"$2 , $4}' merged_$filename   > final_$filename
    rm sorted_$filename
    rm merged_$filename
done

