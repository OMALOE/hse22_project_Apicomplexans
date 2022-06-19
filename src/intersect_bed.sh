#!/bin/bash
cd /home/leonid/PycharmProjects/python/minor_2022/project/supergroup
while IFS= read -r line1 && IFS= read -r line2 <&2 && IFS= read -r line3 <&3; do
  intersectBed -a $line1 -b $line2 -wb > $line3
done < inter_tss.txt 2< inter_z.txt 3< intersected.txt
