#!/bin/bash

while IFS= read -r line1 && IFS= read -r line2 <&3; do
  wget -c $line1 -O $line2
  gunzip $line2
done < GCA_full_files_for_group_project.txt 3< GCA_full_names_for_files_for_group_project.txt
