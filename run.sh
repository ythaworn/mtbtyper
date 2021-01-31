#!/bin/bash

# timer=$SECONDS

# python3 mtbtyper.py -h
python3 mtbtyper.py vcf -o lineage #--final_prediction

# timer=$(($SECONDS-timer))
# printf "Time used: %02d:%02d:%02d\n" "$((timer/3600))" "$((timer/60%60))" "$((timer%60))"

