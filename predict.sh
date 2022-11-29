#!/bin/bash 

file=$1

diamond blastp --db /home/zhengly/PFmap/PFmap/data/cafa4.dmnd --query $1/input-protein.dat --out $1/case.txt

python /home/zhengly/PFmap/PFmap/data_procession/main.py --file_path=$1

nohup python /home/zhengly/PFmap/predict/main.py --file_path=$1 > $1/case.log

