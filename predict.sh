#!/bin/bash 

file=$1
diamond=$2

diamond blastp --db PFmap/data/cafa4.dmnd --query $1/input-protein.dat --out $1/case.txt 

python PFmap/data_procession/main.py --file_path=$1

nohup python main.py --file_path=$1 --with_diamond $2 > $1/case.log

