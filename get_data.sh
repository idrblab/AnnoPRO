#!/bin/bash
# Get the data from the web
wget https://data -o data.zip
mkdir PFmap/data
unzip data.zip -d PFmap/data

wget https://model_param -o model_param.zip
mkdir PFmap/model_param
unzip model_param.zip -d PFmap/model_param