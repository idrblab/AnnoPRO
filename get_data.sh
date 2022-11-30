#!/bin/bash
# Get the data from the web
wget https://promap.oss-cn-hangzhou.aliyuncs.com/data.tar.gz -o data.tar.gz
mkdir PFmap/data
tar -zxvf data.tar.gz -C PFmap/data

wget https://promap.oss-cn-hangzhou.aliyuncs.com/model_param.tar.gz -o model_param.tar.gz
mkdir PFmap/model_param
tar -zxvf model_param.tar.gz -C PFmap/model_param
