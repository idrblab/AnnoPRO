#!/bin/bash
# Get the data from the web
wget https://promap.oss-cn-hangzhou.aliyuncs.com/data.tar.gz
tar -zxvf data.tar.gz -C PFmap

wget https://promap.oss-cn-hangzhou.aliyuncs.com/model_param.tar.gz
tar -zxvf model_param.tar.gz -C PFmap
