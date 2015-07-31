#!/bin/bash

para=1

if [ $# -eq 1 ];then
  para=$1
fi


ls inp|parallel -j $para cpmdRun.sh
