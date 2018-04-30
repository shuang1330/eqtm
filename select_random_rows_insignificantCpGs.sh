#!/bin/bash

export PROJECT_ROOTDIR="/home/shuang/projects/eqtm"


{ IFS= read -r head; echo "$head"; shuf | head -n 20000; } < $PROJECT_ROOTDIR/data/eqtmZscores/2017-12-09-eQTLsFDR-gt0.05-flipped.txt > $PROJECT_ROOTDIR/data/eqtmZscores/random20000-eQTLsFDR-gt0.05-flipped.txt 
