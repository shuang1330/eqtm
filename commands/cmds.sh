#!/usr/bin/env bash


# interactive node
srun -N 1 --cpus-per-task=1 --qos=dev --mem=50gb --time=24:00:00 --pty bash

find . -name "*.txt"
