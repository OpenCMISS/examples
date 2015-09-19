#!/bin/bash
UNAME=$1
lines=$(squeue -u $UNAME -h|grep -v suspended|sed '/^$/d'|wc -l)
while [ $lines -ge 1 ];
do
#echo "Waiting.... $lines"
sleep 1
lines=$(squeue -u $UNAME -h|grep -v suspended|sed '/^$/d'|wc -l)
done
#echo "DONE"
