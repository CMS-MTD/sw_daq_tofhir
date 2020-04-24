#!/bin/bash
for i in `seq 14461 15717`;
#for i in `seq 15600 15717`;
do
	echo $i
	./ConvertTOFPETSinglesToEvents $i 
	echo "\n\n"
done

