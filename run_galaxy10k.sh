#!/bin/bash
ncs=(2 4 8)
steps=(50 75 100)
file="galaxy10k.txt"
input="inputs/${file}"

tuantu_file="${file}_tuantu"
songsong1may_file="${file}_songsong1may"
songsong2may_file="${file}_songsong2may"


echo "" > ${tuantu_file}

for step in ${steps[*]}
do
	echo "-------------step ${step} - ----------------"
	echo "-------------step ${step} - ----------------" >> ${tuantu_file}
	mpirun bh -i ${input} -m s -s ${step} >> ${tuantu_file}
done

echo "" > ${songsong1may_file}

for nc in ${ncs[*]}
do
	for step in ${steps[*]}
	do
		echo "--------------nc ${nc} step ${step} -------------------"
		echo "--------------nc ${nc} step ${step} -------------------" >> ${songsong1may_file}
		mpirun -np ${nc} bh -i ${input} -s ${step} >> ${songsong1may_file}
	done
done

echo "" >> ${songsong2may_file}
for nc in ${ncs[*]}
do
	for step in ${steps[*]}
	do
		echo "--------------nc ${nc} step ${step} -------------------"
		echo "--------------nc ${nc} step ${step} -------------------" >> ${songsong2may_file}
		mpirun -np ${nc} --hostfile hostfile --host localhost,client bh -i ${input} -s ${step} >> ${songsong2may_file}
	done
done
