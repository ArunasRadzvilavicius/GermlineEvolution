#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:50:0 
#$ -l mem=500M
#$ -l tmpfs=50M
#$ -t 1-400
#$ -N GInvasion
#$ -wd /home/ucbprad/scratch/FINAL/GInvasion
fname=$(printf "1D-%.4f-%.4f-$ginv-to-$Gg-$Gt-$Gs-$M-$IM-Q$Q-L$L.evo" $mus $mug)
echo $fname
./delta $Gg $Gt $Gs $M $mus $mug $ginv $IM $ginv $L $Q >> "$fname"
rm GI*
