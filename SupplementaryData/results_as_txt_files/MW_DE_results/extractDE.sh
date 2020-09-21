#/bin/bash

#logFC cutoff
upcutoff=0.3785 #log_2(1.3) fold
dwncutoff=$(bc <<< -1*$upcutoff)
fdr=0.05

outdir="DE_genes_1.3"

mkdir -p $outdir

files=$(ls *_mw.txt)
for f in $files 
do
	echo "$f"
	name=$(basename $f .txt)
	#echo "awk '$4>=$upcutoff && $6<$fdr' $f > "${outdir}/${name}_up.txt""
	awk -v a="$upcutoff" -v b="$fdr" '$4>=a && $6<b' $f > "${outdir}/${name}_up.txt"
	awk -v c="$dwncutoff" -v b="$fdr" '$4<=c && $6<b' $f > "${outdir}/${name}_dwn.txt"
	
done
