#!/bin/bash

mkdir -p collapse
mkdir -p hg19_sam
mkdir -p logs/bowtie1
mkdir -p logs/pirna_map/
mkdir -p logs/cutadapt/
mkdir -p pirna_sam
mkdir -p trimmed
mkdir -p pirna_sam/filter
mkdir -p bed/putative
mkdir -p bed/ncrna_subtract_putative
mkdir -p bed/cluster
mkdir -p bed/TE_origin
mkdir -p TE_profile

#enter the sample prefix in the list. For example if your samples are named SRR1950612.fastq.gz then enter SRR1950612 in list. separate multiple entries with spaces
list=( )
sample_dir_prefix="./trimmed/"
out_dir_prefix="./"
tools_dir="/tools/"

#please enter the number of samples in your study
num_sample=0

#please modify the cutadapt and the corresponding adapters, bowtie index path (assuming that you are using human genome). You can skip this part and do your alignments separately.
#We have created the bowtie indexes of the pirnas and provided them in the tools directory. Change that in the 4th step where you have extracted them. It also assumes that you have samtools in the PATH
for prefix in "${list[@]}"
do
	echo $prefix >> prefix_list
	echo "Adapter trimming"
	cutadapt-1.8.1/bin/cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -m 15 -M 76 -o $out_dir_prefix"trimmed/"$prefix.fastq.gz $out_dir_prefix"fastq/"$prefix.fastq.gz > $out_dir_prefix"logs/cutadapt/"$prefix.txt
	echo "Running Collapse"
	python $tools_dir"collapse.py" -i $sample_dir_prefix$prefix".fastq.gz" -o $out_dir_prefix"collapse/"$prefix".fa" --format fastq
	echo "Mapping to hg19"
	gzip -d -c $out_dir_prefix"collapse/"$prefix".fa.gz" | bowtie -p 8 <hg19 index> -f - -k 300 -l 22 -v 1 -n 1 -S ./hg19_sam/$prefix.sam 2> ./logs/bowtie1/"$prefix"_bowtie1.log
	echo "piRNA alignment"
	gzip -d -c $out_dir_prefix"collapse/"$prefix".fa.gz" | bowtie -p 8 tools/pirna/pirna_uniq --norc -f - -a -v 1 -n 0 -l 26 -S ./pirna_sam/$prefix.sam 2> ./logs/pirna_map/$prefix.log
	echo "Filtering SAM"
	samtools view -F 4 -S $out_dir_prefix"pirna_sam/"$prefix.sam | python $tools_dir"sam_filter.py" > $out_dir_prefix"pirna_sam/filter/"$prefix"_filtered.sam"

done


cut -f10 ./pirna_sam/filter/* | sort | uniq > $out_dir_prefix"known_pirna"

for prefix in "${list[@]}"
do
	echo $prefix
	samtools view -h -F 4 $out_dir_prefix"hg19_sam/"$prefix.sam | python $tools_dir"annotate_sam.py" $out_dir_prefix"known_pirna" | samtools view -bh -S - > $out_dir_prefix"hg19_sam/"$prefix.bam
	rm $out_dir_prefix"hg19_sam/"$prefix.sam
	echo "Making BED files"
	samtools view ./hg19_sam/$prefix.bam | python $tools_dir"extract_read.py" | sort -V -k1,1 -k2,2 > $out_dir_prefix"bed/putative/"$prefix.bed
	#Here we have created a bed file of non-coding RNAs and we are subtracting them out. My apologies that we did not provide them with the supplementary files. You can obtain them from Ensemble website and convert it into a BED format
	bedtools subtract -b <ncrna.bed> -a $out_dir_prefix"bed/putative/"$prefix.bed -s -f 0.5 -A > $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed
	python $tools_dir"pilfer.py" -i $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed > $out_dir_prefix"bed/cluster/"$prefix
done

python union.py prefix_list /bed/cluster/ > cluster_union.tsv
awk 'NR>1{print $1}' cluster_union.tsv | tr ":" "\t" | tr "-" "\t" | sort -V -k1,1 -k2,2 > $out_dir_prefix"bed/clusters.bed"
bedtools merge -i  $out_dir_prefix"bed/clusters.bed" > $out_dir_prefix"bed/merge_clusters.bed"

echo "Merging clusters"
python $tools_dir"merge_cluster.py" $out_dir_prefix"bed/merge_clusters.bed" $out_dir_prefix"cluster_union.tsv" $num_sample > $out_dir_prefix"cluster_union_merged.tsv"

echo "TE profiling"
#We are doing an optional TE profiling here and we have provided the corresponding transposon bed files in the tools directory. Please extract it and change the location to the extracted location
for prefix in "${list[@]}"
do
	echo $prefix
	bedtools intersect -a <path/>retro-transposons.bed -b $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed -s -F 1 -wa -wb | awk 'BEGIN{OFS="\t"}{print $7,$8,$9,$10,$11,$12,$4}' > $out_dir_prefix"bed/TE_origin/"$prefix
done

cut -f4 "$out_dir_prefix"bed/ncrna_subtract_putative/* | sort | uniq | awk 'BEGIN{FS="::"}$2=="PU"{print $1}' | python  $tools_dir"merge.py" | rev | sort | rev | python  $tools_dir"merge.py" | awk '{print ">seq"NR"|PU\n"$1}' > $out_dir_prefix"putative.fa"
cut -f4 "$out_dir_prefix"bed/ncrna_subtract_putative/* | sort | uniq | awk 'BEGIN{FS="::"}$2=="PI"{print $1}' | python  $tools_dir"merge.py" | rev | sort | rev | python  $tools_dir"merge.py" | awk '{print ">seq"NR"|PI\n"$1}' > $out_dir_prefix"known.fa"

cat $out_dir_prefix"putative.fa" $out_dir_prefix"known.fa" > all_pirna.fa

#Similarly we created a bowtie 2 index of the transposons we used in the previous step. You can obtain the fasta using the getfasta module of bedtools and then create an index for it to map
bowtie2 -p 8 -k 100 --local -x <transposon index> -f -U $out_dir_prefix"all_pirna.fa" --nofw --no-unal --no-hd -S $out_dir_prefix"known_putative.sam" 2> $out_dir_prefix"logs/TE_target.txt"
