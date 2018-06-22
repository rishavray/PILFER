#!/bin/bash

mkdir -p collapse
mkdir -p mm10_sam
mkdir -p logs/bowtie1
mkdir -p logs/pirna_map/
mkdir -p logs/cutadapt/
mkdir -p pirna_sam
mkdir -p pirna_sam/filter
mkdir -p bed/putative
mkdir -p bed/ncrna_subtract_putative
mkdir -p bed/cluster
mkdir -p bed/TE_origin
mkdir -p TE_profile

#Enter sample names list without extension. For example, if you have two samples named Ciccio.fastq.gz and Pasticcio.fastq.gz, enter Ciccio and Pasticcio. Separate names with space.
list=(DIV0-1)

sample_dir_prefix="./trimmed/"
out_dir_prefix="./"
tools_dir="./tools/"

#Enter number of samples analyzed.
num_sample=1

for prefix in "${list[@]}"
do
	echo $prefix >> prefix_list
	echo "Skipping Adapter trimming"
	echo "Running Collaps for $prefix"
	python $tools_dir"collapse.py" -i $sample_dir_prefix$prefix".fastq.gz" -o $out_dir_prefix"collapse/"$prefix".fa" --format fastq
	echo "Mapping collapsed $prefix reads to mm10"
	gzip -d -c $out_dir_prefix"collapse/"$prefix".fa" | bowtie -p 8 <index> -f - -k 300 -l 22 -v 1 -n 1 -S ./mm10_sam/$prefix.sam 2> ./logs/bowtie1/"$prefix"_bowtie1.log
	echo "piRNA alignment for $prefix"
	gzip -d -c $out_dir_prefix"collapse/"$prefix".fa.gz" | bowtie -p 8 <piRNA index> --norc -f - -a -v 1 -n 0 -l 26 -S ./pirna_sam/$prefix.sam 2> ./logs/pirna_map/$prefix.log
	echo "Filtering SAM for $prefix"
	samtools view -F 4 -S $out_dir_prefix"pirna_sam/"$prefix.sam | python $tools_dir"sam_filter.py" > $out_dir_prefix"pirna_sam/filter/"$prefix"_filtered.sam"

done


cut -f10 ./pirna_sam/filter/* | sort | uniq > $out_dir_prefix"known_pirna"

for prefix in "${list[@]}"
do
	samtools view -h -F 4 $out_dir_prefix"mm10_sam/"$prefix.sam | python $tools_dir"annotate_sam.py" $out_dir_prefix"known_pirna" | samtools view -bh -S - > $out_dir_prefix"mm10_sam/"$prefix.bam
	rm $out_dir_prefix"mm10_sam/"$prefix.sam
	echo "Making BED files for $prefix"
	samtools view ./mm10_sam/$prefix.bam | python $tools_dir"extract_read.py" | sort -V -k1,1 -k2,2 > $out_dir_prefix"bed/putative/"$prefix.bed
	#SUBTRACT known ncRNA (GENCODE snRNA-snoRNA-miRNA-rRNA)
	bedtools subtract -b <path_to/>gencodeVM17-ncRNA.bed -a $out_dir_prefix"bed/putative/"$prefix.bed -s -f 0.5 -A > $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed
	python $tools_dir"pilfer.py" -i $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed > $out_dir_prefix"bed/cluster/"$prefix
done

python $tools_dir"union.py" prefix_list $out_dir_prefix"/bed/cluster/" > cluster_union.tsv
awk 'NR>1{print $1}' cluster_union.tsv | tr ":" "\t" | tr "-" "\t" | sort -V -k1,1 -k2,2 > $out_dir_prefix"bed/clusters.bed"
bedtools merge -i  $out_dir_prefix"bed/clusters.bed" > $out_dir_prefix"bed/merge_clusters.bed"

echo "Merging clusters for $prefix"
python $tools_dir"merge_cluster.py" $out_dir_prefix"bed/merge_clusters.bed" $out_dir_prefix"cluster_union.tsv" $num_sample > $out_dir_prefix"cluster_union_merged.tsv"

for prefix in "${list[@]}"
do
	echo "REs profiling for $prefix"
	bedtools intersect -a <path_to/>retro-mm10.bed -b $out_dir_prefix"bed/ncrna_subtract_putative/"$prefix.bed -s -wa -wb | awk 'BEGIN{OFS="\t"}{print $7,$8,$9,$10,$11,$12,$4}' > $out_dir_prefix"bed/TE_origin/"$prefix
done

cut -f4 "$out_dir_prefix"bed/ncrna_subtract_putative/* | sort | uniq | awk 'BEGIN{FS="::"}$2=="PU"{print $1}' | python  $tools_dir"merge.py" | rev | sort | rev | python  $tools_dir"merge.py" | awk '{print ">seq"NR"|PU\n"$1}' > $out_dir_prefix"putative.fa"
cut -f4 "$out_dir_prefix"bed/ncrna_subtract_putative/* | sort | uniq | awk 'BEGIN{FS="::"}$2=="PI"{print $1}' | python  $tools_dir"merge.py" | rev | sort | rev | python  $tools_dir"merge.py" | awk '{print ">seq"NR"|PI\n"$1}' > $out_dir_prefix"known.fa"

cat $out_dir_prefix"putative.fa" $out_dir_prefix"known.fa" > all_pirna.fa

#Map piRNA to REs to identify target REs.
bowtie2 -p 8 -k 100 --local -x <bowtie2 index> -f $out_dir_prefix"all_pirna.fa" --nofw --no-unal -S $out_dir_prefix"known_putative.sam" 2> $out_dir_prefix"logs/TE_target.txt"

mkdir ./$prefix
mv all_pirna.fa bed/ cluster_union* collapse/ mm10_sam/ known.fa known_pirna known_putative.sam logs/ pirna_sam/ prefix_list putative.fa TE_profile/ ./$prefix

echo Finished running pipeline for $prefix
