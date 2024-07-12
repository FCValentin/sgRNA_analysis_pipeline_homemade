#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e  ./log/
#$ -o ./log/

mkdir trimgalore
mkdir bwa
mkdir log
mkdir input
conda activate VanessaPipeline
bwa index -a is S014969_sgRNA_rabbit.fasta -p Rabbit_sgRNA #testé aussi avec bwtsw, aucun impact

#############################
#### Alignment pipeline #####
#############################

conda activate CutRun_ChenPipeline
#Remove adapter/low quality reads and export fastqc report
for f in `ls -1 fastq/*_1.fastq.gz | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do trim_galore --fastqc --quality 20 --length 20 -o trimgalore fastq/${f}_1.fastq.gz
done

conda activate ChIP_SeqPipeline
# bwa alignment on sgRNA 
for f in `ls -1 trimgalore/*_1_trimmed.fq.gz  | cut -d "/" -f 2 | sed 's/_1_trimmed.fq.gz//'`
	do bwa mem -k 2 -O 0 -t 8 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/Rabbit_sgRNA trimgalore/${f}_1_trimmed.fq.gz | samtools view -hbS  | samtools sort  > bwa/${f}.sort.bam #testé aussi avec des options de missmatch, aucun impact
done

# Check fragmentsize, and compute alignement data
for f in `ls -1 bwa/*.sort.bam  | cut -d "/" -f 2 | sed 's/.sort.bam//'`
	do 
	samtools flagstat bwa/${f}.sort.bam > log/${f}.txt
done


# Check fragmentsize, and compute alignement data
for f in `ls -1 bwa/*.bam  | cut -d "/" -f 2 | sed 's/.bam//'`
	do 
	samtools index bwa/${f}.bam
	samtools idxstats bwa/${f}.bam > log/${f}.tsv
done

#Merged all tsv samples into one sgRNAmatrix.tsv file => Excel docs "StatsAlignment.xlsx" to save and do stats on sgRNA without reads or distribution
mageck count -l sgRNA_Table.tsv -n Test --sgrna-len 20 --trim-5 0 --sample-label GeCKO_1,GeCKO_2,GeCKO_3,GeCKO_4,GeCKO_5,GeCKO_6,GeCKO_7 --fastq trimgalore/cutadapt/cutadapt_GeCKO-31_S1_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-32_S2_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-33_S3_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-34_S4_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-35_S5_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-36_S6_1.fq.gz trimgalore/cutadapt/cutadapt_GeCKO-37_S7_1.fq.gz

# -t and -c are columns names, -n output prefix
mageck test -k input/sgRNAmatrix.tsv -t D0_None_None -c D8_Algox_EGFP -n input/D0_vs_Day_8_Algox_EGFP --remove-zero both --remove-zero-threshold 0
mageck test -k input/sgRNAmatrix.tsv -t D0_None_None -c D8_Algox_FOLR1 -n input/D0_vs_Day_8_Algox_FOLR1 --remove-zero both --remove-zero-threshold 0
mageck test -k input/sgRNAmatrix.tsv -t D8_Algox_EGFP -c D8_Algox_FOLR1 -n input/Day_8_Algox_EGFP_vs_Day_8_Algox_FOLR1 --remove-zero both --remove-zero-threshold 0

mageck test -k input/sgRNAmatrix.tsv -t D8_Algox_EGFP -c D15_Algox_EGFP -n input/Day_8_Algox_EGFP_vs_Day_15_Algox_EGFP --remove-zero both --remove-zero-threshold 0
mageck test -k input/sgRNAmatrix.tsv -t D8_Algox_FOLR1 -c D15_Algox_FOLR1 -n input/Day_8_Algox_FOLR1_vs_Day_15_Algox_FOLR1 --remove-zero both --remove-zero-threshold 0

mageck test -k input/sgRNAmatrix.tsv -t D8_Algox_EGFP -c D8_UFO_EGFP -n input/Day_8_Algox_EGFP_vs_Day_8_UFO_EGFP --remove-zero both --remove-zero-threshold 0
mageck test -k input/sgRNAmatrix.tsv -t D8_Algox_FOLR1 -c D8_UFO_FOLR1 -n input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1 --remove-zero both --remove-zero-threshold 0

mageck mle --count-table input/sgRNAmatrix.tsv --design-matrix input/SampleAnnot_binary.tsv --norm-method total --output-prefix input/test.mle