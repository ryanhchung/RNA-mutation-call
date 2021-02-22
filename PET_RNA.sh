#!/bin/bash -e

# authors : Heechul Chung, M.S., Eun Jeong Cho, M.S., and Chang Ohk Sung, M.D., Ph.D., Asan Medical Center
# Mutation calling from RNA fastq (somatic)
# Mapping : STAR 2pass
# Calling : Mutect2
# Filtering : SNPiR


#Input : num or threads, output directory, and fastq name prefix.
if [ $# -lt 3 ]
then
        echo usage: $0 [prefix] [thread] [output_home]
        exit 1
fi


prefix=$1
thread=$2
output_home=$3


# reference
ref=/home/super3/Reference/GATK/Bundle/Hg19/ori_ucsc.hg19.fasta
dbsnp=/home/super3/Reference/GATK/Bundle/Hg19/dbsnp_138.hg19.vcf
mills=/home/super3/Reference/GATK/Bundle/Hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
phase1indel=/home/super3/Reference/GATK/Bundle/Hg19/1000G_phase1.indels.hg19.sites.vcf
gnomad=/home/super3/Reference/GATK/Bundle/Hg19/af-only-gnomad.raw.sites.hg19.vcf

# mapping & preprocessing
export PATH=$PATH:/data/ryanchung/programs/STAR-2.7.3a/bin/Linux_x86_64
java=/home/super3/Programs/Java/jdk1.8.0_151/bin/java
picard_path=/data/ryanchung/programs/picard-2.22.4
gatk_path=/data/ryanchung/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
GATK4=/data/ryanchung/programs/gatk-4.1.7.0/gatk

source /data/ryanchung/anaconda3/bin/activate GATK4

#SNPiR & maf path
snpir_path=/data/ryanchung/programs/SNPiR
RepeatMasker=/data/ryanchung/programs/SNPiR/snp_ir/h19.repeatmasker
gene_annotation=/data/ryanchung/programs/SNPiR/snp_ir/gencode.v19.annotation.gtf
#junction=/home/super2/reference/snp_ir/hg19_junctions_95bp.fa
rnaedit=/data/ryanchung/programs/SNPiR/snp_ir/Human_AG_all_hg19_v2.bed
extractvcf=/data/ryanchung/programs/SNPiR/snp_ir/extractvcf.edit2.R
vcf2maf=/home/super3/Programs/vcf2maf/vcf2maf.pl
vep=/home/super3/Programs/VEP/ensembl-vep


export PATH=$PATH:/data/ryanchung/programs/SNPiR:/data/ryanchung/programs/bedtools2/bin

#Use below code only if you have more than 2 fastq directories.
'''
if [ ! -f $fastq_dir/${prefix}_1.fastq ]
then
    cd $fastq_dir
    zcat ${prefix}_1.fastq.gz > ${prefix}_1.fastq & zcat ${prefix}_2.fastq.gz > ${prefix}_2.fastq
fi
'''
#RG='@RG\tID:'$6'\tPL:illumina\tSM:'$6'\tLB:'$6''


#echo "Indexing steps -> to change FASTA reference format suitable for STAR"
mkdir -p $output_home/$prefix
#cd $output_home
#mkdir $prefix
output_dir=$output_home/$prefix
cd $output_dir

mafdir=$output_home/mafdir
mkdir $mafdir
bamdir=$output_home/bamdir
mkdir $bamdir

javatmp=${output_dir}/javatmp
mkdir $javatmp

echo "Indexing steps -> to change FASTA reference format suitable for STAR"

#mkdir -p /data/ryanchung/HCC/PET/STAR_genomeDir 
#cd ${output_dir} 


#STAR --runMode genomeGenerate --genomeDir /data/ryanchung/HCC/PET/STAR_genomeDir --genomeFastaFiles $ref --runThreadN $thread
cd ${output_dir}
STAR --genomeDir /data/ryanchung/HCC/PET/STAR_genomeDir  --readFilesIn /data/HCC_PET_RNAseq/${prefix}_1.fastq /data/HCC_PET_RNAseq/${prefix}_2.fastq --runThreadN $thread

genomedir=${output_dir}/hg19_2pass
mkdir $genomedir


STAR --runMode genomeGenerate --genomeDir $genomedir --genomeFastaFiles $ref --sjdbFileChrStartEnd ${output_dir}/SJ.out.tab --sjdbOverhang 75 --runThreadN $thread


runDir=${output_dir}/STAR_output
mkdir $runDir
cd $runDir

STAR --genomeDir $genomedir --readFilesIn /data/HCC_PET_RNAseq/${prefix}_1.fastq /data/HCC_PET_RNAseq/${prefix}_2.fastq --runThreadN $thread

echo "STAR end, start generateing bam file from sam"

$java -jar ${picard_path}/picard.jar \
    AddOrReplaceReadGroups \
    I=${runDir}/Aligned.out.sam \
    O=${runDir}/rg_added_sorted.bam \
    SO=coordinate \
    RGID=${prefix} RGLB=${prefix} RGPL=illumina RGPU=${prefix} RGSM=${prefix}

$java -jar ${picard_path}/picard.jar \
    MarkDuplicates \
    I=${runDir}/rg_added_sorted.bam \
    O=${runDir}/dedupped.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=output.metrics 

$java -Djava.io.tmpdir=${javatmp} -jar ${gatk_path}/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R $ref \
    -I ${runDir}/dedupped.bam \
    -o ${runDir}/split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


#2.gatk
echo "Realignment and Recalibration steps using picard and GATK3.8. Ref : GATK3.8 best practice"
echo "realignment step 1"
$java -Djava.io.tmpdir=${javatmp} -jar ${gatk_path}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $ref \
    -known $mills \
    -known $phase1indel \
    -I ${runDir}/split.bam \
    -o ${runDir}/split.realigned.intervals

echo "realignment step 2"
$java -Djava.io.tmpdir=${javatmp} -jar ${gatk_path}/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $ref \
    -known $mills \
    -known $phase1indel \
    -I ${runDir}/split.bam \
    -targetIntervals ${runDir}/split.realigned.intervals \
    -o ${runDir}/sorted.dedup.realigned.bam



echo "recalibration step 1"
$java -Djava.io.tmpdir=${javatmp} -jar $gatk_path/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R $ref \
        -I ${runDir}/sorted.dedup.realigned.bam \
        -knownSites $dbsnp \
        -knownSites $mills \
        -knownSites $phase1indel \
        -o ${runDir}/recal.table \



echo "recalibration step 2"
$java -Djava.io.tmpdir=${javatmp} -jar $gatk_path/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R $ref \
        -I ${runDir}/sorted.dedup.realigned.bam \
        -BQSR ${runDir}/recal.table \
        -o ${runDir}/${prefix}.sorted.dedup.realigned.recal.bam \

#Haplotypecaller typically not used in somatic variations

#Mutect2
cd $output_dir
mkdir tmp_mutect2
tmpDir2=$output_dir/tmp_mutect2

echo "Mutect2 -> Call somatic SNVs, Indels via local assembly of haplotypes"

$GATK4 Mutect2 \
    --reference ${ref} \
    --input ${runDir}/${prefix}.sorted.dedup.realigned.recal.bam \
    --output ${tmpDir2}/${prefix}.Mutect2.RNA.vcf \
    --tumor-sample ${prefix} \
    --germline-resource $gnomad \
    --af-of-alleles-not-in-resource 0 \

$GATK4 FilterMutectCalls \
    -V ${tmpDir2}/${prefix}.Mutect2.RNA.vcf \
    -O ${tmpDir2}/${prefix}.Mutect2.RNA.filtered.vcf \
    -R $ref
  
egrep "^#|PASS" ${tmpDir2}/${prefix}.Mutect2.RNA.filtered.vcf > ${tmpDir2}/${prefix}.Mutect2.RNA.PASS.vcf

################################################################################
# Filtering after calling -> SNPiR
################################################################################

filter_dir=${output_dir}/SNPiR

mkdir $filter_dir

cp $tmpDir2/${prefix}.Mutect2.RNA.PASS.vcf $filter_dir/${prefix}.Mutect2.RNA.PASS.vcf

export PERL5LIB=$PERL5LIB:/data/ryanchung/programs/SNPiR

#remove unknown sequence and MT -> extract only Chr1 ~ ChrY using one of codes below. I used 2nd one because 2nd one was faster than 1st one in my case
#cat $filter_dir/${prefix}.Mutect2.RNA.PASS.vcf | egrep "#|chr[0-9X-Y]*[[:space:]]" | sed "s/^chr//g" > $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.vcf 
cat $filter_dir/${prefix}.Mutect2.RNA.PASS.vcf | egrep -v $1"(gl|chrM|^chr[0-9]*_)" > $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.vcf 

echo "### SNPiR 1/7 Converting VCF to SNPiR format ###"
# convert vcf format into custom SNPiR format and filter variants with quality <20
${snpir_path}/convertVCF.sh $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.vcf $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.txt


echo "### SNPiR 2/7 Filtering out mismatches in first 6 bp of reads ###"
# filter mismatches at read ends
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format
${snpir_path}/filter_mismatch_first6bp.pl \
        -infile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.txt \
        -outfile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.txt \
        -bamfile $runDir/${prefix}.sorted.dedup.realigned.recal.bam

echo "### SNPiR 3/7 Remove RepeatMasker sites ###"""
# filter variants in repetitive regions
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.txt | \
        intersectBed -a stdin -b $RepeatMasker -v | \
        cut -f1,3-7 > $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.txt

echo "### SNPiR 4/7 Filtering intronic candidates within 4 bp of splicing junctions ###"
# make sure your gene annotation file is in UCSC text format and sorted by chromosome and 
# transcript start position
${snpir_path}/filter_intron_near_splicejuncts.pl \
        -infile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.txt \
        -outfile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.txt \
        -genefile $gene_annotation


echo "### SNPiR 5/7 Filtering candidates in homopolymers ###"
# filter variants in homopolymers
${snpir_path}/filter_homopolymer_nucleotides.pl \
        -infile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.txt \
        -outfile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.rmhom.txt \
        -refgenome $ref


# filter variants that were caused by mismapped reads
# this may take a while depending on the number of variants to screen and the size of the reference genome
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format


echo "### SNPiR 6/7 run blat for all aligned reads for each editing sites ###"
#Check the format of bamfile, infile, and ref.(Whether the chromosome format is Chr1 or 1). IF you want to add 'chr'to your infile, use below code 
${snpir_path}/BLAT_candidates.pl \
    -infile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.rmhom.txt \
    -outfile $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.rmhom.rmblat.txt \
    -bamfile ${runDir}/${prefix}.sorted.dedup.realigned.recal.bam \
    -refgenome $ref \


echo "### SNPiR 7/7 remove known RNA editing sites ###"
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.rmhom.rmblat.txt | \
    intersectBed -a stdin -b $rnaedit -v  > \
    $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed

echo "extract variants from the raw vcf"
#Make sure to check separator of VCF and BED file.
skipn=$(cat $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.vcf | grep '#' | wc -l)
cat $filter_dir/${prefix}.Mutect2.RNA.GATK3.8.vcf | grep '#' > $filter_dir/${prefix}.Mutect2.RNA.GATK3.8_header.vcf
Rscript $extractvcf $filter_dir $skipn ${prefix}.Mutect2.RNA.GATK3.8
cat $filter_dir/${prefix}.Mutect2.RNA.GATK3.8_header.vcf $filter_dir/${prefix}.Mutect2.RNA.GATK3.8_red_variants.vcf > $filter_dir/${prefix}_final_variants.vcf

#Deactivate GATK4 python3 environment and activate perl environment
source /data/ryanchung/anaconda3/bin/deactivate
source /data/ryanchung/anaconda3/bin/activate VC.Mu1.SI

echo "make annotation & convert to maf"
perl $vcf2maf --input-vcf $filter_dir/${prefix}_final_variants.vcf --output-maf $filter_dir/${prefix}.maf --vep-data /home/super3/.vep --tumor-id ${prefix} --normal-id none --vep-path ${vep} --ref-fasta $ref --vep-forks 4 --filter-vcf /home/super3/.vep/ExAC/ExAC_nonTCGA.r1.sites.vep.vcf.gz --ncbi-build GRCh37 --cache-version 94

mv $filter_dir/${prefix}.maf $output_home/mafdir
mv $runDir/${prefix}.sorted.dedup.realigned.recal.bam $output_home/bamdir

rm -rf ${output_dir}

