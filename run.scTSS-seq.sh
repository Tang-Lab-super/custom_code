f=scTSS.fq.gz
fi=`basename $f`
file=${fi/.fq.gz/}
cutadapt -O 8 --discard-untrimmed -g TTTCTTATATGGG -o ${file}.cuted.fq.gz ${f} > log.cutadapt.${file}.txt

Fastq_R1=${file}.cuted.fq.gz
hisat2=/data/public/software/hisat2-2.1.0/hisat2
samtools="/data/public/software/samtools-1.3.1/samtools"
GENOME=/data/public/refGenome/hisat2_index/hg38/hg38
GENOME_SIZE=/data/public/refGenome/bwa_index/hg38/ChromInfo.txt
${hisat2} -x ${GENOME} -U ${Fastq_R1} -p ${cpu_num} -S ${PREFIX}_R1.cuted.sam 2>${PREFIX}_R1.cuted.bwa.log.txt
cat ${PREFIX}_R1.cuted.sam | \
    perl -lane 'if($_ =~ /^@SQ/ && $F[1] =~ /^SN\:chr[0-9XYAB]+$/){print $_} elsif($F[2] =~ /^chr[0-9XYAB]+$/){print $_}' | \
    ${samtools} view -b -o ${PREFIX}_R1.cuted.bam -q 20 -F 256 -F 2048 -F 1024 - 2>${PREFIX}_R1.cuted.sam2bam.log.txt
${samtools} sort --threads 20 -o ${PREFIX}_R1.cuted.sorted.bam ${PREFIX}_R1.cuted.bam

#########################################################################
scRNA-seq=scRNA-seq.fq.gz
zcat ${scRNA-seq} |\
 awk 'NR%4==1 || NR%4==2{print $1}' | sed -e "N;s/\n/\t/g" -e "s/@//g" |\
 awk '{bc=substr($2,0,16);umi=substr($2,17,10);pr=substr($2,27,13);print $1"\t"bc"\t"umi"\t"pr}' \
 >> ${PREFIX}.UMI.txt

bamToBed -i ${PREFIX}_R1.cuted.sorted.sorted.bam |\
    awk 'NR==FNR{a[$1]=$2"\t"$3}NR>FNR{print $0"\t"a[$4]}' ${PREFIX}.UMI.txt - |\
    awk '{a[$7"\t"$8]=$0}END{for(i in a){print a[i]}}' | sort -k1,1 -k2,2n |\
    bedToBam -i - -g hg38.fa.fai \
    > ${PREFIX}_R1.cuted.merge.sorted.bam


