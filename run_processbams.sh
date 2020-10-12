#!/usr/bin/bash

export PATH=/bioinfo/guests/nriddifo/miniconda2/bin:$PATH
source activate picard

bam_files=/data/kdi_prod/project_result/948/01.00/Analysis/Bwa/visor

for b in $( ls -1 "$bam_files"/*.bam )
do
    file_name=$( basename "${b}" )
    output_base=$(echo $file_name | cut -d '.' -f 1)


    picard AddOrReplaceReadGroups \
        INPUT=${b} \
        OUTPUT=$bam_files/${output_base}.tagged.filt.SC.RG.bam \
        VALIDATION_STRINGENCY=LENIENT \
        RGID=$output_base \
        RGLB=visor \
        RGPL=illumina \
        RGPU=1 \
        RGSM=$output_base
    samtools index $bam_files/${output_base}.tagged.filt.SC.RG.bam

done
