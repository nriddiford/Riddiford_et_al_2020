#!/usr/bin/bash
genome=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa
out_dir=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor

tumour_h1=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/tumour_h1.bed
tumour_h2=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/tumour_h2.bed
tumour_regions=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/tumour_regions.bed

normal_h1=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/normal_h1.bed
normal_h2=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/normal_h2.bed
normal_regions=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/normal_regions.bed


tumour_out=${out_dir}/tumour
normal_out=${out_dir}/normal

tumour_coverage=40
normal_coverage=70


setup() {
    # mkdir -p $out_dir $tumour_out $normal_out

    rm -r $tumour_out
    rm -r $normal_out

    perl -p -i -e 's/ /\t/g' $tumour_h1
    perl -p -i -e 's/ /\t/g' $tumour_h2
    perl -p -i -e 's/ /\t/g' $normal_h1
    perl -p -i -e 's/ /\t/g' $normal_h2

    cat $tumour_h1 $tumour_h2 > $out_dir/data/tumour_svs.bed
    cat $normal_h1 $normal_h2 > $out_dir/data/normal_svs.bed

}


runvisor() {

    export PATH=/bioinfo/guests/nriddifo/miniconda3/bin:$PATH
    source activate visorenv

    echo
    """
    Running for tumour sample

    VISOR HACk -bed $tumour_h1 $tumour_h2 -g $genome -o $tumour_out
    VISOR SHORtS --length 100 --error 0.00001 --coverage $tumour_coverage -s $tumour_out -bed $tumour_regions -g $genome  -o $tumour_out/bam
    """

    VISOR HACk -bed $tumour_h1 $tumour_h2 -g $genome -o $tumour_out

    VISOR SHORtS --length 100 --error 0.00001 --coverage $tumour_coverage -s $tumour_out -bed $tumour_regions -g $genome -o $tumour_out/bam

    mv $tumour_out/bam/sim.srt.bam $tumour_out/bam/visorR01.bam
    samtools index $tumour_out/bam/visorR01.bam

    echo """
    Running for normal sample

    VISOR HACk -bed $normal_h1 $normal_h2 -g $genome -o $normal_out
    VISOR SHORtS --length 100 --error 0.00001 --coverage $normal_coverage -s $normal_out -bed $normal_regions -g $genome  -o $normal_out/bam
    """

    VISOR HACk -bed $normal_h1 $normal_h2 -g $genome -o $normal_out

    VISOR SHORtS --length 100 --error 0.00001 --coverage $normal_coverage -s $normal_out -bed $normal_regions -g $genome -o $normal_out/bam

    mv $normal_out/bam/sim.srt.bam $normal_out/bam/visorR02.bam
    samtools index $normal_out/bam/visorR02.bam

}


cleanup(){
    export PATH=/bioinfo/guests/nriddifo/miniconda2/bin:$PATH
    source activate picard

    bwa_out=/data/kdi_prod/project_result/948/01.00/Analysis/Bwa/visor

    picard AddOrReplaceReadGroups \
        INPUT=$tumour_out/bam/visorR01.bam \
        OUTPUT=$bwa_out/visorR01.tagged.filt.SC.RG.bam \
        VALIDATION_STRINGENCY=LENIENT \
        RGID=visorR01 \
        RGLB=visor \
        RGPL=illumina \
        RGPU=1 \
        RGSM=visorR01

    samtools index $bwa_out/visorR01.tagged.filt.SC.RG.bam

    picard AddOrReplaceReadGroups \
        INPUT=$normal_out/bam/visorR02.bam \
        OUTPUT=$bwa_out/visorR02.tagged.filt.SC.RG.bam \
        VALIDATION_STRINGENCY=LENIENT \
        RGID=visorR02 \
        RGLB=visor \
        RGPL=illumina \
        RGPU=1 \
        RGSM=visorR02

    samtools index $bwa_out/visorR02.tagged.filt.SC.RG.bam
}

setup
runvisor
cleanup
