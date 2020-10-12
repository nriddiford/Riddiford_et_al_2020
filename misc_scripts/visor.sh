tumour_svs=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/tumour_svs.bed
normal_svs=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/normal_svs.bed

bam_regions=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor/data/tumour_regions.bed

genome=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa
out_dir=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor

mkdir -p $out_dir

#rm -r $out_dir/*

# Run for tumour
VISOR HACk -bed $tumour_svs -g $genome -o $out_dir/tumour
VISOR HACk -bed $normal_svs -g $genome -o $out_dir/normal

VISOR SHORtS --error 0.0001 --coverage 40 --length 100 -s $out_dir/tumour -bed $bam_regions -g $genome -o $out_dir/tumour/bam
VISOR SHORtS --error 0.0001 --coverage 40 --length 100 -s $out_dir/normal -bed $bam_regions -g $genome -o $out_dir/normal/bam
