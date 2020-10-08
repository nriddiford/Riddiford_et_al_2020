sim_svs=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/data/sv_sim1.bed
genome=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa
out_dir=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/SV_paper_20/visor

mkdir -p $out_dir

rm -r $out_dir/*

VISOR HACk -bed $sim_svs -g $genome -o $out_dir 
