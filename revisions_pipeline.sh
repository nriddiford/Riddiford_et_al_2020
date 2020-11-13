unmappable=~/iCloud/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed

cnv_calls=/Volumes/perso/Analysis/Analysis/CNV-Seq/visor/results/w_500/
cnv_calls=~/Desktop/visorCNV/

bash ~/iCloud/Desktop/script_test/svParser/script/freecFilt.sh -e $unmappable -d .


# svParser
bash runParser \
  -d ~/Desktop/SV_paper_20/svParser/visor/data/ \
  -o ~/Desktop/SV_paper_20/svParser/visor/ \
  -c $cnv_calls \
  -fma \
  -e $unmappable \
  -g /Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf


# svSupport
bash scripts/run_all.sh \
  -v ~/Desktop/SV_paper_20/svParser/visor/summary/merged/ \
  -p data/tumour_purity.txt \
  -o /Users/Nick/Desktop/SV_paper_20/svSupport \
  -cs


# svStitch
script_bin=/Users/Nick/iCloud/Desktop/script_test/svParser/script

for file in *_svSupport.txt
do
  python $script_bin/svStitch.py -i $file \
  -w 5500 \
  -o $(pwd)
done

# Remove false positives and reannotate breakpoints
# bash runParser -d ~/Desktop/final_analysis/data/ \
#   -o ~/Desktop/SV_paper_20/svParser/visor/ \
#   -e $unmappable \
#   -g /Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf \
#   -r

out_dir=~/Desktop/SV_paper_20/svParser/visor/summary/merged
for file in ${out_dir}/*_stitched.txt
  do
  stem=$(basename "${file}")
  output_base=$(echo $stem | cut -d '_' -f 1)
  echo "Running for sample $output_base"
  if [ -f $out_dir/${output_base}_mh.txt ]
  then
    rm $out_dir/${output_base}_mh.txt
  fi
  python $script_bin/extract_vars.py -v $file -o ${out_dir}/${output_base}_mh.txt --write_genes --write_breakpoints
done

rm all_bps_merged.txt
rm all_genes_merged.txt
cat *_mh.txt > all_bps_merged.txt
cat *_hit_genes.txt > all_genes_merged.txt

python $script_bin/merge_files.py -e _stitched.txt -d $out_dir -o all_samples_merged.txt

rm *_mh.txt
rm *_hit_genes.txt






#### Wholegut

# svSupport
bash scripts/run_all.sh \
  -v /Users/Nick/iCloud/Desktop/parserTest/wholegut_71020/summary/merged \
  -p data/tumour_purity.txt \
  -o /Users/Nick/Desktop/SV_paper_20/svSupport \
  -cs


# svStitch
script_bin=/Users/Nick/iCloud/Desktop/script_test/svParser/script
for file in *_svSupport.txt
do
  python $script_bin/svStitch.py -i $file \
  -w 5500 \
  -o $(pwd)
done


unmappable=~/iCloud/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed

bash runParser -d ~/Desktop/SV_paper_20/svParser/wholegut/ \
  -o ~/Desktop/SV_paper_20/svParser/wholegut \
  -e $unmappable \
  -g /Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf \
  -r

out_dir=~/Desktop/SV_paper_20/svParser/wholegut/summary/merged
alleleFreqs=/Users/Nick/iCloud/Desktop/script_test/alleleFreqs
for f in ${out_dir}/*_reannotated_SVs.txt
do
  python $alleleFreqs/script/freqIn.py -v $f --config ${alleleFreqs}/data/samples.tsv
done

out_dir=~/Desktop/SV_paper_20/svParser/wholegut/summary/merged
for file in ${out_dir}/*_snv_added.txt
  do
  stem=$(basename "${file}")
  output_base=$(echo $stem | cut -d '_' -f 1)
  echo "Running for sample $output_base"
  if [ -f $out_dir/${output_base}_mh.txt ]
  then
    rm $out_dir/${output_base}_mh.txt
  fi
  python $script_bin/extract_vars.py -v $file -o ${out_dir}/${output_base}_mh.txt --write_genes --write_breakpoints
done

rm all_bps_merged.txt
rm all_genes_merged.txt
cat *_mh.txt > all_bps_merged.txt
cat *_hit_genes.txt > all_genes_merged.txt

python $script_bin/merge_files.py -e _snv_added.txt -d $out_dir -o all_samples_merged.txt

rm *_mh.txt
rm *_hit_genes.txt


# Need to annotate wholegut DELS
