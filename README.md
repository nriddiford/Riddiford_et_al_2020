# SV Pipeline

## Variant calling

### 1. Align reads to combined genomes and tag TE-mapped reads
* Run `bwa_submit.sh`
    * Assumes files are named : '${group}${sample}.info.fq'
    * Assumes tumour/normal are sequentially named (tumour: R1, normal R2, tumour R3, normal R4)

### 2. Call CNVs with control-FREEC and CNV-Seq
#### 2.1 control-FREEC (v11.0)

* Run `bash run_controlFreec.sh <group>`. Uses default settings, & selects for sig cnvs

#### 2.1 CNV-Seq

* Run `bash cnvSeq_counts_submit.sh <group>` to produce per-sample hits file (read counts per window)
* Run `bash cnvSeq_submit.sh <group>` which runs CNV-Seq on tumour-normal pairs with three different window sizes (500, 10000, 50000) with `--global-normalization` option


### 3. Call structural variants
#### 3.1 novoBreak (v1.1 (r20151007))
* Run `bash novobreak_submit.sh <group>` to run `run_novoBreak.sh` on tumour-normal pairs using default settings

#### 3.1 lumpy (v0.2.13)
* Run `bash lumpy_full.prep.sh <group>` to generate per-sample discordant and split reads + insert size distribution
* Run `bash lumpy_full.sh <group>` to run lumpy on tumour-normal pairs and generate a list of normals in group `${group}.PON` for svTyper. Specifies unmappable genome as `-e`

#### 3.2 svTyper (svtools 0.3.0)
* Run `bash svTyper.sh <group>` to create a PON and annotate these onto lumpy calls `${tumour_sample}.lumpy.gt_all.vcf`

#### 3.3 Delly (v0.7.8)
* Run `bash delly_call.sh <group>` to run delly. Requires PON from above step. Outputs `${tumour_sample}.delly.vcf`


## Variant filtering & annotation

### 1. Organise files for filtering and install dependencies
```{bash}
mkdir -p data
cd data
mkdir -p lumpy delly novobreak freec cnv
```
* Move all SV calls into corresponding directory in `data/`
* Organise all CNV-Seq files (`*.cnv`) generated with a small window (we use 500bp) into a single directory
* Install dependencies following [instructions on github](https://github.com/nriddiford/svParser#installation)

### 2. Filter, merge and annotate somatic variants (-fmas) using [svParser](https://github.com/nriddiford/svParser)

#### Need to include the freec filtering step: `bash ~/Desktop/parserTest/script/freecFilt.sh *_sig_cnvs.txt`
```{bash}
cd svParser
dir=$(dirname "$0")
script_bin="$dir/script"
out_dir=filtered
```
* Run `runParser` wrapper for `svParser`:
```{bash}
  bash runParser \
  -d data/ \
  -o $out_dir \
  -fmas \
  -c /path/to/cnv-seq/files/
```

* Final out put is a combined, per-sample annotated breakpoints file `${tumour_sample}_annotated_SVs.txt`

### 3. Refine calls, filter, genotype and calculate allele frequency
* Install [svSupport](https://github.com/nriddiford/svSupport)
* To run svSupport for each `${tumour_sample}_annotated_SVs.txt` file
* Run `bash scripts/run_all.sh -v ~/Desktop/final_analysis/filtered/summary/merged/ -o /Users/Nick_curie/local_data/svSupport/ -cs`

### 4. Stitch together complex variants
```{bash}
svSupport_dir=/Users/Nick_curie/Desktop/script_test/svSupport
for file in ${svSupport_dir}/*_svSupport.txt
do
  python $script_bin/svStitch.py -i $file \
  -w 5500 \
  -o $out_dir/summary/merged
done
```
* Output: `${tumour_sample}_stitched.txt`


### 5. Manually inspect calls, annotate uncalled events and mark false positives

```{bash}
bash runParser -d ~/Desktop/final_analysis/data/ \
-o $out_dir \
-r
```
* Input: `${tumour_sample}_stitched.txt`
* Output: `${tumour_sample}_reannotated_SVs.txt`

### 6. Optional - annotate mechanisms called by splitvision

#### Extract bedpe for splitvision
```{bash}
for file in *_stitched.txt
  do
  stem=$(basename "${file}")
  output_base=$(echo $stem | cut -d '_' -f 1)
  echo "Running for sample $output_base"

  python ~/Desktop/script_test/svParser/script/parser2bedpe.py -i $file
done
```

* `cd $out_dir/summary/merged`
```{bash}
bash ~/Desktop/parserTest/script/run_mhannotate.sh -axm -d . -o .
```

### 8. Optional - annotate and filter CNVs with SNP allele frequencies
alleleFreqs=/Users/Nick_curie/Desktop/script_test/alleleFreqs
```{bash}
for f in ~/Desktop/final_analysis/filtered/summary/merged/*_mechanisms.txt
do
  python $alleleFreqs/script/freqIn.py -v $f
done
```

### 9. Extract all vars and merge into R-friendly files

```{bash}
out_dir=~/Desktop/final_analysis/filtered/summary/merged
for file in ${out_dir}/*_snv_added.txt
  do
  stem=$(basename "${file}")
  output_base=$(echo $stem | cut -d '_' -f 1)
  echo "Running for sample $output_base"
  if [ -f $out_dir/${output_base}_mh.txt ]
  then
    rm $out_dir/${output_base}_mh.txt
  fi
  python ~/Desktop/script_test/svParser/script/extract_vars.py -v $file -o ${out_dir}/${output_base}_mh.txt --write_genes --write_breakpoints
done

rm all_bps_merged.txt
rm all_genes_merged.txt
cat *_mh.txt > all_bps_merged.txt
cat *_hit_genes.txt > all_genes_merged.txt

python ~/Desktop/script_test/svParser/script/merge_files.py -e _snv_added.txt -d $out_dir -o all_samples_merged.txt

rm *_mh.txt
rm *_hit_genes.txt
```


# SNV Pipeline

### 1. Align reads to dmel6 and remove clipped reads (~ 10 hours)
* Run `bwa_SNV_submit.sh`
* Dependencies: [samclip](https://github.com/tseemann/samclip)
    * Assumes files are named : '${group}${sample}.${whatever}.fq'
    * Assumes tumour/normal are sequentially named (tumour: R1, normal R2, tumour R3, normal R4)

### 2. Make Mutect2 Panel Of Normals (~ 20 mins)
* Run Mutect2 in tumour-only mode (`run_mutect2_pon.pbs`)
* Run `mutect2_combineNormals.sh` to combine all into `panel_of_normals.vcf.gz`

### 3. Call SNVs

#### 3.1 Run Mutect2 (GATK 4.1.0) in tumour/normal mode
* Run Mutect2 (`run_mutect2_tn.pbs`)
* Uses default `FilterMutectCalls` and outputs SNVs/INDELS that PASS filters

#### 3.2 Run Varscan2 (v2.4.3) in tumour/normal mode

##### 3.2.1 Create merged tumour/normal pileups
* Run `bash makePileup.sh <group>`

##### 3.2.2 Run Varscan
* Run `bash Varscan2_pileup.submit.sh <group>`

#### 3.3 Run Freebayes (v1.2.0-dirty) in tumour/normal mode
* Run `bash freebayes_submit.sh <group>`

#### 3.4 Run Strelka (v2.9.10) in tumour/normal mode
* Run `bash strelka.sh <group>`

#### 3.5 Run Somatic Sniper (v1.0.5.0) in tumour/normal mode
* Run `bash sniper.sh <group>`


### 4. Merge calls

#### 4.1 Run Somatic Seq
* Run `bash somaticSeq.sh <group>`

#### 4.2 Run combinevcf.py to add Freebayes calls
* `cd /Users/Nick_curie/Desktop/script_test/mutationProfiles`
* Run `bash script/run_combine_vcf.sh`
* Outputs `${tumour_name}_merged.vcf`

#### 4.3 Run filtervcf.py to filter against PON & depth < 20
* `cd /Users/Nick_curie/Desktop/script_test/mutationProfiles`
* Run `python script/filtervcf.py -d .`
* Outputs `${tumour_name}_consensus_filt.vcf`

### 5. Annotate calls
#### 5.1 Annotate with SnpEff

##### Move all _consensus_filt.vcf into CombinedVCF/

`bash script/run_snpEff.sh -c`
`bash script/run_snpEff.sh -i`

#### Move all files into mutationProfiles data/raw

`conda activate svParser`
`bash run_trinucs.sh -ea`

# LOH Pipeline

### 1. Run SNV pipeline steps 1, 3.2 & 3.3

### 2.
