# Workshop: Assembly polishing and variant calling

## Hands-on

### Polishing

Polish the assembly you produced with `flye`. Use the filtered long-read data. Do this in two steps: using `racon` followed by `medaka`. This will also involve mapping again the long reads to your calculated assemblies! 

#### Mapping (minimap2)

**You already did this mapping to look at the SAM/BAM file in IGV and/or tablet! If you still have the files, you dont need to redo the following steps!**

```bash
# map
minimap2 -ax map-ont flye_output/assembly.fasta barcode01-filtered.fastq > barcode01-mapping.sam
# first, we need to convert the SAM file into a sorted BAM file to load it subsequently in IGV
samtools view -bS barcode01-mapping.sam | samtools sort -@ 4 > barcode01-mapping.sorted.bam  
samtools index barcode01-mapping.sorted.bam
```

#### Assembly polishing (Racon)

```bash
# run racon, as input you need the reads, the mapping file, and the assembly you want to polish
racon -t 4 barcode01-filtered.fastq barcode01-mapping.sam flye_output/assembly.fasta > barcode01-consensus-racon.fasta

# map to new consensus
minimap2 -t 4 -ax map-ont barcode01-consensus-racon.fasta barcode01-filtered.fastq > barcode01-consensus-mapping.sam

# now look at it in tablet or IGV again. For IGV you have to convert to BAM again and index the mapping file!
igv &

tablet &
# load mapping file as 'primary assembly'
# ->  barcode01-consensus-mapping.sam

# load assembly file as 'Reference/consensus file'
# ->  flye_output/assembly.fasta


```
[Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411768/) | [Code](https://github.com/isovic/racon)

* a common practice is 2x `racon` polishing followed by 1x `Medaka` (see below)
* for practice, 1x `racon` is fine
* with the newest R10 chemistry and recent basecalling models, it seems that `racon` is also not necessary anymore and people switch to only polish via `Medaka`
* with improvements in ONT accuracy, it might be even the case that long-read polishing is becoming more and more obsolete!

### Assembly polishing and final consensus (Medaka)

`Medaka` is not in your current `workshop` environment because it was conflicting with the other tools. That's why we need a separate Conda environment for `Medaka`:

* Make a new environment for `medaka` 
    * `medaka` might have many dependencies that conflict 
* an alternative to `conda` is `mamba`
    * `mamba` can be much faster in solving your environment, e.g. here for the tool `medaka`
    * thus, let us install `mamba` via `conda` and then install `medaka`

```bash
cd /scratch/$USER/nanopore-workshop
mamba create -y -p envs/medaka "medaka>=1.8.0"
conda activate envs/medaka
```

```bash
# Run Medaka
# ATTENTION: it is always good to assign an appropriate Medaka model -m based on 
# the performed basecalling! Here, we use some example model for a R10.4.1 run with 260 bp/s speed and SUP basecalling. 
# This might to be adjusted based on your data! 
# If you are on the RKI HPC: due to restrictions it might be even difficult to run other Medaka models because 
# they need to be downloaded first. 
medaka_consensus -i barcode01-filtered.fastq -d barcode01-consensus-racon.fasta -o barcode01-medaka -m r1041_e82_260bps_sup_v4.0.0 -t 4

# Exercise: look at it in tablet or IGV
# Hint: first need a mapping to the new consensus again to generate the SAM/BAM file!
```
[Code](https://github.com/nanoporetech/medaka)

**Note** that you should usually change the model parameter (`-m`) to whatever is most appropriate for your basecalling. Also note that `medaka_consensus` is not the same thing as `medaka consensus` (underscore vs space) - the former is a convenience script which does the entire process (including read mapping) while the latter is a subcommand of Medaka which only does the polishing step. (thx to [Ryan Wick for this explanation](https://github.com/rrwick/Trycycler/wiki/Polishing-after-Trycycler)).


## Exercise

Now, use another sample, for example one that you sequenced the first day. Remember, how basecalling was done and with which parameters. 

Perform QC, filter the reads, and then calculate an `flye` assembly. Now polish the genome using `racon` and `medaka`. Try to chose an approriate `medaka` model. You can use the following command to list `medaka` models:

```bash
medaka tools list_models | grep -v Default
```

Did the per-base quality improve? Annotate genes (e.g. using `Prokka` or `Bakta`)! How many genes (CDS) do you find now in comparison to the _de novo_ assembly without any polishing? 

Now, we want to call variants for your Nanopore sample in comparison to a reference sequence (**not** using the _de novo_ assembly you calculated!).

* To do so, download a reference genome for your sample from [NCBI](https://www.ncbi.nlm.nih.gov/genome/)
* Map the Nanopore reads of your sample against the reference genome
* use `Medaka` now for variant calling and not directly for consensus calculation, here are some hints (that you can/must adjust! Check also https://github.com/nanoporetech/medaka):

```sh
# Call variants with Medaka
#__Important__: Always use the matching `medaka` model based on how the `guppy` basecalling was done! You can check which `medaka` models are available via:
medaka tools list_models | grep -v Default

# first, use the `medaka consensus` command similar to before
# for the Salmonelle ONT data from 2019 MinION device was used and the FAST model!
medaka consensus <SORTED-BAM-FILE> --model <MODEL> --threads 4 output.consensus.hdf

# actually call the variants
medaka snp <REFERENCE-GENOME> output.consensus.hdf medaka.snp.vcf --verbose
```

* inspect the resulting [VCF file](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/), read about the format and its structure
* how many variants do you find? 
* can you find the called variants that you see in the VCF file also in a genome browser, when you load the mapping file? Do you see any differences/problems?

## Bonus 1

Try different allele frequency cutoffs with `Medaka` for the variant calling: `--threshold 0.1` How do the results change?

## Bonus 2

If you also have short-read Illumina data corresponding to your Nanopore data. Use the high-accuracy short-read data to polish your _best_ long-read assembly again. Use `Polypolish` for that. Install `Polypolish` and [familiarize yourself with the tool](https://github.com/rrwick/Polypolish/wiki). For the necessary mapping of the short reads to the assembly you should use `bwa`. Check the `Polypolish` manual!  

### Quick start
(see https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish for more details!)

```sh
bwa index draft.fasta
bwa mem -t 16 -a draft.fasta reads_1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a draft.fasta reads_2.fastq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish draft.fasta filtered_1.sam filtered_2.sam > polished.fasta
```
