# Workshop: Hands on cheat sheet

**Note**: If internet connection is slow, we can also distribute the example data via an USB stick (ask your instructor ;) ). 

**Example data**: We will work example data on the HPC:

```bash
/scratch/Tausch/2023-RKI-NGS-Workshop/data/ngs_kurs
```

## 1. Prepare environment

* **Note**: Bioinformatics tools are regulary updated and input parameters might change (use `--help` or `-h` to see the manual for a tool!)
* Install most of them into our environment
    * we will already install many tools that we will use over the next days!
 
```bash
#start an interactive bash session using the default ressources
srun --pty bash -i

#start an interactive bash session using 8 CPUs, 40GB RAM, 30GB HDD
srun --cpus-per-task=8 --mem=40GB --gres=local:30 --pty bash -i

#start an interactive bash session using 10 CPUs, 80GB RAM, 50GB HDD, 1GPU
srun --cpus-per-task=10 --mem=80GB --gres=local:50 --gpus=1 --pty bash -i

#IMPORTANT to free the blocked resources after our work is done close the interactive shell via:
exit
```

```bash
# go to workplace folder
cd /scratch/$USER/......
conda activate envs/workshop
# test
NanoPlot --help
flye --version
# Switch to a path on your system where you want to store your data and results
cd /scratch/$USER
# Create new folder
mkdir nanopore-workshop
cd nanopore-workshop
```

## 2. Copy data

```bash
cp -r /scratch/Tausch/2023-RKI-NGS-Workshop/data/ngs_kurs/<folder> /scratch/$USER/nanopore-workshop/data/
# double-check that everything is in place:
ls -lah
# all good? Let's move on to QC!
```

### 2.1 combine data (if needed) 

```bash
cd nanopore-workshop
cat <folder>/barcodeXX/*.fastq.gz > data/<file>.fastq.gz
ls -lah
```

## 3. QC ONT data

```bash
cd nanopore-workshop
NanoPlot -t 4 --fastq <file>.fastq.gz --title "Raw reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/raw
```

```bash
# Note: we use 1 kb as the minimum length cutoff as an example. For your "real" samples other parameters might be better. Do QC before. 
filtlong --min_length 1000 --keep_percent 90 \
    --target_bases 500000000 <file>.fastq.gz > <file>-filtered.fastq

# Check the quality:
NanoPlot -t 4 --fastq <file>-filtered.fastq --title "Filtered reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/clean
```

[Code](https://github.com/rrwick/Filtlong)


## 4. _De novo_ assembly (Flye)

```bash
# run the assembly, this will take a bit time
flye --nano-raw <file>-filtered.fastq -o flye_output_<file> -t <cores> --meta --genome-size 5M
# the final output genome assembly will be in flye_output/assembly.fasta
```

## 5. Mapping (minimap2)

```bash
minimap2 -ax map-ont flye_output_<file>/assembly.fasta <file>-filtered.fastq > <file>-mapping.sam
```
[Publication](https://doi.org/10.1093/bioinformatics/bty191) | [Code](https://github.com/lh3/minimap2)

###5.1 Visualization of the mapping (IGV)

```bash
# first, we need to convert the SAM file into a sorted BAM file to load it subsequently in IGV
samtools view -bS <file>-mapping.sam | samtools sort -@ 4 > <file>-mapping.sorted.bam  
samtools index <file>-mapping.sorted.bam

# start IGV browser and load the assembly (FASTA) and BAM file, inspect the output
igv &

# load mapping file as 'primary assembly'
# ->  <file>-mapping.sorted.bam

# load assembly file as 'Reference/consensus file'
# ->  flye_output_<file>/assembly.fasta
```

## 6. Assembly polishing 

### 6.1 First Racon

```bash
# run racon, as input you need the reads, the mapping file, and the assembly you want to polish
racon -t 4 <file>-filtered.fastq 2<file>-mapping.sam flye_output_<file>/assembly.fasta > <file>-consensus-racon.fasta

# map to new consensus
minimap2 -ax map-ont <file>-consensus-racon.fasta <file>-filtered.fastq > <file>-consensus-mapping.sam

# now look at it in tablet or IGV again
igv &

# load mapping file as 'primary assembly'
# ->  <file>-consensus-mapping.sorted.bam

# load assembly file as 'Reference/consensus file'
# ->  flye_output_<file>/assembly.fasta
```

### 6.2 Second Medaka


`Medaka` is not in your current `workshop` environment because it was conflicting with the other tools. That's why we need a separate Conda environment for `Medaka`:

* Make a new environment for `medaka` 
    * `medaka` might have many dependencies that conflict 
* an alternative to `conda` is `mamba`
    * `mamba` can be much faster in solving your environment, e.g. here for the tool `medaka`
    * thus, let us install `mamba` via `conda` and then install `medaka`

```bash
mamba create -y -p envs/medaka "medaka>=1.8.0"
conda activate envs/medaka
```

```bash
# Run Medaka
# ATTENTION: it is always good to assign an appropriate Medaka model -m based on 
# the performed basecalling! Here, we use some example model for the E. coli 
# data. Adjust that in the following Exercise! 
# If you are on the RKI HPC: due to restrictions it might be even difficult to run other Medaka models because 
# they need to be downloaded first. 
medaka_consensus -i <file>-filtered.fastq -d <file>-consensus-racon.fasta -o <file>-medaka -m r941_min_sup_g507 -t 4

# now look at it in tablet or IGV again
igv &

# compare alle three assemblies
```
[Code](https://github.com/nanoporetech/medaka)

