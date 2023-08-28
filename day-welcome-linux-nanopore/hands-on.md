# Workshop: Linux re-cap and Nanopore QC

**Note**: If internet connection is slow, we can also distribute the example data via an USB stick (ask your instructor ;) ). 

**Example data**: We can work with example data on the HPC:

```bash
/scratch/Tausch/2023-RKI-NGS-Workshop/data/previously-nanopore-sequenced/
```

Please Note: If possible, you can instead also work with your own data that you generated yesterday. We can discuss this live during the workshop.

## Hands-on

### Create a folder for the hands-on work

Below are just example folder names, you can also adjust them and use other folder names! Assuming you are on a Linux system on a local machine (laptop, workstation):

```sh
# Switch to a path on your system where you want to store your data and results (you should be already on this path)
cd /scratch/$USER
# Create new folder
mkdir nanopore-workshop
cd nanopore-workshop
mkdir data
```

### Install and use analysis tools

* **Note**: Bioinformatics tools are regulary updated and input parameters might change (use `--help` or `-h` to see the manual for a tool!)
* Install most of them into our environment
    * we will already install many tools that we will use over the next days!

```bash
cd /scratch/$USER/nanopore-workshop
mkdir envs
mamba create -y -p envs/workshop fastqc nanoplot filtlong flye bandage minimap2 tablet racon samtools igv
conda activate envs/workshop
# test
NanoPlot --help
flye --version
```

__Reminder: You can also install specific versions of a tool!__
* important for full reproducibility
* e.g. `mamba install flye==2.9.0`
* per default, `mamba` will try to install the newest tool version based on your configured channels and system architecture and dependencies to other tools

### Get some example long-read data 

Get some example data. Either from your own sequencing run or, for example,

```bash
cp -r /scratch/Tausch/2023-RKI-NGS-Workshop/data/previously-nanopore-sequenced/230217_GI1-4_Run23-047 /scratch/$USER/nanopore-workshop/data
# double-check that everything is in place:
ls -lah data/
# all good? Let's move on to QC!
```

### ONT basecalled fastq output 

Nanopore data can be basecalled during the run. You can also run basecalling yourself, but it is time consuming and we'll skip it here. The output is sorted into barcodes and in every folder are a list of compressed fastq files. These files need to be combined into one file for further use.

```bash
# make sure that you are located in your workshop folder, if not cd
cat data/230217_GI1-4_Run23-047/fastq_pass/barcode01/*.fastq.gz > data/barcode01.fastq.gz
ls -lah data/
```

**Attention:** Before doing any calculations it makes sense to start an interactive compute session on the HPC! See [Linux Crash Course](../linux.md). If you do so, you need to `conda activate envs/workshop` again!

### Quality control (NanoPlot)

```bash
cd /scratch/$USER/nanopore-workshop
NanoPlot -t 4 --fastq data/barcode01.fastq.gz --title "Raw reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/raw
```
[Publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939) | [Code](https://github.com/wdecoster/NanoPlot)

**Note**: The `\` at the end of a line is only for convenience to write a long command into several lines. It tells the command-line that all lines still belong together although they are separated by "enter" keys. However, if you type all of the command, i.e., paths etc, in one line do not copy/type the backslash at the end of the lines.

### Read filtering (Filtlong)

```bash
# Note: we use 1 kb as the minimum length cutoff as an example. For your "real" samples other parameters might be better. Do QC before! 
filtlong --min_length 1000 --keep_percent 90 \
    --target_bases 500000000 data/barcode01.fastq.gz > barcode01-filtered.fastq

# Check the quality again:
NanoPlot -t 4 --fastq barcode01-filtered.fastq --title "Filtered reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/clean
```
[Code](https://github.com/rrwick/Filtlong)


## Excercise

1) Investigate the content of the FASTQ file. What are the first four lines telling you? What do you need to do to make the content of the file "human readable"? 
2) Try another tool, `FastQC`, on the example FASTQ. Inspect the output. What is the GC content? 

## Bonus 1

1) Use `PycoQC` to generate qc plots for the data set. Install `PycoQC` via Mamba or use an available environment. In difference to `NanoPlot`, `PycoQC` needs as input a file called `sequencing_summary.txt` or similar. This is provided after the basecalling alongside with the FASTQ files.

**Note on installing `PycoQC`**: On my system it was a pain to install `PycoQC`. I finally managed using mamba:

```bash
mamba create -y -p envs/pycoqc pycoqc
```

If that doesn't work you might need to:

* creating a new conda environment and installing `mamba` and Python version 3.7 into it
* activating the new environment
* then installing `PycoQC` and explicitly defining the newest version 2.5.2 (as of 2023-08-07)

Older versions might not work correctly with the input FAST5 data! Maybe you have better luck in installing a newer version of `PycoQC` w/o the hussle... 
