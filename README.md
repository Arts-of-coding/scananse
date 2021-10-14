# single_cell_ANANSE_pipeline

# Introduction
This repository contains all code and command used to analyse a publicly available dataset of the human corneal atlas.

# Installing Miniconda.
Install Miniconda within an existing Linux environment (Ubuntu).

Retrieve the latest Miniconda.
```console
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Generate the script to be executable.
```console
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
```

Run the executable miniconda script.
```console
$ ./Miniconda3-latest-Linux-x86_64.sh
```

# Installing the right channels for downstream analysis in Miniconda.
Install the nessesary channels if Bioconda has not been used before. 
```console
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

# Preprocessing files for single cell integration from public 10X data
You can download publicly fastq data once conda is installed, the samples.tsv and the config.yaml are correctly filled in. Additionally, you need to have the right activated seq2science environment.

```console
nice -30 seq2science run download-fastq -j 8 -k -c /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/seq2science_config/config.yaml
```
Next, the SRA files need to be split, because seq2science can't deal properly with single cell data. This can be accomplished by using sra-tools. In the correct activated sra-tools conda environment run the code below.

```console
nice -30 fastq-dump --split-files --outdir /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/fastq/ /ceph/rimlsfnwi/data/moldevbio/zhou/jsmits/iPSC_LSC_s2s/data_Lako/Lako_data/sra/SRR12386341/SRR12386341/SRR12386341.sra
```

Lastly, the fastq files need to have the correct naming conventions to be used by the tool cellranger. For installing cellranger, please take a look at: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation.

# Generating nessesary scRNA-seq files for SEURAT
With the correct input files, cellranger count can be run.

BASH example code:
```console
cellranger count --id=lako341 --transcriptome=/.... --fastqs=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/fastq/ --sample=SRR12386341 --expected-cells=... --localcores=... --localmem=...
```

# Quality control and visualization for scRNA-seq data with Seurat
To use Seurat you need three files (generated with cellranger count: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.

# Generating nessesary scATAC-seq files for snapATAC
With a correctly installed cellranger-atac environment (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation), run cellranger-atac count.

BASH example code:
```console
nice -30 cellranger-atac count --id=ataclako4 --reference=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/ --fastqs=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/fastq/ --sample=SR12386372 --localcores=24 --localmem=100
```
Go to the output directory above and generate the snap file with samtools:
```console
samtools view possorted_bam.bam
```
Make the header file:
```console
samtools view possorted_bam.bam -H > ataclako4_possorted.header.sam
```
Make the unsorted.snap.bam:
```console
cat <( cat ataclako4_possorted.header.sam ) \
<( samtools view possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -bS - > ataclako4_possorted.snap.bam
```
View the newly created file:

```console
samtools view ataclako4_possorted.snap.bam | cut -f 1 | head 
```
Generate the possorted snap.bam file:
```console
nice -30 samtools sort -n -@ 10 -m 1G ataclako4_possorted.snap.bam -o ataclako4.snap.nsrt.bam
```
Download the good reference genome file in the current environment and directory cellranger atac:
```console
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```

From the cellranger atac directory or a directory of choice, activate the snaptools environment and make the snap file:
```console
nice -30 snaptools snap-pre  \
	--input-file=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako4/outs/ataclako4.snap.nsrt.bam  \
	--output-snap=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako4.snap  \
	--genome-name=hg38  \
	--genome-size=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/hg38.chrom.sizes  \
	--min-mapq=30  \
	--min-flen=50  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=20000  \
	--min-cov=500  \
	--verbose=True
```

Add bins to the final snap object:
```console
nice -30 snaptools snap-add-bmat \
    --snap-file=ataclako4.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True
```

# Adding pmat to snap (nessessary for plotting the differentially accessible regions within snapATAC):
```console
snaptools snap-add-pmat --snap-file /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako2.snap --peak-file /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210604_peaks_combined.bed
```

# Generating bam files from distinct cell populations
From the R file you retrieve several files where all barcodes present in each cell population (for each dataset). From these generated files you can subslice the original bam files and merge them together. First go to the directory where all of these metadata.csv files are stored (see R file) and set this as your directory.

Remove unwanted X column header and '"' characters
```console
sed -i '1d' *_metadata.csv
sed -i ' s/b//g' *_metadata.csv
```
Next install cellranger-dna (see docs Cellranger). Specify for each dataset the barcodes.csv (see template file). Afterwards export directory (see notes below) and run the subslice command in Bash. Do this for each dataset.

```console
nice -30 cellranger-dna bamslice --id testLNPCatac --csv=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210615scRNA_integration/barcodes.csv --bam=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako2/outs/possorted_bam.bam
```

Activate your conda samtools environment. Now you are ready to merge the bamfiles for each cell population.

```console
nice -30 samtools merge mergedCjS234.bam /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210615scRNA_integration/lako2tracks/outs/subsets/CjS.bam /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210616tracks_integration_all/lako3tracks2/outs/subsets/CjS.bam /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210616tracks_integration_all/lako4tracksV/outs/subsets/CjS.bam
```

# Generating narrowPeak files for the combine peaks function
First, install MACS2 peak caller via conda and activate the environment. Then run the command below to generate peak files for each population.

```console
for bam in *.bam
do
macs2 callpeak -t $bam --name $bam --nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR
done
```

# Making bigwig tracks from bdg files
It is possible to also make bigwig tracks from the output of MACS2. Excluding nonsense chromosomes which give errors when compared to the reference genome file.
```console
grep -Ev 'GL000|KI|chrM' $bdg_file > $ bdg_file2
```
Removing unwanted characters from the bed files
```console
nice -30 sed -i ' s/b//g' $bdg_file2 && sed -i " s/[']//g" $bdg_file2
```
Generating the tracks in a bigwig environment
```console
for bdg_file in *bdg; do LC_COLLATE=C sort -k1,1 -k2,2n $bdg_file > conversion_file.bdg ; bedGraphToBigWig conversion_file.bdg /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/hg38.chrom.sizes $(basename -s .bdg $bdg_file).bw; rm conversion_file.bdg; done
```

# Gimme Maelstrom analysis pipeline in bash and Conda
Exclude unwanted chromosomes in bash.
```console
grep -Ev 'GL000|KI|chrM' *.narrowPeak
```

Next input the .narrowPeak files and the corresponding .bam files in the files.tsv, like specified below:
```console
$ head files.tsv

filename        file_location   data_type       cell_type       replicates
#New filepaths
genome_gtf      /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/genome/hg38/hg38.annotation.gtf  gtf_file
genome_path_size        /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/genome/hg38/hg38.fa.sizes        fa_sizes
###peakfiles
Ves_ATAC_peaks  /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_4cells/4cells/outs/subsets/LiCo.bam_peaks.narrowPeak  ATAC_peak LiCo
StC_ATAC_peaks  /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_4cells/4cells/outs/subsets/StCSC.bam_peaks.narrowPeak ATAC_peak StCSC
## ATACseq BAM files;
LiCo_ATAC_BAM   /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_4cells/4cells/outs/subsets/LiCo.bam   ATAC_BAM        LiCo
StCSC_ATAC_BAM  /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_4cells/4cells/outs/subsets/StCSC.bam  ATAC_BAM        StCSC
```
If you want to calculate the significant peaks (for Gimme Maelstrom analysis), then you also need to specify replicate files as well in the files.tsv. The replicates were generated within "Workflow_snapATAC_integration_V5.R".
```console
$ tail files.tsv
## ATACseq BAM replicate files
CSB_ATAC_BAM    /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_repex/reps/outs/subsets/CSB1.bam      ATAC_BAM        CSBrep1 CSBrep1
CSB_ATAC_BAM    /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_repex/reps/outs/subsets/CSB2.bam      ATAC_BAM        CSBrep2 CSBrep2
CSB_ATAC_BAM    /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/split_bam_repex/reps/outs/subsets/CSB3.bam      ATAC_BAM        CSBrep3 CSBrep3
```

To determine motif analysis on significant peaks, run the Python script "Combine_peaks_V2.ipynb" with measured replicates. Next, run the R script "significance_calc_peaks_gimme.R" for determining the significant peaks.

The joined coverage table file can be quantile normalized and subselected for significant regions with the Python script "TPM_Quantile_script_V3.ipynb" . This file can be used as input for Gimme Mealstrom.
```console
$ head quantilepeaks.txt
loc	Ves	StC	MEC	LPCs	LNPCs	IC	FCECs	CjS	CSSCs	CSB	CB
chr1:911401-911601	1.0851503710073496	0.769986796647647	0.963644488152334	2.55757675487907	2.901962254588777	2.5544817899948575	2.438817238123795	2.438817238123795	1.3817937286844038	2.9084297593729525	2.714010079063348
chr1:967906-968106	1.0851503710073496	1.7798851919799248	0.963644488152334	1.673191906190539	2.901962254588777	1.007465278910087	
```

Install and activate your Gimme motifs (Gimme3.yaml) environment in conda and run Gimme Maelstrom.
```console
$ gimme maelstrom quantilepeaks.txt /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/genome/genomehg38gimme/hg38/hg38.fa /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/motif_analysis/all_differential/ -N 6 -F 0.8
```

# Combined analysis of motifs and gene expression
With the script "significance_calc_RNA_DE_ESC.R" you can generate a Z-score matrix of gene expression. You can determine which transcription factors binding to your motifs are also expressed in your single cell populations. Run the R script "complex_heatmap_motifs_V3.R" to generate heatmaps with gene expression and quantile normalization.

# Generating the nessesary files for ANANSE and running ANANSE.
Lastly, ananse can be run with the python script "Run_ananse_V4_JA.ipynb". The files can be further analysed with the R script "ANANSE_analysis_V4.R".

# Notes
Note for setting the correct paths in Cellranger within the terminal:
```console
export PATH=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellranger/cellranger-6.0.1:$PATH
```

```console
export PATH=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0:$PATH
```

```console
export PATH=/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangerdna/cellranger-dna-1.1.0:$PATH
```
