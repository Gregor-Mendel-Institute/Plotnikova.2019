# Plotnikova.2019
Custom software used to align sRNA-seq data to annotated mature miRNAs and generate figures for Plotnikova et al. (2019) MicroRNA Dynamics and Functions During Arabidopsis Embryogenesis, unpublished.


General Description (from manuscript):
Cutadapt (Martin, 2011) was used to trim adapter sequences from sRNA-seq reads and 18–30-base sequences that contained an adapter were retained. The trimmed sequences were aligned to the Arabidopsis thaliana TAIR10 genome (Lamesch et al., 2012) with STAR (Dobin et al., 2013) requiring no mismatches and allowing ≤100 multiple end-to-end alignments. Resulting SAM files were then processed with the readmapIO.py script to re-assign multimappers with a “rich-get-richer” algorithm as previously described (Schon et al., 2018). Output bedFiles were sorted, condensed and normalized for total genome-matching reads. The BEDtools map function (Quinlan and Hall, 2010) was then used to quantify the number of reads mapping to the same strand and overlapping ≥80% of mature miRNAs as annotated in TAIR10 and miRBase21 (Kozomara et al., 2019). Statistical analyses and associated figures were generated with the R statistical computing package (R Core Team, 2018).


Detailed Descriptions for Each Step:

sRNA-seq Alignments
  A PBS script (sRNA.align.git.sh) is included that describes the mapping procedure in eight steps. You will need the following software:
      1. Cutadapt (cutadapt/1.18-foss-2018b-Python-2.7.15)
      2. STAR (rna-star/2.5.2a-foss-2016a)
      3. readmapIO.py (included, from Schon et al. (2018) Genome Research)
      4. fasta_utils.py (included, from Schon et al. (2018) Genome Research)
      4. BEDtools (BEDTools/2.27.1-GCCcore-6.4.0)
      5. bed_collapse.py (included, this study)
  And the following files: 
      1. FASTQ file of sRNA-seq reads (test.fastq included for testing)
      2. FASTA file of reference genome (not included but can be downloaded from ftp://ftp.ensemblgenomes.org/pub/plants/release-       44/fasta/arabidopsis_thaliana/dna/) 
      3. BEDfile of miRNA loci (miRNA_mature.bed included)
  Executing the series of progams in sRNA.align.git.sh should generate a test.quant.on.miRNA_mature.tsv file with miRNAs as rows and the     number of reads for across various sRNA-seq read sizes as columns. Values were normalized for either the number of non-unique alignments   (i.e. raw) or the number of non-unique alignments and reads per million genome-matching reads (i.e. norm). This procedure was performed   for the 33 sRNA-seq datasets: 4 dilution series, 24 wild-type embryos, 3 dcl1-5 globular embryos, 1 leaf and 1 floral buds      
 
Figure 1 and Supplemental Figure 1 (see Figure_1_S1.R)
  

