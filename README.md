# Circulation of mumps virus genotype G in Pakistan during the 2023 outbreak

Reproducible code and workflow that mirror the analysis in the Future Virology research article on the Islamabad 2023 mumps outbreak. This repository covers both partial SH gene analysis and whole genome analysis and follows the toolchain and model choices reported in the paper.

## Program summary
End to end analysis of MuV genotype G using mNGS and Sanger SH sequences. Steps match the study design and are fully scripted so reviewers and collaborators can reproduce results.

Components in one workflow:
- Inputs
  - Paired end FASTQ from buccal or throat swabs sequenced on Illumina MiSeq 2x150 for WGS.
  - Optional FASTA of Sanger SH amplicons for SH only analysis.
- Quality control and trimming
  - FastQC for initial QC.
  - Trimmomatic for adapter and quality trimming with parameters used in the study: ILLUMINACLIP:adapters-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50.
  - Picard MarkDuplicates to identify PCR duplicates.
- De novo assembly and contig validation
  - SPAdes assembly.
  - Contig QC table and N50.
  - BLAST remote against NCBI nt to find the closest MuV matches.
- Reference based mapping and masked consensus
  - Select best reference from BLAST hits or from a preferred list of genotype G accessions.
  - Map with BWA MEM then sort and index.
  - Depth profile then Bed masking of sites with depth lower than the configured threshold.
  - Call variants and create masked consensus with bcftools. Default minimum depth 10.
- SH gene analysis
  - Extract the 316 nt SH region from WGS consensus using alignment to a genotype G reference then build a SH alignment that includes Sanger SH sequences if provided.
  - IQ-TREE maximum likelihood tree for SH with K2+I and 1000 ultrafast bootstraps. Optionally use ModelFinder.
  - Mutation table for SH versus Sheffield genotype G reference and the Jeryl Lynn vaccine strain.
- Whole genome phylogeny
  - Build a whole genome alignment for all masked consensuses plus context genomes.
  - IQ-TREE tree for whole genome with GTR+G and 1000 ultrafast bootstraps. Optionally use ModelFinder.
- HN and F protein mutation analysis
  - Extract HN and F coding regions using alignment to a genotype G reference then translate and compare amino acids versus Jeryl Lynn minor and major components.
  - Report per position differences and a concise summary by gene.
- Outputs
  - Consensus FASTA per sample and a combined FASTA.
  - SH and whole genome alignments and Newick trees.
  - Mutation reports for SH, HN and F.
  - QC, BLAST, and reference selection tables.

Methods and model choices reflect the study. SH trees use K2+I. Whole genome trees use GTR+G. Minimum depth for consensus is 10. Tool versions match the paper where possible. See the paper for details. 

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake for the full pipeline
- BLAST+ with remote enabled and Entrez Direct

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One command end to end run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate mumps-2023-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Minimal example:
```yaml
pairs:
  - sample: MuV_001
    r1: data-private/MuV_001_R1.fastq.gz
    r2: data-private/MuV_001_R2.fastq.gz

# Optional Sanger SH sequences
sh_fasta: data-sh/sh_sequences.fasta   # comment this line if not available

# Reference preferences and context sets
reference_preference:
  genotype_g_pref:
    - ON148331.1   # Sheffield
    - JN012242.1   # Iowa
    - LC633315.1   # Yokohama 2010
    - KY006858.1   # Ontario 2010
    - EU370207.1   # Croatia 2005
    - MW261742.1   # Utrecht 2010
    - AF280799.1   # WHO G
context_wg_acc:
  - ON148331.1
  - JN012242.1
  - JX287385.1
  - LC633315.1
  - KY006858.1
  - EU370207.1
  - MW261742.1
  - LC633320.1
  - LC685516.1
  - LC685525.1

context_sh_acc:
  - AF280799.1   # WHO reference for genotype G
  - FN431985.1   # Jeryl Lynn
  - KY006858.1
  - KF738114.1
  - KY604739.2
  - MT506353.1
  - KX609907.1
  - KX609928.1
  - KX609935.1

params:
  threads: 4
  trim_adapters: env/TruSeq3-PE.fa
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model_sh: K2+I
  iqtree_model_wg: GTR+G
  bootstrap: 1000
  max_blast_hits: 50
  use_model_finder: false
```
Notes
- If `use_model_finder` is set to true the pipeline adds `-m MFP` to IQ-TREE.
- SH extraction script aligns each sample consensus to a genotype G reference to locate the 316 nt SH segment and then writes the sample SH segment to FASTA.

## How to cite
- Paper: Umair M, Haider SA, Jamal Z, Hakim R, Farooq A, Salman M. Circulation of Mumps virus genotype G in Pakistan during 2023 outbreak. Future Virology. 2024. 18(18):1137-1149. https://doi.org/10.2217/fvl-2023-0145
- Software: Haider SA. Mumps genotype G genomic analysis for the 2023 Islamabad outbreak. Version 1.0.0. {year}. GitHub repository. Include commit hash when available.

## References
- Andrews S. 2010. FastQC. Babraham Bioinformatics. 
- Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30:2114-2120.
- Broad Institute. 2019. Picard Toolkit. GitHub repository.
- Bankevich A, et al. 2012. SPAdes: a new genome assembly algorithm and its applications to single cell sequencing. J Comput Biol 19:455-477.
- Li H. 2013. Aligning sequence reads with BWA-MEM. arXiv:1303.3997.
- Li H, et al. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25:2078-2079.
- Danecek P, et al. 2021. Twelve years of SAMtools and BCFtools. GigaScience 10:giab008.
- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7. Mol Biol Evol 30:772-780.
- Minh BQ, et al. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference. Mol Biol Evol 37:1530-1534.
- Larkin MA, et al. 2007. Clustal W and Clustal X version 2.0. Bioinformatics 23:2947-2948.
- Shen W, Xiong J. 2016. SeqKit: a cross platform toolkit for FASTA and FASTQ. PLoS One 11:e0163962.
- Camacho C, et al. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10:421.
- Kans J. 2010 and later. Entrez Programming Utilities Help. NCBI.
- Cock PJ, et al. 2009. Biopython: freely available Python tools for computational molecular biology. Bioinformatics 25:1422-1423.
- Köster J, Rahmann S. 2012. Snakemake: a scalable bioinformatics workflow engine. Bioinformatics 28:2520-2522.


## Contributing
Contributions are welcome. Please read `CONTRIBUTING.md` for workflow and style. Open an issue for questions. Do not commit patient or restricted data.

## License
MIT. See `LICENSE` for details.
