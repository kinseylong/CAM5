# CAM5-IDR Project: Data and Analysis Overview

Code and data for the analysis of intrinsically disordered regions (IDRs) in the Arabidopsis thaliana CAM5 protein and related calmodulin proteins. The project integrates protein sequence data, structural predictions, and climate/genomic variation to explore the relationship between sequence features, disorder, and environmental adaptation.

## Data Overview

- **TAIR10/ArabidopsisUniprot.tsv**: UniProt annotations for Arabidopsis proteins, including AlphaFoldDB IDs, gene names, and sequences.
- **TAIR10/TAIR10proteins.fasta**: FASTA file of all Arabidopsis protein sequences.
- **TAIR10/IDRs.csv**: Output file listing predicted IDRs for each protein, including start/end positions, sequence, and AlphaFold confidence scores.
- **TAIR10/ArabidopsisReducedAlphabet.csv**: Protein sequences reduced to various simplified alphabets for comparative analysis.
- **GreneNetClimate/2029gaccession.csv**: Metadata for Arabidopsis accessions, including geographic and collector information.
- **GreneNetClimate/2029gclimate.csv**: Climate data (e.g., mean annual temperature) for each accession.
- **GreneNetClimate/chr2_11531967_11534358.recode.vcf**: Genomic variant data (VCF) for CAM5 region.
- **GreneNetClimate/pseudogenomes_CAM5.fa**: Pseudogenome sequences for CAM5 region across accessions.

## Notebooks and Scripts

- **ArabidopsisIDRpuller.py**: Script for batch extraction of IDRs from AlphaFoldDB.
- **CAM5_IDR.ipynb**: Main notebook for IDR extraction, reduced alphabet analysis, and Hamming distance calculations.
- **climateaccessions.ipynb**: Notebook for integrating climate, genomic, and sequence data, and for motif/statistical analyses.
- **IDR_distance_calc.py**: Script for calculating sequence distances in reduced alphabets.

## Analysis Pipeline

1. **Protein Sequence Extraction**
   - Extract all protein sequences from TAIR10 FASTA files and write to text for downstream processing.

2. **IDR Prediction**
   - Use AlphaFoldDB confidence categories to identify IDRs (regions of low confidence, D/L) in each protein.
   - Filter IDRs by minimum length and store start/end positions, sequences, and confidence scores in `IDRs.csv`.

3. **Reduced Alphabet Analysis**
   - Map protein sequences to reduced alphabets (4, 8, 10, 12, 15, 18 groups) to study sequence similarity and diversity.
   - Compute Hamming distances between CAM5 IDR and other IDRs in reduced alphabets.

4. **Disorder Prediction**
   - Use `metapredict` to calculate percent disorder for each protein or IDR.

5. **Climate and Genomic Variation Integration**
   - Join accession metadata, climate data, and pseudogenome sequences.
   - Parse VCF files to extract SNPs in the CAM5 intron region and merge with accession data.
   - Analyze sequence motifs (e.g., splice sites, runs of C/T, GT/AG motifs) and their association with climate variables.

6. **Statistical and Visualization Analyses**
   - Calculate motif frequencies, longest runs, and sequence features for each accession.
   - Group accessions by climate percentile and compare sequence features.
   - Visualize relationships using seaborn and matplotlib (e.g., boxplots, catplots).

