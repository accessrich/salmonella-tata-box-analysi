# Salmonella TATA Box Analysis

## Overview
Comprehensive analysis of TATA promoter elements in 100 *Salmonella enterica* genomes.

## Key Findings

### Statistics
- **Total TATA matches:** 2,177,798
- **Unique sequences:** 100
- **Mean matches/sequence:** 21,777.98 (±646.27 SD)
- **Range:** 20,342 to 23,090 matches

### Distribution
- **Strand balance:** 50% forward / 50% reverse (perfect balance)
- **Normal distribution:** Confirmed by Q-Q plot
- **Conservation:** ±1.4% variation across strains

### Spatial Organization
- **Mean spacing:** 441.51 bp between TATA boxes
- **Median spacing:** 219.00 bp (consistent with bacterial gene spacing)
- **Range:** 2 bp to 8,911 bp

## Files

### Visualization Plots (9 PNG images)
- `TATA_distribution_histogram.png` - Frequency distribution + Q-Q plot
- `TATA_strand_distribution.png` - Strand bias (pie chart)
- `TATA_abundance_boxplot.png` - Abundance group comparison
- `TATA_abundance_visualization.png` - Combined box + scatter plot
- `TATA_spacing_histogram.png` - Inter-box distance distribution
- `TATA_top_10_barplot.png` - Top 10 sequences
- `TATA_vs_length.png` - Correlation with sequence length
- `TATA_vs_gc.png` - Correlation with GC content
- `TATA_cumulative_distribution.png` - Cumulative distribution function

### Analysis Data Files
- `TATA_descriptive_stats.txt` - Basic statistics
- `TATA_strand_bias.txt` - Strand distribution
- `TATA_spacing_analysis.txt` - Spacing metrics
- `TATA_abundance_classification.txt` - All sequences grouped by abundance
- `TATA_count_per_seqID.txt` - Count per sequence (sorted)
- `TATA_top_10_sequences.txt` - Top 10 sequences
- `TATA_bottom_10_sequences.txt` - Bottom 10 sequences
- `TATA_with_metadata.txt` - TATA counts with genome properties
- `TATA_Analysis_Summary.txt` - Comprehensive interpretation

### Documentation
- `PLOTS_MANIFEST.txt` - Detailed plot descriptions
- `COMPLETE_ANALYSIS_INDEX.txt` - Full deliverables index
- `README.md` - This file

## Methodology

### Data Source
- **Dataset:** seqdump_cleaned.fasta
- **Organism:** *Salmonella enterica*
- **Sequences:** 100 genomic sequences
- **Total nucleotides:** ~4.87 billion bp

### Analysis Tools
- **Pattern matching:** seqkit locate
- **Statistics:** Python 3
- **Visualization:** R 4.5.0
- **Data processing:** Unix tools

### Analysis Types Performed
✓ Descriptive statistics (mean, median, SD, quartiles)
✓ Distribution analysis (normality testing, Q-Q plots)
✓ Strand bias analysis
✓ Spatial distribution analysis (spacing patterns)
✓ Abundance grouping (High/Medium/Low classification)
✓ Correlation analysis (vs length, vs GC content)
✓ Ranked sequence identification
✓ Publication-quality visualization

## Biological Significance

### TATA Box
- Ancient promoter element found in eukaryotes and bacteria
- Consensus sequence: TATAAT (variable in real genomes)
- Located ~10 bp upstream of transcription start site
- Essential for transcription initiation

### Salmonella Findings
1. **Abundant:** ~22K TATA matches per 4.8 Mb genome
2. **Uniform:** No strand bias (perfect 50/50 distribution)
3. **Conserved:** Consistent spacing (~219 bp median) across strains
4. **Stable:** High similarity across all 100 sequences

## Implications
- Highly stable promoter architecture across *Salmonella* strains
- Consistent gene regulation machinery across pathogenic variants
- Target for evolutionary and comparative genomics studies
- Foundation for understanding virulence factor regulation

## Usage Guide

**Quick Overview (5 min):**
1. Read: `TATA_Analysis_Summary.txt`
2. View: `TATA_distribution_histogram.png`
3. View: `TATA_strand_distribution.png`

**Detailed Analysis (30 min):**
1. Read: `COMPLETE_ANALYSIS_INDEX.txt`
2. Read: `PLOTS_MANIFEST.txt`
3. View all 9 PNG plots
4. Review data files

**Statistical Details:**
- `TATA_descriptive_stats.txt`
- `TATA_spacing_analysis.txt`
- `TATA_strand_bias.txt`

## Author Notes
Analysis completed: 2026-03-23
Tools: seqkit, Python, R, Unix tools
All analysis results are reproducible from source FASTA file

## License
[Add your preferred license here]

## Contact
[Add contact information here]
