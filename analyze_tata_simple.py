#!/usr/bin/env python3
import re
from collections import defaultdict

print("=" * 80)
print("COMPREHENSIVE TATA ANALYSIS")
print("=" * 80)

# 1. DESCRIPTIVE STATISTICS
print("\n1. DESCRIPTIVE STATISTICS")
print("=" * 80)

tata_counts = defaultdict(int)
strand_counts = defaultdict(int)
spacing_data = []
prev_positions = defaultdict(int)

with open('TATA_locate_results.txt') as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        parts = line.strip().split('\t')
        seqid = parts[0]
        strand = parts[3]
        start = int(parts[4])
        
        tata_counts[seqid] += 1
        strand_counts[strand] += 1
        
        # Calculate spacing
        if seqid in prev_positions:
            spacing = start - prev_positions[seqid]
            if spacing > 0:
                spacing_data.append(spacing)
        prev_positions[seqid] = start

counts = list(tata_counts.values())
total_matches = sum(counts)

print(f"Total Matches:        {total_matches}")
print(f"Unique Sequences:     {len(counts)}")
print(f"Mean Matches/Seq:     {sum(counts)/len(counts):.2f}")
print(f"Median Matches/Seq:   {sorted(counts)[len(counts)//2]}")
print(f"Std Dev:              {(sum((x - sum(counts)/len(counts))**2 for x in counts) / len(counts))**0.5:.2f}")
print(f"Min:                  {min(counts)}")
print(f"Max:                  {max(counts)}")

stats_output = [
    "Statistic\tValue",
    f"Total_Matches\t{total_matches}",
    f"Unique_Sequences\t{len(counts)}",
    f"Mean_Matches_Per_Sequence\t{sum(counts)/len(counts):.2f}",
    f"Median_Matches_Per_Sequence\t{sorted(counts)[len(counts)//2]}",
    f"Min\t{min(counts)}",
    f"Max\t{max(counts)}"
]
with open('TATA_descriptive_stats.txt', 'w') as f:
    f.write('\n'.join(stats_output))

# 2. STRAND BIAS
print("\n2. STRAND BIAS ANALYSIS")
print("=" * 80)
print("Strand\tCount\t\tPercentage")
with open('TATA_strand_bias.txt', 'w') as f:
    f.write("Strand\tCount\tPercentage\n")
    for strand in ['+', '-']:
        count = strand_counts.get(strand, 0)
        pct = (count / total_matches * 100) if total_matches > 0 else 0
        print(f"{strand}\t{count}\t\t{pct:.2f}%")
        f.write(f"{strand}\t{count}\t{pct:.2f}%\n")

# 3. SPATIAL DISTRIBUTION
print("\n3. SPATIAL DISTRIBUTION ANALYSIS")
print("=" * 80)
if spacing_data:
    spacing_data.sort()
    print(f"Total spacing measurements: {len(spacing_data)}")
    print(f"Mean spacing:               {sum(spacing_data)/len(spacing_data):.2f} bp")
    print(f"Median spacing:             {spacing_data[len(spacing_data)//2]:.2f} bp")
    print(f"Min spacing:                {min(spacing_data)} bp")
    print(f"Max spacing:                {max(spacing_data)} bp")
    
    spacing_output = [
        "Metric\tValue",
        f"Total_Measurements\t{len(spacing_data)}",
        f"Mean_Spacing_bp\t{sum(spacing_data)/len(spacing_data):.2f}",
        f"Median_Spacing_bp\t{spacing_data[len(spacing_data)//2]}",
        f"Min_Spacing_bp\t{min(spacing_data)}",
        f"Max_Spacing_bp\t{max(spacing_data)}"
    ]
    with open('TATA_spacing_analysis.txt', 'w') as f:
        f.write('\n'.join(spacing_output))

# 4. ABUNDANCE GROUPING
print("\n4. ABUNDANCE GROUPING")
print("=" * 80)
sorted_counts = sorted(tata_counts.items(), key=lambda x: x[1], reverse=True)
third = len(sorted_counts) // 3

high = sorted_counts[:third]
medium = sorted_counts[third:2*third]
low = sorted_counts[2*third:]

print(f"High abundance:   {len(high)} sequences (mean: {sum(c for _, c in high)/len(high):.0f} matches)")
print(f"Medium abundance: {len(medium)} sequences (mean: {sum(c for _, c in medium)/len(medium) if medium else 0:.0f} matches)")
print(f"Low abundance:    {len(low)} sequences (mean: {sum(c for _, c in low)/len(low) if low else 0:.0f} matches)")

abundance_output = ["seqID\tTATA_Count\tAbundance_Group"]
for seqid, count in high:
    abundance_output.append(f"{seqid}\t{count}\tHigh")
for seqid, count in medium:
    abundance_output.append(f"{seqid}\t{count}\tMedium")
for seqid, count in low:
    abundance_output.append(f"{seqid}\t{count}\tLow")
with open('TATA_abundance_classification.txt', 'w') as f:
    f.write('\n'.join(abundance_output))

# 5. TOP 10 SEQUENCES
print("\n5. TOP 10 SEQUENCES WITH MOST TATA")
print("=" * 80)
print("Rank\tseqID\t\t\tCount")
top_output = ["Rank\tseqID\tCount"]
for i, (seqid, count) in enumerate(sorted_counts[:10], 1):
    print(f"{i}\t{seqid:<20}\t{count}")
    top_output.append(f"{i}\t{seqid}\t{count}")
with open('TATA_top_10_sequences.txt', 'w') as f:
    f.write('\n'.join(top_output))

# 6. BOTTOM 10 SEQUENCES
print("\n6. BOTTOM 10 SEQUENCES WITH LEAST TATA")
print("=" * 80)
print("Rank\tseqID\t\t\tCount")
bottom_output = ["Rank\tseqID\tCount"]
for i, (seqid, count) in enumerate(sorted_counts[-10:], len(sorted_counts)-9):
    print(f"{i}\t{seqid:<20}\t{count}")
    bottom_output.append(f"{i}\t{seqid}\t{count}")
with open('TATA_bottom_10_sequences.txt', 'w') as f:
    f.write('\n'.join(bottom_output))

# 7. SEQUENCE PROPERTIES
print("\n7. SEQUENCE PROPERTIES CORRELATION")
print("=" * 80)

props_output = ["seqID\tTATA_Count\tSequence_Length\tGC_Percent"]
import csv
with open('seqdump_metadata.csv') as f:
    reader = csv.DictReader(f)
    metadata = {row['accession']: row for row in reader}

count = 0
for seqid, tata_count in tata_counts.items():
    if seqid in metadata:
        m = metadata[seqid]
        props_output.append(f"{seqid}\t{tata_count}\t{m['length']}\t{m['gc_pct']}")
        count += 1

with open('TATA_with_metadata.txt', 'w') as f:
    f.write('\n'.join(props_output))

# Calculate TATA density
if count > 0:
    lengths = [int(metadata[sid]['length']) for sid in tata_counts if sid in metadata]
    counts_with_meta = [tata_counts[sid] for sid in tata_counts if sid in metadata]
    densities = [c / l * 1000 for c, l in zip(counts_with_meta, lengths)]
    print(f"Sequences with metadata: {count}")
    print(f"TATA Density (matches per Kb):")
    print(f"  Mean:   {sum(densities)/len(densities):.2f}")
    print(f"  Min:    {min(densities):.2f}")
    print(f"  Max:    {max(densities):.2f}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print("\nOutput files generated:")
print("  • TATA_descriptive_stats.txt")
print("  • TATA_strand_bias.txt")
print("  • TATA_spacing_analysis.txt")
print("  • TATA_abundance_classification.txt")
print("  • TATA_top_10_sequences.txt")
print("  • TATA_bottom_10_sequences.txt")
print("  • TATA_with_metadata.txt")
print("=" * 80)

