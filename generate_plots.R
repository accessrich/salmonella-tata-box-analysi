#!/usr/bin/env Rscript

cat("Generating TATA visualization plots with R...\n")

# 1. DISTRIBUTION HISTOGRAM
cat("1. Creating TATA distribution histogram...\n")
counts <- read.table('TATA_count_per_seqID.txt', skip=1)
tata_counts <- counts[,1]

png('TATA_distribution_histogram.png', width=1200, height=600, res=100)
par(mfrow=c(1,2))

# Histogram with mean and median
hist(tata_counts, breaks=30, col='steelblue', border='black', 
     xlab='TATA Matches per Sequence', ylab='Number of Sequences',
     main='Distribution of TATA Box Counts', cex.main=1.5, cex.lab=1.2)
abline(v=mean(tata_counts), col='red', lwd=2, lty=2)
abline(v=median(tata_counts), col='green', lwd=2, lty=2)
legend('topright', legend=c(paste0('Mean: ', round(mean(tata_counts))),
                            paste0('Median: ', median(tata_counts))),
       col=c('red', 'green'), lty=c(2,2), lwd=2)

# Q-Q plot for normality check
qqnorm(tata_counts, main='Q-Q Plot (Normality Check)', cex.main=1.5)
qqline(tata_counts, col='red', lwd=2)

dev.off()
cat("   ✓ TATA_distribution_histogram.png\n")

# 2. STRAND BIAS PIE CHART
cat("2. Creating strand bias pie chart...\n")
png('TATA_strand_distribution.png', width=800, height=600, res=100)
par(mar=c(1,1,3,1))
strand_counts <- c(1088899, 1088899)
pie(strand_counts, labels=c('+', '-'), 
    main='TATA Distribution by Strand', 
    col=c('#66c2a5', '#fc8d62'),
    cex.main=1.5, cex=1.3,
    percentage=TRUE)
dev.off()
cat("   ✓ TATA_strand_distribution.png\n")

# 3. ABUNDANCE GROUP BOX PLOT
cat("3. Creating abundance group box plot...\n")
abundance <- read.table('TATA_abundance_classification.txt', header=TRUE, sep='\t')

png('TATA_abundance_boxplot.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))
boxplot(abundance$TATA_Count ~ factor(abundance$Abundance_Group, levels=c('High', 'Medium', 'Low')),
        col=c('#8dd3c7', '#ffffb3', '#bebada'),
        xlab='Abundance Group', ylab='TATA Count',
        main='TATA Count Distribution by Abundance Group',
        cex.lab=1.2, cex.main=1.5, cex=1.1)
grid(NA, NULL)
dev.off()
cat("   ✓ TATA_abundance_boxplot.png\n")

# 4. SPACING DISTRIBUTION
cat("4. Creating spacing distribution histogram...\n")
spacings <- as.numeric(readLines('spacings.txt'))

png('TATA_spacing_histogram.png', width=1400, height=600, res=100)
par(mfrow=c(1,2), mar=c(5, 5, 4, 2))

# Full range
hist(spacings, breaks=100, col='coral', border='black',
     xlab='Distance between TATA boxes (bp)', ylab='Frequency',
     main='TATA Box Spacing Distribution (Full Range)',
     cex.lab=1.1, cex.main=1.3)
grid(NA, NULL)

# Zoomed (< 500 bp)
spacings_zoom <- spacings[spacings < 500]
hist(spacings_zoom, breaks=50, col='skyblue', border='black',
     xlab='Distance between TATA boxes (bp)', ylab='Frequency',
     main='TATA Box Spacing Distribution (< 500 bp)',
     cex.lab=1.1, cex.main=1.3)
grid(NA, NULL)

dev.off()
cat("   ✓ TATA_spacing_histogram.png\n")

# 5. TOP 10 SEQUENCES BAR PLOT
cat("5. Creating top 10 sequences bar plot...\n")
top10 <- read.table('TATA_top_10_sequences.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

png('TATA_top_10_barplot.png', width=1400, height=600, res=100)
par(mar=c(8, 5, 4, 2))
barplot(top10$Count, names.arg=top10$seqID,
        col='#1f77b4', border='black',
        xlab='Sequence ID', ylab='TATA Count',
        main='Top 10 Sequences by TATA Count',
        cex.lab=1.2, cex.main=1.5, cex.names=0.9, las=2)
grid(NA, NULL)
dev.off()
cat("   ✓ TATA_top_10_barplot.png\n")

# 6. SCATTER PLOT - TATA vs Sequence Length
cat("6. Creating TATA vs sequence length scatter plot...\n")
metadata <- read.table('TATA_with_metadata.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

png('TATA_vs_length.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))
length_mb <- as.numeric(metadata$Sequence_Length) / 1e6
plot(length_mb, metadata$TATA_Count,
     main='TATA Count vs Sequence Length',
     xlab='Sequence Length (Mb)', ylab='TATA Count',
     pch=19, col='#2ca02c', cex=1.2, cex.lab=1.2, cex.main=1.5)
grid()

# Calculate correlation
corr_len <- cor(length_mb, metadata$TATA_Count, use='complete.obs')
p_val_len <- cor.test(length_mb, metadata$TATA_Count)$p.value
legend('topright', legend=c(paste0('r = ', round(corr_len, 3)),
                            paste0('p-value = ', format(p_val_len, digits=2, scientific=TRUE))),
       cex=1.1, box.lwd=1.5)

dev.off()
cat("   ✓ TATA_vs_length.png\n")

# 7. SCATTER PLOT - TATA vs GC Content
cat("7. Creating TATA vs GC content scatter plot...\n")

png('TATA_vs_gc.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))
gc_pct <- as.numeric(metadata$GC_Percent)
plot(gc_pct, metadata$TATA_Count,
     main='TATA Count vs GC Content',
     xlab='GC Content (%)', ylab='TATA Count',
     pch=19, col='#9467bd', cex=1.2, cex.lab=1.2, cex.main=1.5)
grid()

# Calculate correlation
corr_gc <- cor(gc_pct, metadata$TATA_Count, use='complete.obs')
p_val_gc <- cor.test(gc_pct, metadata$TATA_Count)$p.value
legend('topright', legend=c(paste0('r = ', round(corr_gc, 3)),
                            paste0('p-value = ', format(p_val_gc, digits=2, scientific=TRUE))),
       cex=1.1, box.lwd=1.5)

dev.off()
cat("   ✓ TATA_vs_gc.png\n")

# 8. CUMULATIVE DISTRIBUTION
cat("8. Creating cumulative distribution plot...\n")
sorted_counts <- sort(tata_counts)
cumulative <- (1:length(sorted_counts)) / length(sorted_counts) * 100

png('TATA_cumulative_distribution.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))
plot(sorted_counts, cumulative, type='l', lwd=2, col='darkblue',
     main='Cumulative Distribution of TATA Counts',
     xlab='TATA Count per Sequence', ylab='Cumulative Percentage (%)',
     cex.lab=1.2, cex.main=1.5)
polygon(c(sorted_counts, rev(sorted_counts)), 
        c(cumulative, rep(0, length(cumulative))),
        col=rgb(0.3, 0.3, 0.8, 0.3), border=NA)
grid()
dev.off()
cat("   ✓ TATA_cumulative_distribution.png\n")

# 9. COMBINED ABUNDANCE VISUALIZATION
cat("9. Creating combined abundance visualization...\n")
png('TATA_abundance_visualization.png', width=1200, height=600, res=100)
par(mfrow=c(1,2), mar=c(5,5,4,2))

# Box plot
boxplot(abundance$TATA_Count ~ factor(abundance$Abundance_Group, levels=c('High', 'Medium', 'Low')),
        col=c('#8dd3c7', '#ffffb3', '#bebada'),
        xlab='Abundance Group', ylab='TATA Count',
        main='Box Plot by Abundance Group',
        cex.lab=1.1, cex.main=1.3)
grid(NA, NULL)

# Violin-like plot with points
abundance_factor <- factor(abundance$Abundance_Group, levels=c('High', 'Medium', 'Low'))
plot(abundance_factor, abundance$TATA_Count,
     col=c('#8dd3c7', '#ffffb3', '#bebada')[as.numeric(abundance_factor)],
     xlab='Abundance Group', ylab='TATA Count',
     main='Distribution with Individual Points',
     cex.lab=1.1, cex.main=1.3, cex=1.2)
grid(NA, NULL)

dev.off()
cat("   ✓ TATA_abundance_visualization.png\n")

cat("\n")
cat("==============================================================================\n")
cat("All 9 plots generated successfully!\n")
cat("==============================================================================\n")

