cat("Generating Resistance & Virulence Gene Visualizations...\n")

# Create summary data frame
resistance_status <- data.frame(
  Category = c("Beta-lactamases", "Tetracycline", "Aminoglycosides", "DHFR", "Chloramphenicol"),
  Genes_Searched = c(5, 3, 2, 2, 1),
  Genes_Found = c(0, 0, 0, 0, 0)
)

# Create bar plot - Resistance Gene Search Results
png('Salmonella_resistance_search_results.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))
x <- barplot(resistance_status$Genes_Found, 
             names.arg=resistance_status$Category,
             col='steelblue', border='black',
             xlab='Resistance Gene Category', ylab='Number of Genes Found',
             main='Antibiotic Resistance Genes - Search Results\n(Salmonella Reference Strains)',
             ylim=c(0, max(resistance_status$Genes_Searched) + 1),
             cex.lab=1.2, cex.main=1.4)

# Add "Searched" as overlay
points(x, resistance_status$Genes_Searched, pch=21, bg='red', cex=2)
legend('topright', legend=c('Found', 'Searched'),
       pch=c(22, 21), pt.bg=c('steelblue', 'red'), cex=1.2)

dev.off()
cat("✓ Salmonella_resistance_search_results.png\n")

# Create interpretation plot
png('Salmonella_phenotype_interpretation.png', width=1200, height=600, res=100)
par(mfrow=c(1, 2), mar=c(3, 4, 4, 2))

# Left: Antibiotic Susceptibility Profile
antibiotics <- c('Beta-lactams', 'Fluoroquinolones', 
                 'Macrolides', 'Aminoglycosides',
                 'Tetracyclines', 'Chloramphenicol',
                 'Trimethoprim')
susceptibility <- c(1, 1, 1, 1, 1, 1, 1)

barplot(susceptibility, names.arg=antibiotics,
        col='#2ecc71', border='black', ylim=c(0, 1.2),
        main='Expected Antibiotic Susceptibility\n(Wild-type Salmonella)',
        ylab='Susceptibility', cex.main=1.3, cex.lab=1.1)
text(3.5, 1.05, 'No resistance genes detected\n→ Standard regimens work',
     cex=1.0, adj=0.5)

# Right: Virulence potential
virulence_features <- c('T3SS', 'T6SS', 'SPIs', 'Flagella', 'Phase\nVariation')
presence <- c(1, 1, 1, 1, 1)

barplot(presence, names.arg=virulence_features,
        col='#e74c3c', border='black', ylim=c(0, 1.2),
        main='Expected Virulence Features\n(Salmonella enterica)',
        ylab='Likely Present', cex.main=1.3, cex.lab=1.1)
text(3, 0.6, 'Key pathogenic systems\nexpected in all strains',
     cex=1.0, adj=0.5)

dev.off()
cat("✓ Salmonella_phenotype_interpretation.png\n")

# Create comparison visualization
png('Salmonella_resistance_phenotype.png', width=1000, height=600, res=100)
par(mar=c(5, 5, 4, 2))

phenotypes <- c('Antibiotic\nResistant\n(MDR)', 'Antibiotic\nSensitive\n(Wild-type)',
                'Virulent\n(Pathogenic)','Commensal')
observed <- c(0, 1, 1, 0)
colors <- c('#e74c3c', '#2ecc71', '#e74c3c', '#95a5a6')

barplot(observed, names.arg=phenotypes,
        col=colors, border='black', ylim=c(0, 1.3),
        main='Salmonella Phenotype Profile\n(Reference Genomes)',
        ylab='Status', cex.lab=1.2, cex.main=1.4,
        cex.names=1.0)

dev.off()
cat("✓ Salmonella_resistance_phenotype.png\n")

cat("\nVisualization complete!\n")
