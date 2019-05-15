user="Derek"

#Pratchaya paths
if(user=="Pratchaya"){
	microbiome_custom_functions="/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R"
	full.table.sph_spring="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/PseudoSpring_GENOMES_table.txt"
}



#Derek paths
if(user=="Derek"){
	microbiome_custom_functions="/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R"
	full.table.sph_spring="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/PseudoSpring_GENOMES_table.txt"
	full.table.sph_summer="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/PseudoSummer_GENOMES_table.txt"
	location_to_write_combined_SpringSummer_table="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/PseudoSpringSummer_COMBINED_table.txt"
	all_bulk_metagenomes_metadata="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/all_bulk_metagenomes_metadata.txt"
}


source(microbiome_custom_functions)
library(gplots)

#Pseudomonas
#replace column name called 'unique' into S113 and S133 for spring and summer 
full.table.sph_spring=read.table(full.table.sph_spring, header = T, row.names = 1, sep = "\t")
full.table.sph_summer=read.table(full.table.sph_summer, header = T, row.names = 1, sep = "\t")

#create column names to row names to be ready to change 
as.character(colnames(full.table.sph_spring))->rnames.sph.spring
as.character(colnames(full.table.sph_summer))->rnames.sph.summer
gsub("_export.txt", "H103", rnames.sph.spring)->rnames.sph.spring
gsub("_export.txt", "H135", rnames.sph.summer)->rnames.sph.summer
colnames(full.table.sph_spring)=rnames.sph.spring
colnames(full.table.sph_summer)=rnames.sph.summer

as.matrix(full.table.sph_spring)->full.table.sph_spring
as.matrix(full.table.sph_summer)->full.table.sph_summer
combine_taxa_tables(full.table.sph_spring, full.table.sph_summer)->combined


#filter things less than 1000 reads
combined[,which(colSums(combined)>=20000)]->combined

#make everything mapping a 3 column field, separated by underscore
#field 1 = local or RDP or Decoy
#field 2 = genome ID
#field 3 = contig ID

gsub("_", "-", rownames(combined))->rownames(combined)  #change all underscores to hyphen
gsub("^DL133-OTU5", "DL133.OTU5", rownames(combined))->rownames(combined) 
gsub("^DL133-nonOTU5", "DL133.nonOTU5", rownames(combined))->rownames(combined) 
gsub("^nonOTU5-TK165", "TK165.nonOTU5", rownames(combined))->rownames(combined)
gsub("^OTU5-TK165", "TK165.OTU5", rownames(combined))->rownames(combined)
gsub("(^[A-Za-z0-9.]+)-(.*)", "\\1_\\2", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
gsub("(.*)-([A-Za-z0-9.]+$)", "\\1_\\2", rownames(combined))->rownames(combined) #replace last hyphen with underscore 
gsub("(.*)(PsyRun133)-(S[0-9]+)(.*)", "\\1H133\\3\\4", rownames(combined))->rownames(combined) #replace first hyphen with underscore 

write.table(combined, location_to_write_combined_SpringSummer_table, sep = '\t')
metadata.sph=read.table(all_bulk_metagenomes_metadata, sep = "\t", comment.char = "")

as.matrix(metadata.sph)->metadata.sph
pal1 <- colorRampPalette(c("black","#3f324f","#6A1B9A","#E85285", "#FFECB3"))(n=300)


#simplify genome names for color bar labels
rownames(combined)->genomelabels1
gsub("p[0-9]+.[A-Z][0-9]+", "TK.Pseudo", genomelabels1)->genomelabels1
gsub(".*H133.*", "Local_Eyach", genomelabels1)->genomelabels1
gsub(".*S[0-9]+.*", "Local_Roger", genomelabels1)->genomelabels1
gsub("At_ROOT.*", "AtROOT", genomelabels1)->genomelabels1
gsub("At_LSPHERE-Leaf.*", "AtLSPHERE", genomelabels1)->genomelabels1
gsub("At_SOIL-Soil.*", "AtSOIL", genomelabels1)->genomelabels1
gsub(".*NZ.*", "REFSEQ", genomelabels1)->genomelabels1
gsub(".*NC.*", "REFSEQ", genomelabels1)->genomelabels1
cbind(rownames(combined), genomelabels1)->genome_name_replacement_table  #check why ur rownames are the entire table



metadata.sph[match(colnames(combined), metadata.sph[,3], nomatch=0),18]->colnames(combined)



#remove rogue samples
combined[,grep("SmallClover_S59.1_Pseudobulk", colnames(combined), invert=TRUE)]->combined
combined[,grep("SpringKraut_S3.1_Pseudobulk", colnames(combined), invert=TRUE)]->combined


#arrange plants so they are together
combined[,order(colnames(combined))]->combined


combined <- normalize100(combined)

sqrt(sqrt(combined))->combined_transformed


#GENOME COLORS
unique(genome_name_replacement_table[,2])
c("REFSEQ",      "AtLSPHERE",    "AtROOT",     "AtSOIL", "Local_Eyach",    "TK.Pseudo")->  categories
c("goldenrod3", "deepskyblue3", "darkorange3", "gray52", "darkolivegreen", "black")->names(categories)
#add colors to genome name replacement table
cbind(genome_name_replacement_table, names(categories)[match(genome_name_replacement_table[,2], categories)])->genome_name_replacement_table
	#ignore heirarchical clustering and order heatmap rows by categories in the order above.
reorder_vector=order(match(genome_name_replacement_table[,2],categories))
combined_transformed=combined_transformed[reorder_vector,]
#must also reorder the genome replacement table by the same factor? no...? but doesn't hurt
genome_name_replacement_table=genome_name_replacement_table[reorder_vector,]

genome_name_replacement_table[match(rownames(combined_transformed), genome_name_replacement_table[,1]),3]->genome_colors


#SAMPLE SEASON COLORS
colnames(combined_transformed) -> rnametest
gsub(".*S[0-9]+.*", "Summer", rnametest)->rnametest
gsub(".*Og[0-9]+.*", "Spring", rnametest)->rnametest
c("Spring", "Summer")->categories
c("cadetblue3", "firebrick3")->names(categories)
cbind(colnames(combined_transformed), rnametest)->sample_name_replacement_table
#add colors to sample name replacement table
cbind(sample_name_replacement_table, names(categories)[match(sample_name_replacement_table[,2], categories)])->sample_name_replacement_table
	#ignore heirarchical clustering and order heatmap rows by categories in the order above.
reorder_vector=order(match(sample_name_replacement_table[,2],categories))
combined_transformed=combined_transformed[,reorder_vector]
#must also reorder the sample replacement table by the same factor? no...? but doesn't hurt
sample_name_replacement_table=sample_name_replacement_table[reorder_vector,]

sample_name_replacement_table[match(colnames(combined_transformed), sample_name_replacement_table[,1]),3]->season_colors



cols <- rep('black', length(colnames(combined_transformed)))

prefix <- c('Athaliana', 'Moss', 'FuzzyRosemary', 'Soil', 'SmallClover', 'Draba', 'SpringKraut', 'Cardamine', 'Plantago', 'Thistle', 'BigClover', 'Grass', 'Dandelion', '?_S4.2_Pseudobulk', 'Clover', 'NeckarKing')
color_list <- c('darkolivegreen4', 'green', 'darkviolet','darkorange4', 'green4', 'hotpink1', 'darkorchid1', 'orangered2', 'darkorchid4', 'darkslategray', 'forestgreen', 'darkgreen', 'yellow', 'black', 'darkslategray4', 'darkslateblue')
for(i in 1:length(prefix)){
  cols[grepl(prefix[i], colnames(combined_transformed))] <- color_list[i]
}


#pdf(file="/Users/pratchayapramojnaayutthaya/Desktop/Pseudomonas.pdf", width = 10, height = 10)
tmp.pdf<- heatmap.2(t(combined_transformed), density.info = "none",
                    trace = "none", dendrogram="none", colRow = cols,
                    col = pal1, main = "Pseudomonas reads per genome", 
                    margins = c(10,9), Colv=FALSE, Rowv=FALSE, cexRow = 0.41, 
                    keysize = 0.7,breaks = c(seq(0, 1, length=300),1.5),
                    key.title = NA, key.xlab = NA, 
                    ColSideColors= genome_colors, RowSideColors=season_colors,
                    labCol = F, sepwidth=c(0,0),
                    key.par=list(mar=c(2.5,0,0,30))   )
#dev.off()

duplicated(colnames(test.percent.u.sph))

colnames(test.percent.u.sph)[duplicated(colnames(test.percent.u.sph))]

grep("NK_S32_Pseudobulk", rownames(rotated.test.percent.u.sph))







