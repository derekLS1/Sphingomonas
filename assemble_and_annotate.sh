#!/bin/sh
#reserve running with 16 CPUs for this job
#$ -pe parallel 16
#request 25GB of RAM
#$ -l h_vmem=25G
#use /bin/bash to execute this script
#$ -S /bin/bash
#run this job from current working directory
#$ -cwd

#CONFIRM PROGRAM LOCATIONS
#for running the SPAdes programme, identify the directory of the SPAdes programme
spades=/ebio/abt6_projects9/microbiome_analysis/data/software/SPAdes-3.11.0/bin/spades.py
pilon=/ebio/abt6_projects7/bacterial_strain_analysis/code/pilon-1.20.jar
rename_headers=/ebio/abt6_projects7/bacterial_strain_analysis/code/rename_fasta.py
prokka=/ebio/abt6_projects7/bacterial_strain_analysis/code/prokka/bin/prokka
quast=/ebio/abt6_projects7/bacterial_strain_analysis/code/quast-5.0.2/quast-lg.py
BUSCO=/ebio/abt6_projects9/microbiome_analysis/data/software/busco/BUSCO.py
BUSCO_database=/ebio/abt6_projects9/microbiome_analysis/data/software/busco/datasets/proteobacteria_odb9

#set a date variable
day=$(date +'%Y%m%d')

#These are the variables passed in from the qsub command.
raw_data_directory=$1     #path to a folder containing symbolic links to all raw data files (read 1 and read 2 hiseq)
read1=$2                  #the read1 filename (the read 2 filename will be constructed from this)
assembly_directory=$3     #path to a folder in which all the genome assemblies will be stored
batchname=$4            #prefix to give all genome assemblies (for example, the number of the HiSeq run, or any helpful identifier for the batch)

#parse the read1 filename to get useful elements for making new names
lane=$(sed 's/.*S.*_L00\(.*\)_R.*/\1/' <(echo "$read1"))
echo "$lane"
strain=$(sed 's/.*S\(.*\)_L00.*/\1/' <(echo "$read1"))   #every sample in a lane has an S (sample) number. 
echo "$strain"
echo "$read1"
read2=$(sed 's/R1/R2/g' <(echo "$read1"))   #filename for read 2 is simply made by changing R1 -> R2
echo "$read2"

#now create short and unique genome names for processing genome assemblies
genome_name=$(echo S"$strain""$batchname"_L"$lane")

#go to $assembly_directory
cd $assembly_directory

#now run SPAdes
nice -10 $spades -1 "$raw_data_directory""/""$read1" -2 "$raw_data_directory""/""$read2" -o "$assembly_directory""/""$genome_name" --careful -k 21,33,55,77,99,127 --threads 16

#now to get contig files containing all contigs, go to the genome_name folders
cd "$assembly_directory""/""$genome_name"

#remove unncessary files, and keep only contig fasta files
for file in $(ls | grep -v contigs.fasta)   #grep -v keeps everything but what you search for
do
  rm $file
  rm -r $file
done

#now rename the contig files
mv contigs.fasta "$genome_name"_contigs.fasta

#now we do indexing using BWA programmes and the fasta files of contigs
bwa index "$genome_name"_contigs.fasta

#now do mapping of original raw reads back onto the contigs fasta and create bam files
bwa mem -t 10 "$genome_name"_contigs.fasta "$raw_data_directory""/""$read1" "$raw_data_directory""/""$read2" | samtools view -q 30 -b -@ 10 - | samtools sort -@ 10 - > "$genome_name"_contigs.bam

#then run Pilon programme to detect variants of reads and correct indels.  
samtools index "$genome_name"_contigs.bam     
java -Xmx16G -jar $pilon --genome "$genome_name"_contigs.fasta --frags "$genome_name"_contigs.bam

#Now rename contig headnames so that they are short and compatible with Pan Genome Analysis
python $rename_headers pilon.fasta

#remove anything that says unknown description
sed -i -e 's/<unknown description>//g' pilon_renamed.fasta

#finally rename the pilon fasta file name into the renamed pilon fasta files
mv pilon_renamed.fasta "$genome_name"_pilon.fasta

#renome all un-renamed pilon fasta files
rm pilon.fasta

#now annotate the pilon-corrected fasta files using Prokka programme
mkdir Annotations
"$prokka" "$genome_name"_pilon.fasta --force --compliant --outdir Annotations
cp "$assembly_directory"/"$genome_name"/Annotations/PROKKA_*.gbk "$assembly_directory"/"$genome_name"/Annotations/"$genome_name"_"$day".gbk

#Calculate coverage
samtools depth -a "$genome_name"_contigs.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > "$genome_name"_coverage.txt

#Calculate the N50 stats
sh /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/genomeassemblies/N50.sh "$genome_name"_pilon.fasta > "$genome_name"_N50.txt

#combine coverage and N50 stats and gene count into one stat .txt file
N50=$(grep "N50" "$genome_name"_N50.txt)
coverage=$(grep "Average" "$genome_name"_coverage.txt)
genes=$(grep "gene" "$assembly_directory"/"$genome_name"/Annotations/PROKKA*.txt)
echo -e "$g" '\t'  "$N50"  '\t'   "$coverage"  '\t'   "$genes" >> "$genome_name"_qualitystats.txt

#remove the N50 and coverage separate .txt files
rm "$genome_name"_coverage.txt
rm "$genome_name"_N50.txt

#run QUAST
$quast *_pilon.fasta -o QUASToutput

#run BUSCO program to check genome completeness
out=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/BUSCOoutput 
cd $out

$BUSCO -i "$assembly_directory"/"$genome_name"/Annotations/PROKKA*.faa -o "$genome_name".busco -l $BUSCO_database -m prot -f



################################################################################

#RUN THIS ON COMMAND LINE
raw_data_directory=/ebio/abt6_projects7/bacterial_strain_analysis/data/raw_seqDATA/HiSeq3000-Run0133/Sphingomonas
assembly_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
#name can be changed
name_prefix="SphRun133_"
#read1=nextera_S*_L003_R1_001.fastq.gz
#echo "$read1"
for g in $(ls $raw_data_directory | grep _R1_)
do echo "$g"
read1=$g
lane=$(sed 's/.*S.*_L00\(.*\)_R.*/\1/' <(echo "$read1"))
echo "$lane"
strain=$(sed 's/.*S\(.*\)_L00.*/\1/' <(echo "$read1"))
echo "$strain"
genome_name=$(echo "$name_prefix"S"$strain"_L"$lane")
echo "$genome_name"
mkdir -p $assembly_directory/$genome_name

qsub -N Job"$strain" -o "$assembly_directory"/"$genome_name"/out.txt -e "$assembly_directory"/"$genome_name"/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/genome_assembly.sh $raw_data_directory $read1 $assembly_directory $name_prefix
done


#read1=S10_L008_R1_001.fastq.gz
#RUN THIS COMMAND FOR Pseudomonas
#RUN THIS ON COMMAND LINE
raw_data_directory=/ebio/abt6_projects7/bacterial_strain_analysis/data/raw_seqDATA/HiSeq3000-Run0133/Pseudomonas
assembly_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Pseudomonas
#name can be changed
name_prefix="PsyRun133_"
#read1=nextera_S*_L003_R1_001.fastq.gz
#echo "$read1"
for g in $(ls $raw_data_directory | grep _R1_)
do echo "$g"
read1=$g
lane=$(sed 's/.*S.*_L00\(.*\)_R.*/\1/' <(echo "$read1"))
echo "$lane"
strain=$(sed 's/.*S\(.*\)_L00.*/\1/' <(echo "$read1"))
echo "$strain"
genome_name=$(echo "$name_prefix"S"$strain"_L"$lane")
echo "$genome_name"
mkdir -p $assembly_directory/$genome_name

qsub -N Job"$strain" -o "$assembly_directory"/"$genome_name"/out.txt -e "$assembly_directory"/"$genome_name"/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/genome_assembly.sh $raw_data_directory $read1 $assembly_directory $name_prefix
done
