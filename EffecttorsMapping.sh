#!/bin/bash

# Call the programmes used in the following bash pipeline
module load AdapterRemoval/v2.3.1 
module load perl/v5.28.1 
module load bwa/v0.7.15
module load htslib/v1.6
module load samtools/v1.6 
module load bedtools/v2.29.0

# Definition of variables
PROJECT=/groups/hologenomics/sarai/data/temporal
FASTQ=$PROJECT/fastqs
ADAPTER=$PROJECT/adapRem
MAP=$PROJECT/mapping
DEPTH=$PROJECT/depth

# References genomes 
MORYZEFFECTORS=/groups/hologenomics/sarai/data/temporal/references/moryzEffectors_Control2_Known4_Thorsten3_prediction165__extraprediction3_genomicSequence_NOTduplicated_shortNames.fasta


#index in case to be neccessary
bwa index $MORYZEFFECTORS


##### ADAPTERS #####

## Remove the illumina adapters P5 and P7 using AdapterRemoval 
#Input: fastq files

echo "Remove adapters"
# -p: no error if existing 
mkdir -p $ADAPTER && cd $ADAPTER

if [ ! -e .adap.done ]; then
	for f in $FASTQ/*.fastq.gz 
		# Extract the names of the samples
		do
		bn=$(basename $f .fastq.gz)
    		
 		## Run Adapter removal 
		# collapse: Combined into a single consensus sequence pair aligments
		# Output: output_paired.collapsed containing merged reads,
		# output_paired.collapsed.truncated containing merged reads that have been trimmed due to the --trimns or --trimqualities options.
		# The sequencing primers are not specified, since AdapterRemoval knows the Illumina sequencing primers.
		echo "AdapterRemoval  --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
			--collapse --basename ${bn}_noAdap --file1 $f --file2 ${f/R1/R2}"
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J adap --
 	touch .adap.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy


# To avoid remove relevant reads, we will merge collapsed and uncollapsed reads
cd $ADAPTER
for f in $(ls *_noAdap.collapsed.gz)
do bn=$(basename $f _noAdap.collapsed.gz); cat ${bn}*.gz > ${bn}_merged_NoAdap.fastq.gz
done

#Wait until the jobs above finish 
echo "Waiting"
read dummy



##### MAPPING #####
# Map to the effectors fasta file using BWA mem
# Use mem so that we can soft clip the primer sequences in the beginning of the read
# Sort the mapping by coordinates and mark duplicates using samtools 

echo "Map to effectors data base"
# -p: no error if existing 
mkdir -p $MAP && cd $MAP

if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/*_merged_NoAdap.fastq.gz;
  		do
    		bn=$(basename $f *_merged_NoAdap.fastq.gz )
    		# Run bwa mem and then sort then sort the bam by coordinates, to be able to mark duplicates 
			# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $MORYZEFFECTORS $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MapTOmoryzEffectors.markdup.bam   )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy



##### DEPTH OF COVERAGE #####
# Get the depth of coverage of the efectors in every samples

echo "Calculate the depth of coverage"
# -p: no error if existing 
mkdir -p $DEPTH && cd $DEPTH

if [ ! -e .depth.done ]; then
	## Run samtools
  	for f in $MAP/*_MapTOmoryzEffectors.markdup.bam;
  		do
    		bn=$(basename $f *_MapTOmoryzEffectors.markdup.bam  )
    		echo "(samtools index $f; samtools depth -aa $f > ${bn}_MapTOmoryzEffectors_depth.txt  )"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J map -R --max-array-jobs=10 --
  touch .depth.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy


