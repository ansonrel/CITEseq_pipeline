#!/bin/bash

# use lattest salmon
alias salmon="/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon"

# Path to seq folder with fasts
FOLDER="/home/Shared_sherborne/data/seq/"

# Number of samples
SAMPLES=8

#-------------------------------------------------------------------------------
#create ADT index. tab del file is first column IDs, second column sequences.
# alevin_ADT.csv and alevin_HTO.csv : column 1 = ADT/HTO name, column 2 = barcode

salmon index -t alevin_ADT.csv -i adt_index --features -k7

salmon index -t alevin_HTO.csv -i hto_index --features -k7

# generating an ref, if needed
#salmon index -i index -k 31 --gencode -p 4 -t gencode.v33.transcripts.fa.gz

for i in $(seq 1 $SAMPLES)

do 
	mkdir -p alevin/RNA/W$i
	salmon alevin -l ISR -i /home/Shared_taupo/data/annotation/Human/gencode_GRCh38.33/cdna/index/ \
		-o alevin/RNA/W$i \
		--tgMap /home/Shared_taupo/data/annotation/Human/gencode_GRCh38.33/gtf/txp2gene.tsv \
		--dumpFeatures \
		--expectCells 20000 \
		-1 $FOLDER/"$i"D_S1_L001_R1_001.fastq.gz \
		-2 $FOLDER/"$i"D_S1_L001_R2_001.fastq.gz \
		--chromiumV3 \
		-p 20


#quantify ADTs
# may want to remove barcodeLength, umiLength and end, as take as extra protocol
# other parameters are pretty default
	mkdir -p alevin/ADT/W$i
	salmon alevin -l ISR -i adt_index \
		-1 $FOLDER/"$i"-ADT_R1_001.fastq.gz \
		-2 $FOLDER/"$i"-ADT_R2_001.fastq.gz \
		-o alevin/ADT/W$i \
		-p 16 \
		--citeseq \
		--featureStart 0 \
		--featureLength 15 \
		--expectCells 20000 \
		--naiveEqclass

		#--end 5 \
		#		--umiLength 10 \
		#--barcodeLength 16 \


	mkdir -p alevin/HTO/W$i

	salmon alevin -l ISR -i hto_index \
		-1 $FOLDER/"$i"-HTO_R1_001.fastq.gz \
		-2 $FOLDER/"$i"-HTO_R2_001.fastq.gz \
		-o alevin/HTO/W$i \
		-p 16 \
		--citeseq \
		--featureStart 0 \
		--featureLength 15 \
		--expectCells 20000 \
		--naiveEqclass

		#--end 5 \
		#--umiLength 10 \
		#--barcodeLength 16 \

done
