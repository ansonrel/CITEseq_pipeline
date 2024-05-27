#!/bin/bash

# use lattest salmon
alias salmon="/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon"

FOLDER="/home/Shared_sherborne/data/seq/krieg_REAPseq/aug2021"

#create ADT index. tab del file is first column IDs, second column sequences.

/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon index -t alevin_ADT.csv -i adt_index --features -k7

/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon index -t alevin_HTO.csv -i hto_index --features -k7

#I'm using the old ref I made back in march 2020

#salmon index -i index -k 31 --gencode -p 4 -t gencode.v33.transcripts.fa.gz

END=8
for i in $(seq 1 $END)
#for i in 2 3 5 7

do 
	#mkdir -p alevin/RNA/W$i
	#/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon alevin -l ISR -i /home/Shared_taupo/data/annotation/Human/gencode_GRCh38.33/cdna/index/ \
	#	-o alevin/RNA/W$i \
	#	--tgMap /home/Shared_taupo/data/annotation/Human/gencode_GRCh38.33/gtf/txp2gene.tsv \
	#	--dumpFeatures \
	#	--expectCells 20000 \
	#	-1 $FOLDER/"$i"D_S1_L001_R1_001.fastq.gz \
	#	-2 $FOLDER/"$i"D_S1_L001_R2_001.fastq.gz \
	#	--chromiumV3 \
	#	-p 20


#quantify ADTs
# may want to remove barcodeLength, umiLength and end, as take as extra protocol
	mkdir -p alevin/ADT/W$i
	/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon alevin -l ISR -i adt_index \
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

/home/asonrel/soft/salmon-1.5.2_linux_x86_64/bin/salmon alevin -l ISR -i hto_index \
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
