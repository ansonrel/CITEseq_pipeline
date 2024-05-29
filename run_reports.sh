#!/bin/sh -x

# Number of samples
BEGIN=1
END=8

# R libraries
RLIB=R/x86_64-pc-linux-gnu-library/4.0/
RPATH=/usr/local/R/R-4.0.0/bin/R

# output of step 1
INPUTFOLDER=alevin/HTO

# file with patient/ sample info and (optional) HTO debarcoding matrix. 
# example file in metadata/
METADATA=metadata/metadata.csv

# ---------------------------------------
# Compile the reports

mkdir -p sce/

for i in $(seq $BEGIN $END)

do 
  
	mkdir -p reports//W$i
		
	R_LIBS=$RLIB \
	$RPATH CMD BATCH \
	--no-restore --no-save \
	"--args inputFolder='$INPUTFOLDER/W$i/alevin/quants_mat.gz' fileName='W$i' metadataFile='$METADATA' rmdtemplate='1_HTO_diagnostics.Rmd' outputdir='reports/W$i' outputfile='HTO_diagnostics_W$i.html'" scripts/run_render.R logs/HTO_diagnostics_W$i.log

	R_LIBS=$RLIB $RPATH CMD BATCH --no-restore --no-save "--args inputFolder='$INPUTFOLDER/W$i/alevin/quants_mat.gz' fileName='W$i' rmdtemplate='2_HTO_debarcoding.Rmd' outputdir='reports/W$i' outputfile='HTO_debarcoding_W$i.html'" scripts/run_render.R logs/HTO_debarcoding_W$i.log
  
done

rmarkdown::render("3_RNA_diagnostics.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("4_Ab_diagnostics.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("5_Ab_seuratScaling.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("6_use_celltypes.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("7_markers_celltypes.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("8_Monocle_pseudotime.Rmd", output_format="html_document", output_dir = 'reports')
rmarkdown::render("9_DS_DA_analysis.Rmd", output_format="html_document", output_dir = 'reports')
