#!/bin/sh -x

BEGIN=1
END=3
cell=TB

mkdir -p sce/$cell

for i in $(seq $BEGIN $END)

do 
  
	mkdir -p reports/$cell/W$i

	R_LIBS=/home/sorjuela/R/x86_64-pc-linux-gnu-library/4.0/ \
	/usr/local/R/R-4.0.0/bin/R CMD BATCH \
	--no-restore --no-save \
	"--args cellType='$cell' inputFolder='/home/Shared_s3it/sorjuela/CITEseq/round8/alevin/HTO/W$i/alevin/quants_mat.gz' fileName='W$i' rmdtemplate='1_HTO_diagnostics.Rmd' outputdir='reports/$cell/W$i' outputfile='HTO_diagnostics_W$i.html'" scripts/run_render.R logs/HTO_diagnostics_W$i.log

	R_LIBS=/home/sorjuela/R/x86_64-pc-linux-gnu-library/4.0/ /usr/local/R/R-4.0.0/bin/R CMD BATCH --no-restore --no-save "--args cellType='$cell' inputFolder='/home/Shared_s3it/sorjuela/CITEseq/round8/alevin/HTO/W$i/alevin/quants_mat.gz' fileName='W$i' rmdtemplate='2_HTO_debarcoding.Rmd' outputdir='reports/$cell/W$i' outputfile='HTO_debarcoding_W$i.html'" scripts/run_render.R logs/HTO_debarcoding_W$i.log
  
done

#rmarkdown::render("1_HTO_diagnostics_fullsamples.Rmd", output_format="html_document", output_dir = 'reports/M')
rmarkdown::render("3_RNA_diagnostics_TB.Rmd", output_format="html_document", output_dir = 'reports/TB')
rmarkdown::render("4_Ab_diagnostics.Rmd", output_format="html_document", output_dir = 'reports/TB')
rmarkdown::render("5_use_ADT_celltypes.Rmd", output_format="html_document", output_dir = 'reports/TB')
rmarkdown::render("6_DS_analysis_TB.Rmd", output_format="html_document", output_dir = 'reports/TB')
