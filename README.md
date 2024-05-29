# CITEseq_pipeline

Templates for the analysis of CITE-seq data in the Robinson Lab. 

This is a collection of scripts and analysis reports to run a CITE-seq analysis from `fastq` files up to differential expression/ pathway analysis.

This repo is **not a standardized** pipeline but a collection of scripts and reports that can be adapted to other data. 

The steps covered by the pipeline are : 

0) [Generate a count table from a FASTQ](#0-fastq-to-count-matrix)

1) (if HTO) [Diagnose the output HTO](#1-hto-diagnostic)

2) (if HTO) [Perform sample demultiplexing using HTO debarcoding](#2-hto-debarcoding)

3) [RNA QC and processing](#3-rna-diagnostics)

4) [ADT QC and processing](#4-antibodies-diagnostics-and-processing)

5) (optional) [ADT normalization using Seurat scaling]()

6) [Examine cell-types annotations](#6-look-at-cell-type-annotation)

7) [Examine the expression of markers among celltypes](#7-markers-expression-among-celltypes)

8) (optional) [Monocle pseudotime analysis](#8-optional-monocle-pseudotime-trajectory-analysis)

8)  [DS and DA analysis using pseudobulk](#9-differential-state-abundance-analysis)


-------------------------------------------------------

Steps 1 to X can be executed using the [`run_reports.sh`](run_reports.sh) master script.


## 0) FASTQ to count matrix

**Script**:

[alevin/alevin_run.sh](alevin/alevin_run.sh)

**Description**: 

Performs the mapping/ quantification for RNA, ADT and HTO (optional, depends on your experiment). 

**Adaptation**: 

- path to a salmon version: `salmon` 

- path to the folder with the FASTQ file: `FOLDER`

- number of samples: `SAMPLES`

- create an `alevin_ADT.csv` and `alevin_HTO.csv` (optional) index files with column 1 = ADT/HTO name, column 2 = barcodes. 


## 1) HTO diagnostic

**Script:**

[`analysis_rmd/1_HTO_diagnostics.Rmd`](analysis_rmd/1_HTO_diagnostics.Rmd)

**Description**: 

Run some diagnostics plots to ensure the quality of the HTO expression and that it can be used for debarcoding. 

In some cases where there are bad diagnostics for some samples, it might be worth to remove them. 

**Adaptation**: 

- The name of the samples : `fileName`

- The input folder: `inputFolder`

- The metadata file and the related variables that are used in the report: `metadataFile`

## 2) HTO debarcoding

**Script:**

[`analysis_rmd/2_HTO_debarcoding.Rmd`](analysis_rmd/2_HTO_debarcoding.Rmd)

**Description**: 

Performs sample demultiplexing on RNA/ ADT based on the HTO information. 

Requires HTO metadata in [`metadata/metadata.csv`](metadata/metadata.csv). 

**Adaptation**: 

- The name of the samples : `fileName`

- The input folder: `inputFolder`

- The metadata file,, HTO information, and the related variables that are used in the report: `metadataFile`

## 3) RNA diagnostics and processing

**Script**:

[`analysis_rmd/3_RNA_diagnostics.Rmd`](analysis_rmd/3_RNA_diagnostics.Rmd)


**Description**: 

Performs RNA QC and processing: 

- QC plots

- empty droplets

- RNA/ ADT/ cell filtering

- normalization

- Dimension reduction and plotting

- Clustering

- Quick automated annotation (singleR)


**Adaptation**:

- The metadata file `metadataFile`

- The folder containing the input data `folder`

- The reference used by singleR (`celldex::`)

- Parameters used for filtering

## 4) Antibodies diagnostics and processing 

**Script**:

[`analysis_rmd/4_Ab_diagnostics`](analysis_rmd/4_Ab_diagnostics.Rmd)


**Description**: 

Performs ADT QC and processing: 

- CLR normalization

- ADT expression QC 

- Dimension reduction based on ADTs

- Clustering based on ADTs


**Adaptation**:

- Resolution for the clustering

- Normalization method


## 5) (optional) ADT normalization using Seurat scaling


**Script**:

[`analysis_rmd/5_ab_seuratScaling.Rmd`](analysis_rmd/5_Ab_seuratScaling.Rmd)


**Description**: 

As an alternative to CLR normalization on ADT, performs Seurat normalization, followed by QC plots on the resulting expression.


**Adaptation**:

- Markers to plot to verify the validity of the normalization

- Metadata variables to plot

- Margin to normalize on (ADT-wise or cell-wise)


## 6) Look at cell-type annotation

**Script**:

[`analysis_rmd/6_use_celltypes.Rmd`](analysis_rmd/6_use_celltypes.Rmd)


**Description**: 

At this point, some automatic/ annotation should have be performed. It can be through `singleR` (step 3) or from manual annotation by looking at marker expression. The annotation used in this step should be stored in `metadata/annotation.csv`. 

This report look at the celltypes on reduced dimensions and at celltype-specific markers. Heatmaps of markers expression are saved in `reports/HEATMAPS`. 

**Adaptation**:

- Metadata variable names

- Processing parameters (e.g. normalization assay)



## 7) Markers expression among celltypes

**Script**:

[`analysis_rmd/7_markers_celltypes.Rmd`](analysis_rmd/7_markers_celltypes.Rmd)


**Description**: 

Examine cluster-specific marker expression based on ADT-defined celltypes. The results can help annotation or provide cluster-specific biomarkers. 

The marker discovery is based on `scran::findMarkers` and the report is divided into 2 parts: a 'lenient' and a 'stringent' approach to select the markers. 

The `Rmd` will also save tables of the celltype-specific markers in `reports/TABLES`. 

**Adaptation**:

- metadata variables: conditions, samples, etc. 

- processing parameters: which normalization to use for visualitation, etc. 

- parameters to find the markers.


## 8) (optional) Monocle pseudotime trajectory analysis

**Script**:

[`analysis_rmd/8_Monocle_pseudotime.Rmd`](analysis_rmd/8_Monocle_pseudotime.Rmd)


**Description**: 

Performs a pseudotime analysis using Monocle to examine the connectivity and trajectory between cell types. Also analyzes the population density across this pseudotime and how it changes across conditions. 

**Adaptation**:

- metadata variables: conditions, samples, etc. 

- Parameters to run Monocle 

## 9) Differential state/ abundance analysis

**Script**:

[`analysis_rmd/9_DS_DA_analysis.Rmd`](analysis_rmd/9_DS_DA_analysis.Rmd)


**Description**: 

Performs: 

- Differential state analysis using pseudobulk method (Muscat) 

- Differential abundance analysis

- Examine DS genes in Monocle pseudotime (optional, depending on step 8)

- Gene set enrichment (GO)

- Summarizes the DS genes found among the different contrasts

**Adaptation**:

- metadata variables: conditions to test, samples, etc. 

- Parameters to run the DS/ DA/ Gene enrichment analysis



