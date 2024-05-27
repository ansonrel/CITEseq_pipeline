args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

## Mandatory arguments
print(inputFolder)
print(rmdtemplate)
print(outputdir)
print(outputfile)


source("scripts/generate_report.R")

generateReport(inputFolder = inputFolder, fileName = fileName, 
               rmdTemplate = rmdtemplate, 
               outputDir = outputdir, outputFile = outputfile, ncores = ncores,
               forceOverwrite = TRUE, showCode = TRUE)
