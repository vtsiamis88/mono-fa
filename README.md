### MONO FA: Factor analysis for absolute protein quantification in single sample SWATH-MS experiments.

### What is MONO FA?

MONO FA is an R tool implemented to perform absolute protein quantification for single sample SWATH-MS experiments. The tool uses a scoring schema that combines both the variation and co-variation of ions to assess their quality, filters the ion measurements and summarizes the absolute protein expressions. The quantification report is exported in .xlsx format.

### Input

The tool accepts ion peak area or intensity reports in .xlsx format.

(A template of the input can be found in Data/NIST data/NIST mAb_40min_2000ppm_8ug_3xIDA__FDR_4xSWexpep.xlsx)

### Installation and use
Download the present repository to your computer and install the following R packages:
```R
install.packages(c("readxl", "writexl", "crayon", "seqinr"))
```
The tool's functions are in MONO_FA_Functions.R file.

Open Perform_MONO_FA.R file and run the "source("MONO_FA_Functions.R")" to include the functions.

Run the absolute.mono.fa() function.

Instructions on how to use the function are provided by comments in the Perform_MONO_FA.R file.

### Contact
For software issues and general questions, please submit an issue.

### License
Apache-2.0
