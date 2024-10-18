setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("MONO_FA_Functions.R")

#Example analysis
#####################################################################################

#Function arguments
#
#file: the path to the excel file that contains the ion intensities.
#
#sheet: the sheet number in the file that corresponds to the ions (1 means first sheet, 2 means second sheet).
#
#starting_column: the column in the excel sheet that the ion intensities start. In principle it should be the column after "Residue" column.
#
#normalization_method: the method used to normalize ion intensities per replicate. It could be either by "median" or "average".
#
#fasta: the path to the fasta file.
#
#standard_proteins: a string vector with the names of the protein standards. With the standard analysis in Alphalyze the standards are always
#                   the 7 proteins. This argument should be used only if new standards are used.
#
#standards_ppm: a numeric vector that correspond to the known ppm of the proteins in standard_proteins. This argument should be used only
#               if new standards are used.
#
#ppm_factor: the factor that determines the amount of the standard proteins in the sample for the formula: ppm * concentration * ppm_factor.
#
#concentration: the concentration of the standard proteins.
#
#signal_multiplier: a multiplier that determines the signal used during the factor analysis based on the standard deviation of the ions per protein (probe).
#
#score_threshold: ions with score higher than the score_threshold value are discarded. The score_threshold is applied only to peptides without any missing value.
#
#cv_threshold: cv_threshold is used for ions with at least one missing value. Ions with CV values higher than the threshold are not considered
#              for the protein summarization.
#
#max_na: assists a missing value filtering. Ions with higher number of values than max_na are not considered in the summarization.

#####################################################################################
##### DATASET NIST data (Example)
#####################################################################################

result_NIST <- absolute.mono.fa(file = "Data/NIST data/NIST mAb_40min_2000ppm_8ug_3xIDA__FDR_4xSWexpep.xlsx",
                                sheet = 1, 
                                starting_column = 10, 
                                normalization_method = "median", 
                                fasta = "Data/NIST data/A0046_Scanning SWATH samples_Swissprot_22-Jun-2020.fasta",
                                ppm_factor = 1, 
                                concentration = 1, 
                                signal_multiplier = 0.75, 
                                score_threshold = 0.05, 
                                cv_threshold = 0.5, 
                                max_na = 3)

#Access Ions
NIST_ions <- result_NIST$Ions

#Access HCP protein expressions
NIST_HCP <- result_NIST$HCP

#Access STD protein expressions
NIST_STD <- result_NIST$STDS

#####################################################################################




