#Alphalyse Project - Absolute protein quantification assisted by factor analysis
#June 2020
###################################
#Functions for MONO FA
###################################

#install.packages("readxl")
library(readxl)
#install.packages("crayon")
library(crayon)
#install.packages("writexl")
library(writexl)
#install.packages("seqinr")
library(seqinr)

absolute.mono.fa <- function(file, sheet, starting_column, normalization_method, fasta, standard_proteins, standards_ppm, ppm_factor, concentration, signal_multiplier, score_threshold, cv_threshold, max_na, hcp_prefixes){
  
  
  input_ions <- import.And.Preprocess.Ions(file = file, sheet = sheet, starting_column = starting_column, normalization_method = normalization_method)
  
  cat(crayon::green("__________________________________________________\n"))
  cat(crayon::green("1. Ion data successfully imported and normalized.\n"))
  cat(crayon::green("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"))
  
  data_FA <- perform.Factor.Analysis(data = input_ions$Ions, labels = input_ions$Labels, signal_multiplier = signal_multiplier)
  
  cat(crayon::green("__________________________________________________\n"))
  cat(crayon::green("2. Factor analysis is complete.\n"))
  cat(crayon::green("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"))
  
  data_summary <- summarize.Proteins(data = data_FA, labels = input_ions$Labels, score_threshold = score_threshold, cv_threshold = cv_threshold, max_na = max_na)
  
  cat(crayon::green("__________________________________________________\n"))
  cat(crayon::green("3. Protein abundance summarization is complete.\n"))
  cat(crayon::green("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"))
  
  data_standardized <- standardize.Proteins(data = data_summary, labels = input_ions$Labels, ppm_factor = ppm_factor, concentration = concentration)
  
  cat(crayon::green("__________________________________________________\n"))
  cat(crayon::green("4. Protein standardization is complete.\n"))
  cat(crayon::green("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"))
  
  output <- export.Report(proteins = data_standardized, ions = data_summary$Ions, file = file, fasta = fasta)
  
  cat(crayon::green("__________________________________________________\n"))
  cat(crayon::green("5. Analysis report is exported successfully.\n"))
  cat(crayon::green("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"))
  
  return(output)
}

import.And.Preprocess.Ions <- function(file, sheet = 1, starting_column, normalization_method = "median"){
  
  if(missing(file)){
    
    stop(crayon::bgBlack(crayon::white("Error: file name argument is missing!")), call. = FALSE)
    
  }
  
  #Check if the file can be imported
  error <- try(as.data.frame(readxl::read_xlsx(file, guess_max = 5000, sheet = sheet)), silent = TRUE)
  
  if(class(error) != "try-error"){
    
    ions <- as.data.frame(readxl::read_xlsx(file, guess_max = 5000, sheet = sheet))
    
  } else {
    
    stop(crayon::bgBlack(crayon::white("Error: the file path or the selected excel sheet is not correct! Input a valid file or sheet argument.")), call. = FALSE)
    
  }
  
  if(!all(apply(ions[,starting_column:ncol(ions)], 2, class) == "numeric")){
    
    stop(crayon::bgBlack(crayon::white("Error: the selected columns are not numeric! Input a valid starting_column argument.")), call. = FALSE)
    
  }
  
  if( (ncol(ions) - starting_column + 1) < 3 ){
    
    stop(crayon::bgBlack(crayon::white("Error: at least 3 replicates are required.")), call. = FALSE)
  
  } else if ( (ncol(ions) - starting_column + 1) == 3 ){
    
    #Additional column for the average of each ion if the replicates are 3
    ions$Average <- rowMeans(ions[,starting_column:ncol(ions)], na.rm = T)
    ions$Average[is.nan(ions$Average)] <- NA
    
  }
  
  #Rescale each replicate to the total intensity
  df <- data.frame(t(t(ions[,starting_column:ncol(ions)])/apply(ions[,starting_column:ncol(ions)], 2, sum, na.rm = TRUE)))
  
  #Log-2 transformation
  df <- log2(df)
  
  #Zero center normalization of each replicate
  if(normalization_method == "median"){
    
    df <- data.frame(t(t(df) - apply(df, 2, median, na.rm = TRUE)))
    
  } else if(normalization_method == "average") {
    
    df <- data.frame(t(t(df) - apply(df, 2, mean, na.rm = TRUE)))
    
  } else {
    
    stop(crayon::bgBlack(crayon::white("Error: invalid normalization method! Input a valid normalization method argument (median or average).")), call. = FALSE)
    
  }
  
  if(any(colnames(ions) == "Average")){
    
    colnames(df) <- c(paste0("Rep_", 1:(ncol(df) - 1)), "Rep_Average")
    
  } else {
    
    colnames(df) <- c(paste0("Rep_", 1:ncol(df)))
    
  }
  
  ions <- cbind(ions, df)
  
  labels <- data.frame(matrix(colnames(ions)[starting_column:ncol(ions)], ncol = 2, byrow = F))
  colnames(labels) <- c("Absolute", "Log")
  
  return(list(Ions = ions, Labels = labels))
  
}

fast.Farms <- function(probes, weight, mu, max_iter, force_iter, min_noise, fill_nan){
  
  if(missing(weight)){
    
    weight <- 0.5
    
  }	
  
  if(missing(mu)){
    
    mu <- 0
    
  }	
  
  if(missing(max_iter)){
    
    max_iter <- 1000
    
  }	
  
  if(missing(force_iter)){
    
    force_iter <- FALSE
    
  }
  
  if(missing(min_noise)){
    
    min_noise <- 0.0001
    
  }
  
  if(missing(fill_nan)){
    
    fill_nan <- 0
    
  }
  
  readouts <- as.matrix(probes)
  
  readouts[is.na(readouts)] <- fill_nan
  
  #Normalize and transform X
  X <- t(readouts)
  X <- t(t(X) - colMeans(X, na.rm = T))
  xsd <- apply(X, 2, function(x) sd(x, na.rm = T) * sqrt((length(x) - 1) / length(x)))
  xsd[xsd < min_noise] <- 1
  X <- t(t(X)/xsd)
  X[!is.finite(X)] <- 0
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  C <- crossprod(X, X)/n_samples
  
  #Positive definite
  C <- (C+t(C))/2
  C[which(C < 0)] <- 0
  
  #Robustness
  SVD <- svd(C)
  U <- SVD$u
  s <- SVD$d
  V <- t(SVD$v)
  s[s<min_noise] <- min_noise
  C <- U %*% diag(s) %*% V
  diag(C) <- abs(diag(C))
  
  #Initiation
  lamda <- sqrt(0.75*diag(C))
  psi <- diag(C) - lamda^2
  old_psi <- psi
  alpha <- weight * n_features
  E <- 1
  
  for(i in 1:max_iter){
    #E Step
    phi <- (1/psi)*lamda
    a <- as.vector(1+crossprod(lamda,phi))
    eta <- phi/a
    zeta <- C %*% eta
    E <- 1 - as.vector(eta) %*% lamda + as.vector(eta) %*% zeta
    
    #M Step
    lamda = zeta/(c(E) + as.vector(psi)*alpha)
    psi <- diag(C) - as.vector(zeta)[1] * lamda + psi * alpha * lamda * (mu - lamda)
    psi[psi < min_noise^2] <- min_noise^2
    
    if(!force_iter){
      
      if(max(abs(psi-old_psi))/max(abs(old_psi)) < min_noise/10){
        
        break
        
      }
    }
    
    old_psi <- psi
    
  }
  
  loading <- as.vector(sqrt(E))*lamda
  phi <- loading/psi
  weights <- loading/max(loading)
  weights <- round(weights,5)
  noise <- 1/as.vector(1+crossprod(loading,phi))
  
  loading.noise <- list("loadings" = weights, "noise" = noise)
  
  return(loading.noise)
  
}

perform.Factor.Analysis <- function(data, labels, standard_proteins, signal_multiplier){
  
  if(missing(data)){
    
    stop(crayon::bgBlack(crayon::white("Error: data argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(labels)){
    
    stop(crayon::bgBlack(crayon::white("Error: missing labels argument!")), call. = FALSE)
    
  }
  
  if(missing(standard_proteins)){
    
    standard_proteins <- c("a01|P68082|MYG_HORSE",
                           "a02|G5E5H7|LACB_BOVIN",
                           "a03|P00698|LYSC_CHICK",
                           "a09|P01966|HBA_BOVIN",
                           "a10|P02070|HBB_BOVIN",
                           "a12|P00709|LALBA_HUMAN",
                           "a13|P02787|TRFE_HUMAN")
    
  }
  
  #Check which standard proteins are not found in the data
  not_found <- setdiff(standard_proteins, unique(data[,1]))
  
  if(length(not_found) != 0){
    
    not_found_accessions <- unlist(lapply(not_found, function(x) unlist(strsplit(x = x, split = "\\|"))[2]))
    
    not_found_ids <- lapply(not_found_accessions, function(x) which(grepl(pattern = x, x = data[,1])) )
    
    if(sum(lengths(not_found_ids) == 0) > 0){
      
      stop(crayon::bgBlack(crayon::white("Error: not all standard proteins are found in the ion report.")), call. = FALSE)
      
    }
    
    for (i in 1:length(not_found)) {
      
      data[unlist(not_found_ids[[i]]), 1] <- not_found[i]
      
    }
    
  }
  
  #Create permutations of replicates and signal patterns
  total_combinations <- expand.grid(lapply(1:nrow(labels), function(x) 1:nrow(labels)))
  permutations <- total_combinations[apply(total_combinations, 1, function(x) {(length(unique(x)) == nrow(labels)) & 
                                                                               (!all(x == 1:nrow(labels))) & 
                                                                               (!all(x == nrow(labels):1))}),]
  
  permutations <- permutations[sample(1:nrow(permutations), size = 3, replace = F), ]
  permutations <- t(cbind(t(data.frame(matrix(c(1:nrow(labels), nrow(labels):1), nrow = 2, byrow = T))), t(permutations)))
  
  signal <- (1:nrow(labels) * signal_multiplier) - mean(1:nrow(labels) * signal_multiplier)
  
  signals <- t(apply(permutations, 1, function(x){signal[x]}))
  
  #Find the unique proteins in data and the indices of corresponding ions
  unique_proteins <- unique(data[,1])

  ions_per_protein_ids <- lapply(unique_proteins, function(x) which(data[,1] == x))
  ions_per_protein_data <- lapply(ions_per_protein_ids, function(x) data[x,])

  for (i in 1:length(ions_per_protein_data)) {

    weights_df <- data.frame(matrix(NA, ncol = nrow(signals), nrow = nrow(ions_per_protein_data[[i]])))
    colnames(weights_df) <- paste0("Weights_", 1:nrow(signals))

    if(length(unique(ions_per_protein_data[[i]][,2])) > 1){

      if(any(standard_proteins == ions_per_protein_data[[i]][1,1])){

        prot_df <- ions_per_protein_data[[i]][,labels$Log]
        non_na <- which(colSums(apply(prot_df, 1, is.na)) == 0)

        if(length(non_na) < nrow(prot_df)){

          probe <- prot_df[non_na, ]

        } else {

          probe <- prot_df

        }
        
        if(!is.vector(probe)){

          if(nrow(probe) >= 2){

            prot_farms <- fast.Farms(probe, weight = 0.9, mu = 0)
            weights_df[non_na,1] <- prot_farms$loadings

          }
        }

      } else {

        for (j in 1:nrow(signals)) {

          prot_df <- ions_per_protein_data[[i]][,labels$Log]
          ion_sd <- sd(as.numeric(unlist(prot_df)), na.rm = T)
          prot_df <- t(apply(prot_df, 1, function(x){x + (signals[j,] * ion_sd)}))

          non_na <- which(colSums(apply(prot_df, 1, is.na)) == 0)

          if(length(non_na) < nrow(prot_df)){

            probe <- prot_df[non_na, ]

          } else {

            probe <- prot_df

          }
          
          if(!is.vector(probe)){
            
            if(nrow(probe) >= 2){

              prot_farms <- fast.Farms(probe, weight = 0.9, mu = 0)
              weights_df[non_na,j] <- prot_farms$loadings

            }
          }

        }

      }

    }

    prot_df <- cbind(ions_per_protein_data[[i]], weights_df)

    if(any(standard_proteins == ions_per_protein_data[[i]][1,1])){

      prot_df$Average_norm <- rowMeans(2^ions_per_protein_data[[i]][, setdiff(labels$Log[2:nrow(labels)], "Rep_Average")], na.rm = T)
      prot_df$Average_norm[is.nan(prot_df$Average_norm)] <- NA
      prot_df$SD_norm <- apply(2^ions_per_protein_data[[i]][, setdiff(labels$Log[2:nrow(labels)], "Rep_Average")], 1, sd, na.rm = T)
      prot_df$CV_norm <- prot_df$SD_norm/prot_df$Average_norm


    } else {

      prot_df$Average_norm <- rowMeans(2^ions_per_protein_data[[i]][, setdiff(labels$Log, "Rep_Average")], na.rm = T)
      prot_df$Average_norm[is.nan(prot_df$Average_norm)] <- NA
      prot_df$SD_norm <- apply(2^ions_per_protein_data[[i]][, setdiff(labels$Log, "Rep_Average")], 1, sd, na.rm = T)
      prot_df$CV_norm <- prot_df$SD_norm/prot_df$Average_norm

    }

    prot_df$Min_weight <- apply(prot_df[, colnames(weights_df)], 1, function(x){ if(all(is.na(x))){ return(NA) } else { return(min(x, na.rm = T))} })
    prot_df$Inverse_weight <- 1 - prot_df$Min_weight
    prot_df$Score <- prot_df$Inverse_weight * prot_df$CV_norm

    ions_per_protein_data[[i]] <- prot_df

  }

  return(ions_per_protein_data)
  
}

summarize.Proteins <- function(data, labels, score_threshold = 0.05, cv_threshold = 0.5, max_na){
  
  if(missing(data)){
    
    stop(crayon::bgBlack(crayon::white("Error: data argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(labels)){
    
    stop(crayon::bgBlack(crayon::white("Error: labels argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(max_na)){
    
    stop(crayon::bgBlack(crayon::white("Error: max_na argument is missing!")), call. = FALSE)
    
  }
  
  proteins <- vector(mode = "list", length = length(data))
  
  for (i in 1:length(data)) {
    
    #Filter ions based on Score
    data[[i]]$Accept <- data[[i]]$Score <= score_threshold
    
    #Filter ions based on number of missing values
    na_accept <- apply(data[[i]][is.na(data[[i]]$Score), setdiff(labels$Absolute, "Average")], 1, function(x){ sum(is.na(x)) }) <= max_na
    
    #Filter ions based on CV when Score does not exist
    cv_accept <- data[[i]]$CV_norm[is.na(data[[i]]$Score)] < cv_threshold
    cv_accept[is.na(cv_accept)] <- FALSE
    
    data[[i]]$Accept[is.na(data[[i]]$Score)] <- na_accept & cv_accept
    
    #Summarize protein expression if there are more than one unique peptides
    if(length(unique(data[[i]][,2])) != 1 & sum(data[[i]]$Accept) > 1 ){
      
      prot <- aggregate(data[[i]][data[[i]]$Accept, setdiff(labels$Absolute, "Average")], by = list(data[[i]][data[[i]]$Accept, 1]), FUN = sum, na.rm = T)
      colnames(prot)[1] <- "Protein"
      prot$Peptides <- length(unique(data[[i]][data[[i]]$Accept, 2]))
      prot$Removed_Peptides <- length(unique(data[[i]][ ,2])) - prot$Peptides
      prot$Ions <- sum(data[[i]]$Accept)
      prot$Removed_Ions_By_FA <- sum(data[[i]]$Score > score_threshold, na.rm = T)
      prot$Total_Removed_Ions <- nrow(data[[i]]) - prot$Ions
      prot <- prot[,c("Protein", "Peptides", "Removed_Peptides", "Ions", "Removed_Ions_By_FA", "Total_Removed_Ions", setdiff(labels$Absolute, "Average"))]
      
    } else {
      
      prot <- NA
      
    }
    
    proteins[[i]] <- prot
    
  }
  
  return(list(Ions = data, Proteins = proteins))
  
}

standardize.Proteins <- function(data, standard_proteins, standards_ppm, ppm_factor, concentration, labels, hcp_prefixes){
  
  if(missing(data)){
    
    stop(crayon::bgBlack(crayon::white("Error: data argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(standard_proteins)){
    
    standard_proteins <- c("a01|P68082|MYG_HORSE",
                           "a02|G5E5H7|LACB_BOVIN",
                           "a03|P00698|LYSC_CHICK",
                           "a09|P01966|HBA_BOVIN",
                           "a10|P02070|HBB_BOVIN",
                           "a12|P00709|LALBA_HUMAN",
                           "a13|P02787|TRFE_HUMAN")
    
  }
  
  if(missing(standards_ppm)){
    
    standards_ppm <- c(2014.81481481481,
                       2017.77777777778,
                       2007.93650793651,
                       1954.60317460317,
                       2059.04761904762,
                       2017.46031746032,
                       2001.0582010582)
  }
  
  if(missing(hcp_prefixes)){
    
    hcp_prefixes <- c("sp", "tr")
    
  }
  
  standards_ppm <- standards_ppm * ppm_factor * concentration
  
  protein_names <- unlist(lapply(data$Ions, function(x) x[1,1]))
  
  #Standard proteins
  standard_indices <- sapply(standard_proteins, function(x) which(protein_names == x))
  
  standards <- do.call(rbind, data$Proteins[standard_indices])
  standards$Average <- rowMeans(standards[, setdiff(labels$Absolute[2:nrow(labels)], "Average")], na.rm = T)
  standards$CV <- apply(standards[, setdiff(labels$Absolute[2:nrow(labels)], "Average")], 1, sd, na.rm = T)/standards$Average
  standards$Ratio <- standards$Average/standards[,labels$Absolute[1]]
  standards$ppm <- standards_ppm
  standards$Reference <- standards$Average/standards$ppm
  standards$Calculated_ppm <- standards$Average/median(standards$Reference, na.rm = T)
  standards$Accuracy <- 100*standards$Calculated_ppm/standards$ppm
  rownames(standards) <- 1:nrow(standards)
  
  #Continue here ################################
  
  #Standard impurity proteins
  standard_impurity_indices <- which(rowSums(sapply(hcp_prefixes, function(x) grepl(x, protein_names))) == 0)
  standard_impurity_indices <- setdiff(standard_impurity_indices, standard_indices)
  
  if(length(standard_impurity_indices > 0)){
    
    standard_impurities <- do.call(rbind, data$Proteins[standard_impurity_indices])
    
    standard_impurities <- standard_impurities[rowSums(is.na(standard_impurities)) != ncol(standard_impurities),]
    standard_impurities$Average <- rowMeans(standard_impurities[,setdiff(labels$Absolute, "Average")], na.rm = T)
    standard_impurities$CV <- apply(standard_impurities[,setdiff(labels$Absolute, "Average")], 1, sd, na.rm = T)/standard_impurities$Average
    standard_impurities$Ratio <- standard_impurities$Average/standard_impurities[,labels$Absolute[1]]
    standard_impurities$ppm <- standard_impurities$Average/median(standards$Reference, na.rm = T)
    rownames(standard_impurities) <- 1:nrow(standard_impurities)
    
  }
  
  #Host cell proteins
  hcp_indices <- setdiff(1:length(protein_names), c(standard_indices, standard_impurity_indices))
  
  if(length(hcp_indices) > 0){
    
    hcps <- do.call(rbind, data$Proteins[hcp_indices])
    hcps <- hcps[rowSums(is.na(hcps)) != ncol(hcps),]
    hcps$Average <- rowMeans(hcps[,setdiff(labels$Absolute, "Average")], na.rm = T)
    hcps$CV <- apply(hcps[,setdiff(labels$Absolute, "Average")], 1, sd, na.rm = T)/hcps$Average
    hcps$Ratio <- hcps$Average/hcps[,labels$Absolute[1]]
    hcps$ppm <- hcps$Average/median(standards$Reference, na.rm = T)
    rownames(hcps) <- 1:nrow(hcps)
    
  }
  
  return(list(Standards = standards, Standard_Impurities = standard_impurities, HCPS = hcps))
  
}

export.Report <- function(proteins, ions, file, fasta){
  
  if(missing(proteins) | missing(ions)){
    
    stop(crayon::bgBlack(crayon::white("Error: protein or ion data argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(file)){
    
    stop(crayon::bgBlack(crayon::white("Error: file argument is missing!")), call. = FALSE)
    
  }
  
  if(missing(fasta)){
    
    stop(crayon::bgBlack(crayon::white("Error: fasta argument is missing!")), call. = FALSE)
    
  }
  
  #Check if the file can be imported
  error <- try(seqinr::read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, legacy.mode = TRUE), silent = TRUE)
  
  if(class(error) != "try-error"){
    
    fasta <- seqinr::read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, legacy.mode = TRUE)
    
  } else {
    
    stop(crayon::bgBlack(crayon::white("Error: the fasta file path is not correct!")), call. = FALSE)
    
  }
  
  #Prepare STDS_FA sheet
  standards_sheet <- do.call(rbind, list(proteins$Standards, NA, NA))
  column_diff <- setdiff(colnames(standards_sheet), colnames(proteins$Standard_Impurities))
  to_add <- data.frame(matrix(NA, nrow = nrow(proteins$Standard_Impurities), ncol = length(column_diff) ))
  colnames(to_add) <- column_diff
  standard_impurities_add <- cbind(proteins$Standard_Impurities, to_add)
  standards_sheet <- rbind(standards_sheet, standard_impurities_add)
  
  #Prepare Ions_Area_FA sheet
  ions_sheet <- do.call(rbind, ions)
  
  #Prepare HCP_FA sheet
  hcps_sheet <- proteins$HCPS
  hcps_sheet <- hcps_sheet[order(hcps_sheet$ppm, decreasing = T),]
  rownames(hcps_sheet) <- 1:nrow(hcps_sheet)
  
  p_annotations <- sub(".+? ", "", unlist(seqinr::getAnnot(fasta)))
  p_names <- names(fasta)
  fasta <- data.frame(Name = p_names, Annotation = p_annotations)
  annotation_ids <- sapply(hcps_sheet$Protein, function(x) which(fasta$Name == x))
  hcps_sheet$Description <- fasta$Annotation[annotation_ids]
  
  #Output path
  file_split <- unlist(strsplit(file, split = "\\."))
  file <- paste0(file_split[1], "_FA.", file_split[2], collapse = "", sep = "")
  
  writexl::write_xlsx(x = list(Ions_Area_FA = ions_sheet, HCP_FA = hcps_sheet, STDS_FA = standards_sheet), path = file)
  
  return(list(Ions = ions_sheet, HCP = hcps_sheet, STDS = standards_sheet))
  
}


