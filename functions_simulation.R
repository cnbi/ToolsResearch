########################## Functions for simulation #######################

# Collect results in a matrix---
collect_results <- function(design_matrix, results_folder, finding, pair, b = 1) {
    rows <-  seq(nrow(design_matrix))
    
    if (finding == "N2") {
        results_name <- "/ResultsN2Row"
        file_name <- "final_results_findN2"
    } else if (finding == "N1") {
        results_name <- "/ResultsN1Row"
        file_name <- "final_results_findN1"
    }
    
    if (pair == 2) { # When evaluating set 2 of hypotheses (H1vsH2)
        new_matrix <- matrix(NA, ncol = 7, nrow = nrow(design_matrix))
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            median.BF12 <- median(stored_result[[4]][, "BF.12"])
            median.BF21 <- median(stored_result[[4]][, "BF.21"])
            mean.PMP1 <- mean(stored_result[[4]][, "PMP.1"])
            mean.PMP2 <- mean(stored_result[[4]][, "PMP.2"])
            eta.BF12 <- stored_result$Proportion.BF12
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            new_matrix[row_design, ] <- c(median.BF12, mean.PMP1,
                                          median.BF21, mean.PMP2,
                                          eta.BF12, n2, n1)
        }
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF12",
                                  "mean.PMP1", "median.BF21", "mean.PMP2",
                                  "eta.BF12", "n2.final", "n1.final")
        saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, "_set2.RDS")))
        
    } else if (pair == 1) { # When evaluating set 1 of hypotheses (H0vsH1)
        new_matrix <- matrix(NA, ncol = 10, nrow = nrow(design_matrix))
        for (row_design in rows) {
            stored_result <- readRDS(paste0(results_folder, results_name, row_design, ".RDS"))
            median.BF01 <- median(stored_result[[b]][[6]][, "BF.01"])       # 6: data_H0
            median.BF10 <- median(stored_result[[b]][[7]][, "BF.10"])       # 7: data_H1
            mean.PMP0.H0 <- mean(stored_result[[b]][[6]][, "PMP.0"])        # 6: data_H0
            mean.PMP1.H0 <- mean(stored_result[[b]][[6]][, "PMP.1"])        # 6: data_H0
            mean.PMP0.H1 <- mean(stored_result[[b]][[7]][, "PMP.0"])        # 7: data_H1
            mean.PMP1.H1 <- mean(stored_result[[b]][[7]][, "PMP.1"])        # 7: data_H1
            n2 <- stored_result[[b]]$n2
            eta.BF01 <- stored_result[[b]]$Proportion.BF01
            eta.BF10 <- stored_result[[b]]$Proportion.BF10
            n1 <- stored_result[[b]]$n1
            new_matrix[row_design, ] <- c(median.BF01, mean.PMP0.H0,
                                                   mean.PMP1.H0, median.BF10, mean.PMP0.H1,
                                                   mean.PMP1.H1, eta.BF01, eta.BF10, n1, n2)
        }
        
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(names(design_matrix), "median.BF01",
                                           "mean.PMP0.H0", "mean.PMP1.H0", "median.BF10",
                                           "mean.PMP0.H1", "mean.PMP1.H1", "eta.BF01",
                                           "eta.BF10", "n1.final", "n2.final")
        saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, "_set1.RDS")))
    }
    return(new_matrix)
}

# Collect times in a matrix ----
collect_times <- function(design_matrix, rows = 0, pair, finding, results_folder) {
    rows <- ifelse(rows = 0, seq(nrow(design_matrix)), rows)
    new_matrix <- matrix(NA, nrow = nrow(design_matrix), ncol = 1)
    if (finding == "N2") {
        results_name <- "/timeN2Row"
        file_name <- "final_times_findN2"
    } else if (finding == "N1") {
        results_name <- "/timeN1Row"
        file_name <- "final_times_findN1"
    }
    # Pair of hypotheses to compare
    if (pair == 1) {
        type <- "_set1.RDS"
    } else if (pair == 2) {
        type <- "_set2.RDS"
    }
    
    for (row_result in rows) {
        stored_result <- readRDS(paste0(results_folder, results_name, row_result, ".RDS"))
        new_matrix[row_result, ] <- stored_result
    }
    new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
    colnames(new_matrix) <- c(names(design_matrix), "total.time")
    saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, type)))
    
    return(new_matrix)
}

# Run simulation ----
## Function to run the simulation for hypotheses set 2 -----
run_simulation_set2 <- function(Row, design_matrix, ndatasets, Max, batch_size, results_folder) {
    message("Starting simulation for row:", Row)
    Fixed <- design_matrix[Row, 6]
    Fixed <- toupper(Fixed)
    
    if (Fixed == "N1") {
        finding <- "N2"
    } else if (Fixed == "N2") {
        finding <- "N1"
    }
    
    # Start time
    start_time <- Sys.time()
    
    # Actual simulation
    if (Fixed == "N1") {
        ssd_results <- SSD_crt_inform(eff_size = design_matrix[Row, 2], 
                                      n1 = design_matrix[Row, 5], ndatasets = ndatasets,
                                      rho = design_matrix[Row, 1], 
                                      BF_thresh = design_matrix[Row, 3], eta =  design_matrix[Row, 4],
                                      fixed = as.character(design_matrix[Row, 6]), 
                                      max = Max,
                                      batch_size = batch_size)
    } else if (Fixed == "N2") {
        ssd_results <- SSD_crt_inform(eff_size = design_matrix[Row, 2], 
                                      n2 = design_matrix[Row, 5], ndatasets = ndatasets,
                                      rho = design_matrix[Row, 1], 
                                      BF_thresh = design_matrix[Row, 3], eta =  design_matrix[Row, 4],
                                      fixed = as.character(design_matrix[Row, 6]), 
                                      max = Max,
                                      batch_size = 1000)
    }
    
    # End time and save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0("Results", finding, "Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("time", finding, "Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    message("Finished row: ", Row)
    
    # Clean up memory
    rm(ssd_results)
    NULL
    gc()
}

# Function to run the simulation for hypotheses set 1 ----
run_simulation_set1 <- function(Row, design_matrix, ndatasets, Max, batch_size, results_folder, b) {
  message("Starting simulation for row:", Row)
  Fixed <- design_matrix[Row, 6]
  Fixed <- toupper(Fixed)
  
  if (Fixed == "N1") {
      finding <- "N2"
  } else if (Fixed == "N2") {
      finding <- "N1"
  }
  
  # Start time
  start_time <- Sys.time()
  
  # Actual simulation
  if (Fixed == "N1") {
    ssd_results <- SSD_crt_null(eff_size = design_matrix[Row, 2], 
                                n1 = design_matrix[Row, 5], ndatasets = ndatasets,
                                rho = design_matrix[Row, 1], 
                                BF_thresh = design_matrix[Row, 3], eta = design_matrix[Row, 4],
                                fixed = as.character(design_matrix[Row, 6]), b_fract = b, 
                                max = Max,
                                batch_size = batch_size)
  } else if (Fixed == "N2") {
    ssd_results <- SSD_crt_null(eff_size = design_matrix[Row, 2], 
                                n2 = design_matrix[Row, 5], n1 = 500,
                                ndatasets = ndatasets,
                                rho = design_matrix[Row, 1], 
                                BF_thresh = design_matrix[Row, 3], eta = design_matrix[Row, 4],
                                fixed = as.character(design_matrix[Row, 6]), b_fract = b, 
                                max = Max,
                                batch_size = batch_size)
  }
  
  # End time and save results
  end_time <- Sys.time()
  file_name <- file.path(results_folder, paste0("Results", finding, "Row", Row, ".RDS"))
  saveRDS(ssd_results, file = file_name)
  
  # Save running time
  running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
  time_name <- file.path(results_folder, paste0("time", finding, "Row", Row, ".RDS"))
  saveRDS(running_time, file = time_name)
  
  message("Finished row: ", Row)
  
  # Clean up memory
  rm(ssd_results)
  NULL
  gc()
}

# Check whether all conditions reached the power criterion -----
# Returns matrix with underpowered conditions and its results
reached_condition <- function(results_matrix, pair) {
    if (pair == 1) {
        unpowered <- results_matrix[!((results_matrix[, "eta.BF01"] > results_matrix[, "Probability"]) & (results_matrix[, "eta.BF10"] > results_matrix[, "Probability"])), ]
        
    } else if (pair == 2) {
        unpowered <- results_matrix[!(results_matrix[, "eta.BF12"] > results_matrix[, "Probability"]), ]
    }
    return(unpowered)
}

# Filter underpowered conditions: Return indexes -----
filter_underpowered <- function(results_matrix, pair){
    if (pair == 1) {
        run_again <- which(!((results_matrix[, "eta.BF01"] > results_matrix[, "Probability"]) & (results_matrix[, "eta.BF10"] > results_matrix[, "Probability"])))
    } else if (pair == 2) {
        run_again <- which(!(results_matrix[, "eta.BF12"] > results_matrix[, "Probability"]))
    }
    return(run_again)
}

