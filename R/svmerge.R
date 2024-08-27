#' Checks if the start1, end1, start2, end2 are within the threshold
#' 
#' This is just a macro for checking if
#' (abs(start1 - start2) <= threshold & abs(end1 - end2) <= threshold)
#' is true or not. If it's true, we return True.
#' 
#' @param start1 position 1 for first sample
#' @param end1 position 2 for first sample
#' @param start2 position 1 for second sample
#' @param end2 position 2 for second sample
#' @param threshold how far away is acceptable for start1 and start2
within_threshold = function(start1, end1, start2, end2, threshold) {
  if (abs(start1 - start2) <= threshold & abs(end1 - end2) <= threshold) {
    return(T)
  } else {
    return(F)
  }
}

#' Goes through all SVs and assigns names to SV
#' 
#' If two SVs, x and y, in the sample_list have pos1 and pos2 within threshold of each
#' other, then the two SVs are given the same name. If a third SV, z, is within
#' threshold of y, but not x, it is still given the same name. Names are 
#' assigned sequentially based on the order they appear in sample_list
#' 
#' @param sample_list A list where the indices of the list correspond to the
#' sample. This can be generated using the make_sample_list function
#' @param threshold The threshold for which SVs should be called the same
#' @return A list including all_svs_per_sample, master_list, and sample_list
call_svs_per_sample = function(sample_list, threshold) {
  # Master list of names for all identified SVs
  master_sv = data.frame(name=NULL, start=NULL, end=NULL)
  
  # Counter for giving unique names to each identified SV
  sv_counter = 1 
  
  # Big list to keep track of what SVs were detected in which samples
  all_svs_per_sample = list()
  
  # Run through each tissue sample or cell-free sample
  for (list_name in names(sample_list)) {
    sample_sv = sample_list[[list_name]]
    # Running list of SVs detected in this sample
    sample_sv_names = c()
    # Add a name column to the sample_sv
    sample_list[[list_name]]$SV = "unassigned"
    # Run through all the SVs in sample
    for (i in 1:nrow(sample_sv)) {
      sample_pos1 = sample_sv$pos1[i]
      sample_pos2 = sample_sv$pos2[i]
      sample_chr1 = sample_sv$chr1[i]
      sample_chr2 = sample_sv$chr2[i]
      found = F
      if (nrow(master_sv) >= 1) {
        for (j in 1:nrow(master_sv)) {
          master_pos1 = master_sv$pos1[j]
          master_pos2 = master_sv$pos2[j]
          master_chr1 = master_sv$chr1[j]
          master_chr2 = master_sv$chr2[j]
          # If the sample matches something in the master list, take the name from the master list
          # Also double check that the chromosomes match
          if (within_threshold(sample_pos1, sample_pos2, master_pos1, master_pos2, threshold) && 
              sample_chr1 == master_chr1 &&
              sample_chr2 == master_chr2) {
            sv_name = master_sv$name[j]
            sample_list[[list_name]]$SV[i] = sv_name
            found = T
          } 
        }
      }
      # If the sample is not found in the master list, make a new name
      if (!found) {
        sv_name = sprintf("SV%04d", sv_counter)
        sample_list[[list_name]]$SV[i] = sv_name
        sv_counter = sv_counter + 1
      }
      sample_sv_names = c(sample_sv_names, sv_name)
      master_sv = rbind(master_sv, data.frame(pos1=sample_pos1, pos2=sample_pos2, chr1 = sample_chr1, chr2 = sample_chr2, name=sv_name, sample=list_name))
    }
    all_svs_per_sample[[list_name]] = sample_sv_names
  }
  ret = list(all_svs_per_sample = all_svs_per_sample, master_sv = master_sv, sample_list = sample_list)
  return(ret)
}

#' A function to create a table with all of the SVs within one patient
#' 
#' Takes a sample list, calls the call_svs_per_sample function, and then 
#' calls each SV as either clonal if it has cancer cell fraction (CCF) > 0.9 
#' (by default) or subclonal otherwise on a sample level.
#' On a patient-level, the SVs are called as the following.
#' Founder clone: clonal in all samples
#' Private clonal: clonal in exactly 1 sample, subclonal in 0.
#' Shared subclonal: clonal in no samples, subclonal in more than 1.
#' Private subclonal: clonal in no samples, subclonal in exactly 1 sample.
#' Shared clonal/subclonal: otherwise
#' 
#' @param sample_list A sample_list from the make_sample_list function. 
#' @param threshold The base pair threshold for calling SVs the same
#' @param CCF_cutoff The cutoff for CCF to consider an SV as clonal
make_merged_table = function(sample_list, threshold = 500, CCF_cutoff=0.9) {
  ret = call_svs_per_sample(sample_list, threshold=threshold)
  max_sv = max(as.integer(substring(ret$master_sv$name, 3)))
  df = data.frame(SV = sprintf("SV%04d", 1:max_sv))
  for (sample_name in names(ret$sample_list)) {
    sample = ret$sample_list[[sample_name]]
    CCF = rep(NA, max_sv)
    clonality = rep(NA, max_sv)
    sv_ids = as.integer(substring(sample$SV, 3))
    for (i in 1:nrow(sample)) {
      sv_id = sv_ids[i]
      CCF[sv_id] = sample$CCF[i]
      if (sample$CCF[i] >= CCF_cutoff) {
        clonality[sv_id] = "Clonal"
      } else {
        clonality[sv_id] = "Subclonal"
      }
    }
    df[paste0(sample_name, "CCF")] = CCF
    df[paste0(sample_name, "clonality")] = clonality
  }
  patient_clonality_func = function(i) {
    x = df[i, ]
    n_sample = (length(x) - 1) / 2
    if (sum(x == "Clonal", na.rm=T) == n_sample && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Founder Clone")
    } else if (sum(x == "Clonal", na.rm=T) > 1 && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Shared Clonal/Subclonal")
    } else if (sum(x == "Clonal", na.rm=T) == 1 && sum(x == "Subclonal", na.rm=T) == 0) {
      return("Private Clonal")
    } else if (sum(x == "Clonal", na.rm=T) > 0 && sum(x == "Subclonal", na.rm=T) > 0) {
      return("Shared Clonal/Subclonal")
    } else if (sum(x == "Subclonal", na.rm=T) > 1 && sum(x == "Clonal", na.rm=T) == 0) {
      return("Shared Subclone")
    } else if (sum(x == "Subclonal", na.rm=T) == 1 && sum(x == "Clonal", na.rm=T) == 0) {
      return("Private Subclone")
    } else {
      stop("THIS SHOULDNT HAPPEN")
    }
    return(sum(x == "AHHH THIS IS AN ERROR", na.rm=T) + sum(x == "Subclonal", na.rm=T))
  }
  df$patient_clonality = sapply(1:nrow(df), patient_clonality_func)
  
  all_samples = dplyr::bind_rows(ret$sample_list, .id = "column_label")
  classification = dplyr::group_by(all_samples, SV) 
  classification = dplyr::summarise(classification, min_pos=min(min(pos1), min(pos2)),
              max_pos=max(max(pos1), max(pos2)),
              chr1=paste(unique(chr1), collapse="|"),
              chr2=paste(unique(chr2), collapse="|"),
              classification_list=paste(classification, collapse="|"))
  df = merge(df, classification)
  l = list()
  l$merged_list = df
  l$annotated_samples = all_samples
  return(l)
}

#' Finds all the samples in the folder with sample_name and reads them into a sample_list format
#' 
#' The sample list format is a list of samples that binds together the
#' cluster certainty file which is output from SVClone, and the additional information
#' then makes the index of the list the name of the sample.
#' 
#' @param sample_ids A vector of the identifying names for each of the samples
#' @param cluster_certainty_fns A vector of the file names (fns) for the cluster
#' certainty files which are output by SVClone. For example
#' (IDENTIFIER)_cluster_certainty.txt
#' @param cluster_ccf_fns A vector of the file names (fns) for the cancer cell
#' fraction (CCF) files from SVClone. For example
#' (IDENTIFIER)_subclonal_structure_with_CCF_alan.txt
#' @param sv_info_fns A vector of the file names (fns) for the additional info
#' files from SVClone. For example
#' (IDENTIFIER)__filtered_svs.tsv
#' @return a list of samples that binds together the
#' cluster certainty file which is output from SVClone, and the additional information
#' then makes the index of the list the name of the sample.
make_sample_list = function(sample_ids, cluster_certainty_fns, cluster_ccf_fns, sv_info_fns) {
  if (length(cluster_certainty_fns) != length(cluster_ccf_fns) || length(cluster_ccf_fns) != length(sv_info_fns)) {
    stop("Make sure cluster_certainty_fns, cluster_ccf_fns, sv_info_fns all have the same length")
  } 
  sample_list = list()
  for (i in 1:length(cluster_certainty_fns)) {
    cluster_certainty = read.table(cluster_certainty_fns[i], header=T)
    sv_info = read.table(sv_info_fns[i], header=T)
    sv_info = sv_info[, 8:ncol(sv_info)]
    cluster_ccf = read.table(cluster_ccf_fns[i], header=T)
    sample_id = sample_ids[i]
    data = cbind(cluster_certainty, sv_info)
    data$CCF = cluster_ccf$CCF[data$most_likely_assignment]
    sample_list[[sample_id]] = data
  }
  return(sample_list)
}


#' Takes SVClone output from different tissues in the same patient and merges them
#' 
#' This function takes in the the outputs from SVClone and combines them. It 
#' requires that the function to generate the cancer cell fraction (CCF) and 
#' put it into the SVClone output files has already been run.
#' @export
write_merged_tables = function(sample_ids, patient_ids, cluster_certainty_fns, sv_info_fns, cluster_ccf_fns, merged_table_output_folder, threshold=500) {
  if (substr(merged_table_output_folder, nchar(merged_table_output_folder), nchar(merged_table_output_folder)) != '/') {
    merged_table_output_folder = paste0(merged_table_output_folder, "/")
  }
  if (length(sample_ids) != length(patient_ids) || 
      length(patient_ids) != length(cluster_certainty_fns) ||
      length(cluster_certainty_fns) != length(sv_info_fns) ||
      length(sv_info_fns) != length(cluster_ccf_fns)) {
    stop("All arguments must have the same length")
  }
  for (patient in unique(patient_ids)) {
      cur_cluster_certainty_fns = cluster_certainty_fns[patient_ids == patient]
      cur_sv_info_fns = sv_info_fns[patient_ids == patient]
      cur_cluster_ccf_fns = cluster_ccf_fns[patient_ids == patient]
      sample_list = make_sample_list(sample_ids, cur_cluster_certainty_fns, cur_cluster_ccf_fns, cur_sv_info_fns)
      ret = make_merged_table(sample_list, threshold=threshold)
      merged_list = ret$merged_list
      annotated_samples = ret$annotated_samples
      write.csv(merged_list, paste0(merged_table_output_folder, "matched_sv_", patient, ".csv"), row.names = F)
      write.csv(annotated_samples, paste0(merged_table_output_folder, "orig_sv_", patient, ".csv"), row.names = F)
  }
}


#' Function to create a clonality plot
#' 
#' Creates a plot where each point is a structural variant, and the x-axis is
#' the number of samples for which the SV was called subclonal, and the y-axis
#' is the number of samples for which the SV was called clonal
#' 
#' @param matched_sv_fns A vector of file names for the matched output. 
#' @export
clonality_plot = function(matched_sv_fns) {
  big_plot_df = data.frame()
  for (file in matched_sv_fns) {
    dat = read.csv(file)
    clonality = dat[, grepl("X.*clonality", names(dat))]
    clonal = rowSums(as.matrix(clonality) == "Clonal", na.rm = T)
    subclonal = rowSums(as.matrix(clonality) == "Subclonal", na.rm = T)
    n_sample = sum(grepl("X.*clonality", names(dat)))
    plot_df = data.frame(clonal, subclonal, n_sample, sample=tools::file_path_sans_ext(basename(file)))
    radius = 0.4
    noise1 = runif(min = -radius, max = radius, n=nrow(plot_df))
    noise2 = runif(min = -sqrt(radius^2 - noise1^2), max = sqrt(radius^2 - noise1^2), n=nrow(plot_df))
    plot_df$clonal = plot_df$clonal + noise1
    plot_df$subclonal = plot_df$subclonal + noise2
    
    big_plot_df = rbind(big_plot_df, plot_df)
  }
  
  n_sample_df = dplyr::group_by(big_plot_df, sample)
  n_sample_df = dplyr::summarise(n_sample_df, n_sample = unique(n_sample))
  circle_df = data.frame()
  for (i in 1:nrow(n_sample_df)) {
    n = n_sample_df$n_sample[i]
    sample = n_sample_df$sample[i]
    circles = expand.grid(0:n, 0:n)
    circles = dplyr::filter(circles, (Var1 + Var2) <= n)
    circles$sample = sample
    
    circle_df = rbind(circle_df, circles)
  }
  return(ggplot2::ggplot(big_plot_df) + ggplot2::geom_point(ggplot2::aes(x=subclonal, y=clonal), size = 0.5) + ggplot2::geom_abline(ggplot2::aes(intercept=n_sample, slope=-1), data=n_sample_df) + ggforce::geom_circle(ggplot2::aes(x0=Var1, y0=Var2, r=0.42), data=circle_df, color="#69bcf0") + ggplot2::facet_wrap(ggplot2::vars(sample)))
}

#' Function to motify the output of SVClone to include the CCF into text files
#' 
#' This function reads in files that are stored as .RData (_ccube_sv_results.RData)
#' and pulls out the cancer cell fraction (CCF) values and then appends them into
#' _subclonal_structure.txt files, then writes a new file with this information
#' at _subclonal_structure_with_CCF_alan.txt files. 
#' 
#' @param base_dir The directory where the subfolders are kept
#' @param samples A vector of identifiers for the samples which are used in the 
#' folder structure of SVClone 
#' @param test_run A boolean value which, if True, prevents the function from 
#' writing any files.
#' @param snv A boolean which, if True, tells the function to use the sv version 
#' of the files
add_ccf_to_output = function(base_dir, samples, test_run=T, snv=F) {
  clonal_props = c()
  snvstring = ""
  if (sv) {
    snvstring = "snvs/"
  }
  for (sample in samples) {
    print(paste0(base_dir, sample, "/ccube_out/", snvstring, sample,"_ccube_sv_results.RData"))
    
    # This gives you a file called doubleBreakPtsRes
    load(paste0(base_dir, sample, "/ccube_out/", snvstring, sample,"_ccube_sv_results.RData"))
    
    # This is the file that SVClone gave us that should have the cluster, the number of SVs, and the proportion
    sv_csv = read.table(paste0(base_dir, sample, "/ccube_out/", snvstring, sample,"_subclonal_structure.txt"), header = T)
    sv_csv$CCF = doubleBreakPtsRes$res$full.model$ccfMean
    
    # Quality control check that proportion/CCF is always the same
    tum_frac = sv_csv$proportion / sv_csv$CCF
    if (length(tum_frac) > 1) {
      first_ele_only = rep(tum_frac[1], length(tum_frac))
    } else {
      first_ele_only = tum_frac
    }
    if (!all.equal(tum_frac, first_ele_only)) {
      stop("Proportion should be the CCF times tumor fraction")
    }
    
    # Add the number of SVs that are clonal. We're defining this as CCF >= 0.9
    clonal_prop = sum(sv_csv$n_ssms[sv_csv$CCF >= 0.9])/sum(sv_csv$n_ssms)
    clonal_props = c(clonal_props, clonal_prop)
    if (!test_run) {
      write.table(sv_csv, file=paste0(base_dir, sample, "/ccube_out/", snvstring, sample,"_subclonal_structure_with_CCF_alan.txt"), row.names = F)
    } else {
      warning(paste("TEST RUN: ", paste0(base_dir, sample, "/ccube_out/", snvstring, sample,"_subclonal_structure_with_CCF_alan.txt","\nUse the test_run=F flag to write")))
    }
  }
}


