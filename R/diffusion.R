CalcLigandDist = function(lidardata, dist_threshold = 10, distance_bin = NULL, RealDist = F, imp = T, include_self = TRUE){
  if (isTRUE(RealDist)) stop("RealDist = TRUE is not implemented; use distance_bin.")
  if (length(dist_threshold) != 1L || !is.finite(dist_threshold) || dist_threshold < 1 || dist_threshold != as.integer(dist_threshold)) stop("dist_threshold must be a positive integer.")
  dist_threshold = dist_threshold + 1

  pos = lidardata@Image@pos
  if (imp) {
    ligand_exp = lidardata@ImputationData@ligand_exp
  } else {
    ligand_exp = lidardata@Data@ligand_exp
  }
  ligandname = lidardata@DB@ligandname
  cell_data = lidardata@cellname
  nligand = lidardata@DB@nligand
  ncell = lidardata@ncell

  dist_pos = as.matrix(dist(pos)) %>% round(3)
  if (!is.null(distance_bin)) {
    bin <- distance_bin
  } else {
    bin <- min(dist_pos[dist_pos > 0])
  }
  if (length(bin) != 1L || !is.finite(bin) || bin <= 0) stop("distance_bin must be positive and finite.")
  dist_pos_bin_all = dist_pos%/%bin + 1
  dist_pos_bin_all = Matrix::Matrix(dist_pos_bin_all, sparse = T)

  dist_pos_bin = dist_pos_bin_all
  dist_pos_bin[dist_pos_bin>dist_threshold] = 0
  dist_pos_bin = Matrix::Matrix(dist_pos_bin, sparse = T)

  ligand_dist = array(0, dim = c(nligand, ncell, dist_threshold))

  #250826 change
  start_bin = ifelse(include_self, 1, 2)
  for (num in start_bin:dist_threshold) {
    binary_matrix <- Matrix((dist_pos_bin == num), sparse = TRUE)
    gene_sums = as.matrix(ligand_exp %*%  binary_matrix)
    ligand_dist[, , num] = gene_sums
  }
  #250826 change over


  dimnames(ligand_dist) = list(
    Ligand = ligandname,
    Cell = cell_data,
    Distance = 1:dist_threshold
  )

  lidardata@Image@dist = dist_pos
  lidardata@Image@dist_bin_all = dist_pos_bin_all

  lidardata@Diffusion@ligand_dist = ligand_dist
  lidardata@Diffusion@dist_threshold = dist_threshold
  lidardata@Diffusion@bin = bin
  lidardata@Diffusion@include_self = include_self

  return(lidardata)
}

#Solve diff----
## Wrapped function to optimize ligand diffusion----
OptLigandDiff <- function(lidardata, CalcLigandConc = TRUE, imp = TRUE, ligand = NULL) {



  ligand_act = lidardata@Diffusion@ligand_act
  ligand_dist = lidardata@Diffusion@ligand_dist
  dist_threshold = lidardata@Diffusion@dist_threshold
  if (is.null(ligand) || length(ligand) == 0) {
    ligandname = lidardata@DB@ligandname
  } else {
    ligandname = ligand
  }
  receptorname = lidardata@DB@receptorname
  lr_list = lidardata@DB@lr_list
  if (imp) {
    receptor_exp = lidardata@ImputationData@receptor_exp
  } else {
    receptor_exp = lidardata@Data@receptor_exp
  }
  nligand = lidardata@DB@nligand
  dist <- 0:(dist_threshold - 1)

  ligand_diff_list <- list()
  ligand_conc_list <- list()
  failed_ligand <- c()
  prog = (1:10)/10
  n_prog = 1

  for (i in seq_along(ligandname)) {
    current_ligand <- ligandname[i]
    #print(i)

    # Progress
    if (nligand > 20 && i / nligand >= prog[n_prog]) {
      message(sprintf("Processing %d%%", as.integer(prog[n_prog] * 100)))
      n_prog <- n_prog + 1
    }

    # Filter receptors for the current ligand
    receptor <- filter(lr_list, from %in% current_ligand) %>% dplyr::select(to) %>% unlist() %>% unique()

    # Get receptor expression
    if(length(receptor) > 1) {
      one_receptor_exp <- receptor_exp[match(receptor, receptorname),] %>% as.array() %>% colSums()
    } else {
      one_receptor_exp <- receptor_exp[match(receptor, receptorname),]
    }

    # Extract ligand activity and dissociation
    one_ligand_act <- ligand_act[, colnames(ligand_act) == current_ligand] %>% as.vector()
    one_ligand_dist <- ligand_dist[current_ligand, one_receptor_exp > 0, ]
    one_ligand_act <- one_ligand_act[one_receptor_exp > 0]

    # Check activity
    if (sum(one_ligand_act > 0) <= 10 || var(one_ligand_act) == 0) {
      failed_ligand <- c(failed_ligand, current_ligand)
      next
    }

    # Fit ligand model with timeout and error handling
    par <- tryCatch({
      withTimeout({
        FitLigandModel(current_ligand, one_ligand_act, one_ligand_dist, dist, dist_threshold)
      }, timeout = 3000)  # Timeout after 3000 seconds
    }, error = function(e) {
      message(sprintf("Ligand %s failed FitLigandModel: %s", current_ligand, e$message))
      return(NULL)
    }, TimeoutException = function(e) {
      message(sprintf("Ligand %s timed out (>3000 s)", current_ligand))
      return(NULL)
    })

    if (is.null(par)) {
      failed_ligand <- c(failed_ligand, current_ligand)
      next
    }

    # Compute dis and predict values
    one_ligand_diff <- 1 / (1 + ((dist / (dist_threshold - 1)) / par[["par"]][1])^par[["par"]][2])
    one_ligand_diff <- maxmin(one_ligand_diff)
    predict_value_r <- one_ligand_dist %*% one_ligand_diff

    if (sum(predict_value_r > 0) <= 10 || var(predict_value_r) == 0) {
      failed_ligand <- c(failed_ligand, current_ligand)
      next
    }

    # Get ligand diffusion info
    ligand_diff <- GetLigandDiff(current_ligand, one_ligand_act, predict_value_r, par)
    ligand_diff_list[[current_ligand]] <- ligand_diff

    # Calculate ligand concentration
    if (CalcLigandConc) {
      one_ligand_diff_nor <- maxmin(one_ligand_diff)
      one_ligand_dist_full <- ligand_dist[current_ligand, , ]
      predict_value_full <- one_ligand_dist_full %*% one_ligand_diff_nor
      ligand_conc <- PredLigandConc(predict_value_full)
      ligand_conc_list[[current_ligand]] <- ligand_conc
    }
  }

  # Combine results
  if (length(ligand_diff_list) > 0) {
    combined_ligand_diff_df <- do.call(rbind, ligand_diff_list)
  } else {
    combined_ligand_diff_df <- data.frame(Ligand = character(), ki = numeric(), n = numeric(), r = numeric(), p.value = numeric())
  }

  if (length(ligand_conc_list) > 0) {
    combined_ligand_conc <- do.call(cbind, ligand_conc_list)
    colnames(combined_ligand_conc) = names(ligand_conc_list)
  } else {
    combined_ligand_conc <- Matrix::Matrix(matrix(numeric(), nrow = lidardata@ncell, ncol = 0, dimnames = list(lidardata@cellname, NULL)), sparse = TRUE)
  }

  lidardata@Diffusion@ligand_diff = combined_ligand_diff_df
  lidardata@Diffusion@failed_ligand = as.character(failed_ligand)
  if (CalcLigandConc) lidardata@Diffusion@ligand_conc = combined_ligand_conc

  if (length(failed_ligand) > 0) {
    message("The following ligands failed calculations: ", paste(failed_ligand, collapse = ", "))
  } else {
    message("All ligands successfully calculated.")
  }


  return(lidardata)
}


## Fit bounded diffusion parameters----
FitLigandModel <- function(current_ligand, one_ligand_act, one_ligand_dist, dist, dist_threshold) {

  # Objective function
  objective_function <- function(par_value) {
    ki <- par_value[1]
    n <- par_value[2]
    one_ligand_diff <- 1 / (1 + ((dist / (dist_threshold - 1)) / ki)^n)
    one_ligand_diff = maxmin(one_ligand_diff)
    r <- -cor(one_ligand_act, as.vector(one_ligand_dist %*% one_ligand_diff))
    return(r)
  }
  # Change 0906 start
  # Initialize parameters
  par_initial <- c(0.01, 0.01)

  lower_bounds <- c(0.01, 0.01)
  upper_bounds <- c(1, 10)
  # Change 0906 over
  # Optimize parameters
  par <- optim(par = par_initial,
               fn = objective_function,
               method = "L-BFGS-B",
               lower = lower_bounds,
               upper = upper_bounds,
               control = list(trace = 0))
  #  par <- optim(par = par_initial,
  #               fn = objective_function,
  #               method = "SANN",
  #               control = list(maxit = 10000, temp = 10, tmax = 10))

  return(par)
}

## Function to get ligand diffusion information----
GetLigandDiff <- function(current_ligand, one_ligand_act, predict_value, par) {
  test <- cor.test(one_ligand_act, as.vector(predict_value))
  p.value <- test[["p.value"]]
  r <- test[["estimate"]]
  ligand_diff <- data.frame(Ligand = current_ligand, ki = par[["par"]][1], n = par[["par"]][2], r = r, p.value = p.value)
  return(ligand_diff)
}


## Function to calculate predictions----
PredLigandConc <- function(predict_value) {
  ligand_conc <- Matrix(predict_value, sparse = TRUE)
  return(ligand_conc)
}

maxmin = function(x){
  span <- max(x) - min(x)
  if (!is.finite(span) || span == 0) return(x * 0)
  (x-min(x))/span
}

#Lidarscore=====
