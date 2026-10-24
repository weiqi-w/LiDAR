CompLidarScore = function(lidardata, mode = c("aggregate","pairwise"),
                          ligand_list = NULL,pval = 0.01, imp = T, include_self = T) {
  mode <- match.arg(mode)
  switch(
    mode,
    aggregate = {
      ligand_act = lidardata@Diffusion@ligand_act
      if (imp) {
        ligand_exp = lidardata@ImputationData@ligand_exp
      } else {
        ligand_exp = lidardata@Data@ligand_exp
      }
      ligand_conc = lidardata@Diffusion@ligand_conc
      receptorname = lidardata@DB@receptorname
      lr_list = lidardata@DB@lr_list
      if (imp) {
        receptor_exp = lidardata@ImputationData@receptor_exp
      } else {
        receptor_exp = lidardata@Data@receptor_exp
      }
      nreceptor = lidardata@DB@nreceptor

      #succ_ligand = lidardata@Diffusion@ligand_diff %>% filter(p.value < pval) %>% pull(Ligand)
      if (!is.null(ligand_list)) {
        succ_ligand <- ligand_list
      } else {
        succ_ligand <- lidardata@Diffusion@ligand_diff %>%
          filter(p.value < pval) %>%
          pull(Ligand)
      }
      if (!all(succ_ligand %in% lidardata@Diffusion@ligand_diff$Ligand)) stop("Requested ligands must have fitted diffusion results.")
      if (!length(succ_ligand)) stop("No ligands pass the score selection; adjust pval or supply ligand_list.")
      n_succ_ligand = length(succ_ligand)

      ligand_act = ligand_act[,succ_ligand, drop = FALSE]
      ligand_act_divmax <- Matrix(apply(ligand_act, 2, divmax)*10, sparse = T)

      ligand_exp = ligand_exp[succ_ligand, , drop = FALSE]
      ligand_conc = ligand_conc[,succ_ligand, drop = FALSE]
      ligand_conc_divmax <- Matrix(apply(ligand_conc, 2, divmax)*10, sparse = T)


      #for (i in 1:n_succ_ligand){
      #ligand_exp_max = max(ligand_exp[i,])
      #ligand_exp_min = min(ligand_exp[i,])
      #  ligand_conc[,i] = 10*(ligand_conc[,i]- ligand_exp_min)/(ligand_exp_max-ligand_exp_min)
      #}



      #l_r adjacent----
      lr_match = matrix(0, nrow = n_succ_ligand, ncol = nreceptor,
                        dimnames = list(succ_ligand,receptorname))
      lr_index <- cbind(match(lr_list$from, succ_ligand), match(lr_list$to, receptorname))
      lr_index <- lr_index[stats::complete.cases(lr_index), , drop = FALSE]
      lr_match[lr_index] = 1
      lr_match = Matrix(lr_match, sparse = T)

      ##ligand_cores_receptor_exp----
      #l_r adjacent %*% receptor expr = The expr of receptor corresponding to the ligand
      ligand_cores_receptor_exp = t(lr_match %*% receptor_exp)
      ligand_cores_receptor_exp_divmax <- Matrix(apply(ligand_cores_receptor_exp, 2, divmax)*10, sparse = T)

      aggregated_score = ligand_conc_divmax * ligand_act_divmax * ligand_cores_receptor_exp_divmax

      lidardata@Lidarscore@aggregated_score = aggregated_score
    },
    pairwise = {
      ligand_act = lidardata@Diffusion@ligand_act
      ligand_dist = lidardata@Diffusion@ligand_dist
      dist_pos_bin_all = lidardata@Image@dist_bin_all
      dist_threshold = lidardata@Diffusion@dist_threshold
      ligandname = lidardata@DB@ligandname
      receptorname = lidardata@DB@receptorname
      lr_list = lidardata@DB@lr_list
      if (imp) {
        ligand_exp = lidardata@ImputationData@ligand_exp
      } else {
        ligand_exp = lidardata@Data@ligand_exp
      }
      if (imp) {
        receptor_exp = lidardata@ImputationData@receptor_exp
      } else {
        receptor_exp = lidardata@Data@receptor_exp
      }
      nligand = lidardata@DB@nligand
      ncell = lidardata@ncell
      cellname = lidardata@cellname

      ligand_diff = lidardata@Diffusion@ligand_diff

      dist_threshold = lidardata@Diffusion@dist_threshold
      dist_x = 0:(dist_threshold-1)

      #succ_ligand = lidardata@Diffusion@ligand_diff %>% filter(p.value < pval) %>% pull(Ligand)
      if (!is.null(ligand_list)) {
        succ_ligand <- ligand_list
      } else {
        succ_ligand <- lidardata@Diffusion@ligand_diff %>%
          filter(p.value < pval) %>%
          pull(Ligand)
      }
      if (!all(succ_ligand %in% lidardata@Diffusion@ligand_diff$Ligand)) stop("Requested ligands must have fitted diffusion results.")
      if (!length(succ_ligand)) stop("No ligands pass the score selection; adjust pval or supply ligand_list.")
      n_succ_ligand = length(succ_ligand)
      sender_score = Matrix(0,nrow = n_succ_ligand, ncol = ncell, dimnames = list(succ_ligand, cellname), sparse = T)
      receiver_score = Matrix(0,nrow = n_succ_ligand, ncol = ncell, dimnames = list(succ_ligand, cellname), sparse = T)

      pairwise_lidarscore = list()

      for (i in seq_along(succ_ligand)) {
        current_ligand <- succ_ligand[i]
        print(i)

        # Progress

        #if (nligand > 20 && i / nligand >= prog[n_prog]) {
        #  message(sprintf("Processing %d%%", as.integer(prog[n_prog] * 100)))
        #  n_prog <- n_prog + 1
        #}

        # Get ligand expression
        one_ligand_exp <- ligand_exp[current_ligand,] %>% as.vector()
        one_ligand_exp_lidar <- Matrix(divmax(one_ligand_exp)*10, sparse = T)

        # Filter receptors for the current_ligand
        receptor <- filter(lr_list, from %in% current_ligand) %>% dplyr::select(to) %>% unlist() %>% unique()

        # Get receptor expression
        if(length(receptor) > 1) {
          one_receptor_exp <- receptor_exp[receptor,] %>% as.array() %>% colSums()
        } else {
          one_receptor_exp <- receptor_exp[receptor,]
        }
        one_receptor_exp_lidar <- Matrix(t(divmax(one_receptor_exp)*10), sparse = T)

        # Extract ligand activity and dissociation
        one_ligand_act <- ligand_act[, colnames(ligand_act) == current_ligand] %>% as.vector()
        one_ligand_act_lidar <- Matrix(t(divmax(one_ligand_act)*10), sparse = T)

        # Diffusion
        ki = ligand_diff %>% filter(Ligand == current_ligand) %>% pull(ki)
        n = ligand_diff %>% filter(Ligand == current_ligand) %>% pull(n)
        one_ligand_diff <- 1 / (1 + ((dist_x / (dist_threshold - 1)) / ki)^n)
        diffuion_lidar = round(maxmin(one_ligand_diff), 3)
        diffuion_lidar[diffuion_lidar < 0.01] <- 0

        dist_pos_bin_all[dist_pos_bin_all>dist_threshold] = 0
        dist_pos_bin = Matrix::Matrix(dist_pos_bin_all, sparse = T)
        dist_pos_diff = dist_pos_bin

        #250826 change
        for (j in 1:dist_threshold){
          if (j == 1 && !include_self) {
            dist_pos_diff[dist_pos_diff == j] = 0
          } else {
            dist_pos_diff[dist_pos_diff == j] = diffuion_lidar[j]
          }
        }
        #250826 change over
        pairwise_lidarscore[[current_ligand]] = sweep(sweep(dist_pos_diff, 1, one_ligand_exp_lidar, FUN = "*"), 2, one_receptor_exp_lidar*one_ligand_act_lidar, FUN = "*")
        sender_score[current_ligand,] = rowSums(pairwise_lidarscore[[current_ligand]])
        receiver_score[current_ligand,] = colSums(pairwise_lidarscore[[current_ligand]])
      }
      lidardata@Lidarscore@pairwise_score = pairwise_lidarscore
      lidardata@Lidarscore@sender_score = sender_score
      lidardata@Lidarscore@receiver_score = receiver_score
    }
  )
  return(lidardata)
}

divmax = function(x){
  m <- max(x)
  if (!is.finite(m) || m == 0) return(x * 0)
  x/m
}
