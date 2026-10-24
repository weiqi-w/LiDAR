InferLigandAct <- function(lidardata, normalize = TRUE, parallel = FALSE) {
  cell_data = lidardata@cellname
  ligandname = lidardata@DB@ligandname
  n_list = lidardata@Coembedding@n_list
  nets = lidardata@Coembedding@nets
  lr_list = lidardata@DB@lr_list
  receptor_exp = lidardata@Data@receptor_exp
  receptorname = lidardata@DB@receptorname

  ligand_seed = GetLigandSeeds(lidardata)

  rwrs = list()
  for (i in seq_along(nets)) {
    net  = nets[[i]]
    if (!requireNamespace("dnet", quietly = TRUE)) stop("Install dnet to infer ligand activity.")
    rwr = dnet::dRWR(net, setSeeds = ligand_seed, parallel = parallel, verbose = FALSE, restart = 0.75)

    #Extract RWR results for target cells
    rownames(rwr) = V(net)$name
    rwr_cell = rwr[cell_data, , drop = FALSE]
    colnames(rwr_cell) = colnames(ligand_seed)
    rwrs[[i]] = rwr_cell
  }

  rwr_best_nn = SelectBestNN(rwrs, n_list, ligandname, cell_data, lr_list, receptor_exp, receptorname, normalize)
  lidardata@Diffusion@ligand_act = rwr_best_nn

  return(lidardata)
}

#GetLigandSeeds----
GetLigandSeeds <- function(lidardata) {

  # Extract weight_seed from the ltm matrix (only for ligands with receptors)
  nets = lidardata@Coembedding@nets
  weight_seed = lidardata@DB@ltm
  nligand = lidardata@DB@nligand
  ligandname = lidardata@DB@ligandname

  # Create a zero matrix for seeds (cells x ligands)
  zero_seed = matrix(0, nrow = length(V(nets[[1]])$name),
                     ncol = nligand,
                     dimnames = list(V(nets[[1]])$name, ligandname))

  # Fill in the weight_seed into the zero_seed matrix

  zero_seed[match(rownames(weight_seed), V(nets[[1]])$name), ] = weight_seed
  weight_seed = zero_seed

  # Select top 1% target genes as seed nodes (apply a threshold)
  weight_seed = apply(weight_seed, 2, function(x) {
    threshold <- quantile(x, 0.99)
    x[x < threshold] <- 0
    return(x)
  })

  weight_seed <- matrix(weight_seed, ncol = nligand, dimnames = list(V(nets[[1]])$name, ligandname))
  message(paste("Number of Candidate ligands:", nligand))
  return(weight_seed)
}



#Select best nn for each ligand----
SelectBestNN = function(rwrs, n_list, ligandname, cell_data, lr_list, receptor_exp, receptorname, normalize = TRUE){
  rwr_best_nn = matrix(0, nrow = nrow(rwrs[[1]]), ncol = ncol(rwrs[[1]]))
  r = matrix(0, nrow = length(n_list), ncol = length(ligandname))
  #best_n = list()
  for (i in 1:length(ligandname)){
    ligand = ligandname[i]
    for (j in 1:length(n_list)){
      split_rwr = split_rwr_by_receptor(rwrs,ligand, lr_list, receptor_exp, receptorname,j)
      if (sd(split_rwr$rwr) == 0 || sd(split_rwr$receptor) == 0) {
        r[j, i] = 0
      } else {
        r[j, i] = cor(split_rwr$rwr, split_rwr$receptor)
      }
    }
    best_nn = which.max(r[,i])
    #best_n[[i]] = best_nn
    rwr_best_nn[,i] = rwrs[[best_nn]][,ligand]
  }
  if (any(!is.finite(rwr_best_nn)) || any(rwr_best_nn < 0)) {
    stop("RWR activity must contain finite, nonnegative values.")
  }
  if (normalize == T){
    #rwr_best_nn = log((rwr_best_nn / colSums(rwr_best_nn))* 10000 + 1)
    totals <- colSums(rwr_best_nn)
    zero_activity <- totals == 0
    if (any(zero_activity)) {
      warning("No RWR activity reached cells for ligand(s): ",
              paste(ligandname[zero_activity], collapse = ", "),
              "; keeping their activity at zero.", call. = FALSE)
    }
    totals[zero_activity] <- 1
    rwr_best_nn = log1p(sweep(rwr_best_nn, 2, totals, "/") * 10000)
  }
  rownames(rwr_best_nn) = cell_data
  colnames(rwr_best_nn) = ligandname
  return(rwr_best_nn)
}

##Pearson test at different nn----
#need add no receptor = no cell
split_rwr_by_receptor = function(rwrs, ligand, lr_list, receptor_exp, receptorname, i){
  rwr = rwrs[[i]]
  rwr = rwr[,ligand]
  receptor = filter(lr_list, from %in% ligand ) %>% dplyr::select(to) %>% unlist() %>% unique()
  if(length(receptor) >1){
    one_receptor_exp = receptor_exp[match(receptor, receptorname),] %>% colSums()
  }else{
    one_receptor_exp = receptor_exp[match(receptor, receptorname),]
  }
  cell_with_receptor = names(one_receptor_exp)[one_receptor_exp>0]
  cell_no_receptor = names(one_receptor_exp)[one_receptor_exp==0]
  split_rwr = data.frame(cell = names(rwr), rwr = rwr, receptor = rep(0,length(rwr)))
  split_rwr$receptor = ifelse(split_rwr$cell %in% cell_with_receptor, 1, 0)
  return(split_rwr)
}

#Ligand_dis----
##Need add realdis = T----
