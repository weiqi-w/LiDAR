methods::setClassUnion("LiDARMatrix", c("matrix", "Matrix"))
.Need <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) 
        stop("Required package is not installed: ", pkg, call. = FALSE)
}
methods::setClass("Data", slots = c(data = "LiDARMatrix", ligand_exp = "LiDARMatrix", 
    receptor_exp = "LiDARMatrix"))
methods::setClass("ImputationData", slots = c(data = "LiDARMatrix", 
    ligand_exp = "LiDARMatrix", receptor_exp = "LiDARMatrix"))
methods::setClass("Image", slots = c(pos = "data.frame", dist = "LiDARMatrix", 
    dist_bin_all = "Matrix"))
methods::setClass("Coembedding", slots = c(cell_embedding = "LiDARMatrix", 
    gene_embedding = "LiDARMatrix", nets = "list", n_list = "numeric"))
methods::setClass("Diffusion", slots = c(ligand_act = "LiDARMatrix", 
    ligand_dist = "array", ligand_diff = "data.frame", ligand_conc = "LiDARMatrix", 
    dist_threshold = "numeric", bin = "numeric", failed_ligand = "character"))
methods::setClass("Lidarscore", slots = c(aggregated_score = "LiDARMatrix", 
    pairwise_score = "list", sender_score = "LiDARMatrix", receiver_score = "LiDARMatrix"))
methods::setClass("DB", slots = c(ltm = "LiDARMatrix", lr_list = "data.frame", 
    ligandname = "character", receptorname = "character", nligand = "integer", 
    nreceptor = "integer"))
methods::setClass("LidarData", slots = c(Data = "Data", ImputationData = "ImputationData", 
    Image = "Image", Coembedding = "Coembedding", Diffusion = "Diffusion", 
    Lidarscore = "Lidarscore", DB = "DB", cellname = "character", 
    genename = "character", ncell = "numeric", ngene = "numeric"))
