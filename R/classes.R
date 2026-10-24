methods::setClassUnion("LidarMatrix", c("matrix", "Matrix"))

setClass("Data",
         slots = c(
           data = "LidarMatrix",
           ligand_exp = "LidarMatrix",
           receptor_exp = "LidarMatrix"
         ))

setClass("ImputationData",
         slots = c(
           data = "LidarMatrix",
           ligand_exp = "LidarMatrix",
           receptor_exp = "LidarMatrix"
         ))

setClass("Image",
         slots = c(
           pos = "data.frame",
           dist = "LidarMatrix",
           dist_bin_all = "Matrix"
         ))


setClass("Coembedding",
         slots = c(
           cell_embedding = "LidarMatrix",
           gene_embedding = "LidarMatrix",
           nets = "list",
           n_list = "numeric"
         ))

setClass("Diffusion",
         slots = c(
           ligand_act = "LidarMatrix",
           ligand_dist = "array",
           ligand_diff = "data.frame",
           ligand_conc = "LidarMatrix",
           dist_threshold = "numeric",
           bin = "numeric",
           include_self = "logical",
           failed_ligand = "character"
         ))

setClass("Lidarscore",
         slots = c(
           aggregated_score = "LidarMatrix",
           pairwise_score = "list",
           sender_score = "LidarMatrix",
           receiver_score = "LidarMatrix"
         ))

setClass("DB",
         slots = c(
           ltm = "LidarMatrix",
           lr_list = "data.frame",
           ligandname = "character",
           receptorname = "character",
           nligand = "integer",
           nreceptor = "integer"
         ))

setClass("LidarData",
         slots = c(
           Data = "Data",
           ImputationData = "ImputationData",
           Image = "Image",
           Coembedding = "Coembedding",
           Diffusion = "Diffusion",
           Lidarscore = "Lidarscore",
           DB = "DB",
           cellname = "character",
           genename = "character",
           ncell = "numeric",
           ngene = "numeric"
         ))
