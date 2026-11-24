library(LiDAR)
local({
# Usage: Rscript --vanilla test_publication.R [path/to/Method_publication.R]
checks <- 0L
check <- function(ok) { stopifnot(isTRUE(ok)); checks <<- checks+1L }
bad <- function(expr) inherits(tryCatch({force(expr);NULL},error=function(e)e),"error")
set.seed(41)
cells <- paste0("c",seq_len(120)); ligands <- c("L1","L2"); receptors <- c("R1","R2")
pos <- data.frame(x=seq(0,595,5),y=0,row.names=cells)
lig <- matrix(rexp(240),2,dimnames=list(ligands,cells))
rec <- matrix(1,2,120,dimnames=list(receptors,cells))
z <- .EmptyMatrix()
make_object <- function() methods::new("LidarData",
  Data=methods::new("Data",data=rbind(lig,rec),ligand_exp=lig,receptor_exp=rec),
  ImputationData=methods::new("ImputationData",data=z,ligand_exp=lig,receptor_exp=rec),
  Image=methods::new("Image",pos=pos,dist=z,dist_bin_all=Matrix::Matrix(z,sparse=TRUE)),
  Coembedding=methods::new("Coembedding",cell_embedding=matrix(rnorm(240),120,dimnames=list(cells,NULL)),
    gene_embedding=matrix(1,2,2,dimnames=list(c("g1","g2"),NULL)),nets=list(),n_list=numeric()),
  Diffusion=methods::new("Diffusion",ligand_act=.EmptyMatrix(cells),ligand_dist=array(0,c(0,0,0)),
    ligand_diff=.EmptyFit(),ligand_conc=.EmptyMatrix(cells),dist_threshold=numeric(),bin=numeric(),failed_ligand=character()),
  Lidarscore=.ScoresEmpty(cells),DB=methods::new("DB",ltm=matrix(1,2,2,dimnames=list(c("g1","g2"),ligands)),
    lr_list=data.frame(from=ligands,to=receptors,database="curated"),ligandname=ligands,receptorname=receptors,nligand=2L,nreceptor=2L),
  cellname=cells,genename=c(ligands,receptors),ncell=120,ngene=4)

# Coordinates: precedence, row ordering, aliases and nonnumeric rejection.
p <- data.frame(barcode=rev(cells),X=rev(pos$x),Y=0,unrelated=999)
coord <- .Coords(p,cells)
check(identical(rownames(coord),cells) && identical(coord$x,pos$x) && ncol(coord)==2L)
check(bad(.Coords(transform(p,X=as.character(X)),cells)))
check(bad(.Coords(p,cells,c("X","missing"))))

# Eleven bins, self contribution, boundary behavior, automatic minimum distance.
x <- BinLigandExpr(make_object(),bin=50,n_bins=10,imp=FALSE)
check(dim(x@Diffusion@ligand_dist)[3]==11)
check(all(as.numeric(x@Image@dist_bin_all[1,c(1,10,11,101,111)]) == c(1,1,2,11,0)))
expected <- sum(lig["L1",pos$x<50])
check(abs(x@Diffusion@ligand_dist["L1","c1",1]-expected)<1e-12)
auto <- BinLigandExpr(make_object(),imp=FALSE)
check(auto@Diffusion@bin==5)
dup <- make_object(); dup@Image@pos$x <- 0
check(bad(BinLigandExpr(dup,imp=FALSE)))
check(!bad(BinLigandExpr(dup,bin=50,imp=FALSE)))

# Stable continuous decay and analytical distances.
for (pars in list(c(.3,2),c(.01,.01),c(1,10))) {
  d <- seq(0,500,length.out=101); w <- .DecayWeights(pars,d,500)
  old <- 1/(1+((d/500)/pars[1])^pars[2]); old <- (old-min(old))/(max(old)-min(old))
  check(max(abs(w-old))<1e-10 && w[1]==1 && w[101]==0 && all(diff(w)<=0))
  rr <- CalcRanges(pars[1],pars[2],500)
  check(rr$r50<=rr$r01 && rr$r01<=500 && rr$r50>=0)
  if (rr$r50>0) check(abs(.DecayWeights(pars,rr$r50,500)-.5)<1e-9)
  if (rr$r01>0) check(abs(.DecayWeights(pars,rr$r01,500)-.01)<1e-9)
}
check(identical(ScaleMax(c(0,0)),c(0,0)))
check(bad(ScaleMax(c(0,NA))))
check(bad(.DecayWeights(c(0,1),0:10,10)))

# Gram-based Pearson correlation equals direct correlation.
X <- .LigandBins(x@Diffusion@ligand_dist,"L1",cells)
truew <- .DecayWeights(c(.4,2),seq(0,500,50),500)
y <- as.numeric(CalcExposure(X,truew)) + runif(120,0,.05)
prepared <- .PrepCoupling(X,y)
for (i in seq_len(20)) {
  w <- .DecayWeights(c(runif(1,.01,1),runif(1,.01,10)),seq(0,500,50),500)
  fast <- sum(w*crossprod(prepared$Xc,prepared$yc))/sqrt(as.numeric(crossprod(w,prepared$A%*%w))*prepared$ss)
  check(abs(fast-cor(y,as.numeric(X%*%w)))<1e-10)
}
check(bad(.PrepCoupling(X,rep(0,120))))

# Candidate zero response excludes only that network; nonfinite excludes ligand.
r1 <- matrix(0,120,2,dimnames=list(cells,ligands)); r2 <- r1
r1[,"L2"] <- 1; r2[,"L1"] <- seq_len(120); r2[,"L2"] <- NA
choice <- SelectNetwork(list(r1,r2),c(5,10),ligands,cells,x@DB@lr_list,rec)
check(identical(colnames(choice),"L1"))
check(attr(choice,"network_selection")$n==10)

# Imputation preserves nonzeros and replaces only zeros, including sparse input.
nbrs <- list(nn.idx=matrix(c(2,3,1),3,1),nn.weights=matrix(1,1,3),cells=c("a","b","c"))
ex <- matrix(c(0,2,4),1,dimnames=list("g",nbrs$cells))
check(identical(as.numeric(ImputeZeros(nbrs,ex)),c(2,2,4)))
check(identical(as.numeric(ImputeZeros(nbrs,Matrix::Matrix(ex,sparse=TRUE))),c(2,2,4)))

# Full downstream workflow: optional permutations, exposure, ranges and scores.
x@Diffusion@ligand_act <- matrix(y,120,1,dimnames=list(cells,"L1"))
f <- FitPropagation(x,imp=FALSE)
check(nrow(f@Diffusion@ligand_diff)==1L)
check(!"p.empirical" %in% names(f@Diffusion@ligand_diff))
check(all(c("r50","r01") %in% names(f@Diffusion@ligand_diff)))
check(f@Diffusion@ligand_diff$coupling>.99)
set.seed(555); rng <- .Random.seed
fp <- FitPropagation(x,imp=FALSE,permutation_p=TRUE,n_perm=9,perm_seed=123)
check(identical(rng,.Random.seed))
fp2 <- FitPropagation(x,imp=FALSE,permutation_p=TRUE,n_perm=9,perm_seed=123)
pr <- attr(fp@Diffusion@ligand_diff,"permutation_results")[["L1"]]
pr2 <- attr(fp2@Diffusion@ligand_diff,"permutation_results")[["L1"]]
check(identical(pr$null_cor,pr2$null_cor))
check(pr$n_perm_valid==9 && pr$p.empirical==(1+sum(pr$null_cor>=pr$observed))/10)
f <- ComputeScores(f,"aggregate",ligand_list="L1",imp=FALSE)
f <- ComputeScores(f,"pairwise",ligand_list="L1",imp=FALSE)
check(identical(dim(f@Lidarscore@aggregated_score),c(120L,1L)))
pair <- f@Lidarscore@pairwise_score$L1
check(max(abs(as.numeric(Matrix::rowSums(pair))-as.numeric(f@Lidarscore@sender_score["L1",])))<1e-10)
check(max(abs(as.numeric(Matrix::colSums(pair))-as.numeric(f@Lidarscore@receiver_score["L1",])))<1e-10)
check(all(diag(as.matrix(pair))>0))
check(pair[1,101]==0 && pair[1,111]==0)
check(max(pair)<=1000+1e-10 && max(f@Lidarscore@aggregated_score)<=1000+1e-10)
logged <- ComputeScores(f, "aggregate", imp=FALSE, Log=TRUE)
logged <- ComputeScores(logged, "pairwise", imp=FALSE, Log=TRUE)
for (slot_name in c("aggregated_score", "sender_score", "receiver_score")) {
  raw <- as.matrix(methods::slot(f@Lidarscore, slot_name))
  observed <- as.matrix(methods::slot(logged@Lidarscore, slot_name))
  check(isTRUE(all.equal(observed, log10(1 + raw))))
}
check(isTRUE(all.equal(as.matrix(logged@Lidarscore@pairwise_score$L1), log10(1 + as.matrix(pair)))))
for (boundary in c("p.value", "coupling")) {
  excluded <- f
  excluded@Diffusion@ligand_diff[[boundary]] <- if (boundary == "p.value") 0.001 else 0
  excluded <- ComputeScores(excluded, ligand_list="L1", imp=FALSE)
  check(ncol(excluded@Lidarscore@aggregated_score) == 0L)
  restored <- ComputeScores(excluded, ligand_list="L1", imp=FALSE, pval=NULL, coupling=NULL)
  check(ncol(restored@Lidarscore@aggregated_score) == 1L)
}
# Exact pairwise formula, and receiver sums consistent with unscaled exposure.
cfg <- attr(f@Diffusion@ligand_diff,"fit_config"); fit <- f@Diffusion@ligand_diff[1,]
weights <- .DecayWeights(c(fit$ki,fit$n),cfg$distance,cfg$B)
check(abs(pair[1,12] - ScaleMax(lig["L1",])[1]*10*weights[2]*10*ScaleMax(y)[12]*10)<1e-9)
check(bad(ComputeScores(f,"aggregate",imp=TRUE)))

# Failure paths: valid empty results, constant response, NA and missing receptor.
failed <- x; failed@Diffusion@ligand_act[,] <- 1
ff <- FitPropagation(failed,imp=FALSE,permutation_p=TRUE,n_perm=2)
check(nrow(ff@Diffusion@ligand_diff)==0L && ncol(ff@Diffusion@ligand_conc)==0L)
check(methods::validObject(ff))
failed <- x; failed@Diffusion@ligand_act[1,1] <- NA
check(nrow(FitPropagation(failed,imp=FALSE)@Diffusion@ligand_diff)==0L)
failed <- x; failed@DB@lr_list <- failed@DB@lr_list[FALSE,,drop=FALSE]
check(nrow(FitPropagation(failed,imp=FALSE)@Diffusion@ligand_diff)==0L)
zero <- f; zero@Data@ligand_exp["L1",] <- 0
zero <- ComputeScores(zero,"pairwise",ligand_list="L1",imp=FALSE)
check(all(zero@Lidarscore@pairwise_score$L1==0))
empty <- ComputeScores(f,"pairwise",ligand_list=character(),imp=FALSE)
check(length(empty@Lidarscore@pairwise_score)==0L && nrow(empty@Lidarscore@sender_score)==0L)
tmp <- tempfile(fileext=".rds"); saveRDS(fp,tmp); restored <- readRDS(tmp); unlink(tmp)
check(identical(attr(restored@Diffusion@ligand_diff,"permutation_results"),attr(fp@Diffusion@ligand_diff,"permutation_results")))
# Explicit failed-draw handling: never substitute another permutation or a zero.
mock_env <- new.env(parent=asNamespace("LiDAR"))
mock_env$.FitDecayFast <- function(...) stop("test failure")
permuter <- .PermuteCoupling
environment(permuter) <- mock_env
set.seed(99)
failure_test <- suppressWarnings(permuter(prepared,.5,seq(0,500,50),500,3L,300))
check(is.na(failure_test$p.empirical) && failure_test$n_perm_failed==3L &&
      all(is.na(failure_test$null_cor)))

if (requireNamespace("RANN",quietly=TRUE) && requireNamespace("igraph",quietly=TRUE)) {
  # Duplicate embeddings cannot leave self in the imputation neighborhood.
  embedding <- matrix(c(0,0,1,2,0,0,0,0),4,2,dimnames=list(letters[1:4],NULL))
  nb <- FindNeighbors(embedding,1)
  check(identical(dim(nb$nn.idx),c(4L,1L)) && all(nb$nn.idx[,1]!=seq_len(4)))
  check(all(abs(colSums(nb$nn.weights)-1)<1e-12))
  check(bad(FindNeighbors(embedding,4)))
  netobject <- BuildNetworks(make_object(),c(1,3))
  check(all(vapply(netobject@Coembedding@nets,igraph::is_simple,logical(1))))
  check(all(vapply(netobject@Coembedding@nets,function(g) !igraph::is_directed(g),logical(1))))
  seeds <- SelectSeeds(netobject)
  check(all(seeds[cells,]==0) && all(seeds[c("g1","g2"),]==1))
  # Two positive target weights: only the top survives, regardless of 120 cell zeros.
  netobject@DB@ltm[,"L1"] <- c(1,2)
  seeds <- SelectSeeds(netobject)
  check(seeds["g1","L1"]==0 && seeds["g2","L1"]==2)
  im <- ImputeExpr(netobject,nn=1)
  check(identical(as.matrix(im@ImputationData@ligand_exp),lig))
} else message("SKIP network/seed tests: RANN and/or igraph not installed")
# Installed-package objects retain a valid schema.
check(methods::validObject(fp))
cat("PASS:",checks,"checks\n")
# Run from the directory containing Method_publication.R, LiDAR_plots.R and
# test_publication.R: Rscript --vanilla test_plots.R
plot_checks<-0L
pc<-function(ok){stopifnot(isTRUE(ok));plot_checks<<-plot_checks+1L}
pb<-function(expr)inherits(tryCatch({force(expr);NULL},error=function(e)e),"error")
pc(length(lidar_colors)==120L)
pc(length(select_distinct_colors(0))==0L)
pc(identical(select_distinct_colors(1),lidar_colors[1]))
pc(length(unique(select_distinct_colors(10)))==10L)
pc(pb(select_distinct_colors(-1))&&pb(select_distinct_colors(1.5))&&pb(select_distinct_colors(121)))
pc(identical(order_colors_by_similarity("red"),"red")&&length(order_colors_by_similarity(character()))==0)
pc(pb(order_colors_by_similarity("notacolor")))
pc(inherits(theme_Lidar1(),"theme")&&length(theme_Lidar1(ratio=1))==2L)
pc(inherits(theme_Lidar2(border=TRUE)$panel.border,"element_rect"))
pc(inherits(theme_LiDAR(show_xy=TRUE)$axis.line,"element_line"))

selected_cells<-rev(cells[1:20])
coords<-GetCoords(f,selected_cells)
pc(identical(coords$cell,selected_cells))
pc(pb(GetCoords(f,c("missing"))))
pc(identical(rownames(GetResponse(f,"L1",selected_cells)),selected_cells))
pc(pb(GetResponse(f,"missing")))
pc(identical(dim(GetAvailability(f,"L1",cells=selected_cells,imp=FALSE)),c(20L,1L)))
pc(all(GetAvailability(f,"L1",cells=selected_cells,imp=FALSE)==1))
# Raw versus imputed expression must differ when their inputs differ.
different<-f;different@ImputationData@ligand_exp<-different@Data@ligand_exp*2
features<-GetFeatures(different,features=c("Ligand_exp","Ligand_exp_imp"))
pc(all(features$LigandExprImp==2*features$LigandExpr))
pc(pb(GetFeatures(f,features="unknown")))
pc(identical(GetDecay(f,"L1")$Ligand,"L1"))
pc(pb(GetDecay(f,pval=.05,p_column="p.empirical")))
pc(identical(SelectLigands(fp,pval=.2,p_column="p.empirical"),"L1"))
pc(pb(SelectLigands(f,pct=.1)))
legacy <- f
legacy_fit <- legacy@Diffusion@ligand_diff
names(legacy_fit)[names(legacy_fit) == "coupling"] <- "r"
legacy@Diffusion@ligand_diff <- legacy_fit
pc(identical(GetDecay(legacy), GetDecay(f)))
pc(identical(GetDecay(legacy, r=c(0,1)), GetDecay(f, coupling=c(0,1))))
pc(identical(SelectLigands(legacy, r=0), SelectLigands(f, coupling=0)))
pc(nrow(GetDecay(f,coupling=c(-1,0)))==0L)

values<-setNames(as.numeric(GetResponse(f,"L1")[,"L1"]),cells)
p<-PlotFeature(values[rev(seq_along(values))],GetCoords(f),show_border=TRUE,zscore=FALSE)
pc(identical(p$data$feature,as.numeric(values)))
pc(inherits(p$theme$panel.border,"element_rect"))
pc(identical(formals(ExportPlots)$ratio, 0.8))
pc(isTRUE(all.equal(PlotFeature(values, GetCoords(f), ratio=1.2)$coordinates$ratio, 1.2)))
pc(pb(ExportPlots(f, ratio=0)))
pz<-PlotFeature(values,GetCoords(f))
pc(abs(mean(pz$data$feature))<1e-12&&abs(sd(pz$data$feature)-1)<1e-12)
palette_at <- function(feature) {
  panel <- PlotFeature(values, GetCoords(f), feature=feature)
  panel$scales$get_scales("colour")$palette(c(0, 0.5, 1))
}
pc(identical(palette_at("Score"), palette_at("ReceiverScore")))
pc(identical(palette_at("SenderScore"), palette_at("Exposure")))
pc(!identical(palette_at("LigandExpr"), palette_at("Response")))
pc(identical(formals(BuildNetworks)$n, quote(c(500, 700))))
pc(pb(PlotFeature(unname(values),GetCoords(f))))
constant<-suppressWarnings(PlotFeature(setNames(rep(3,length(cells)),cells),GetCoords(f),scale="zscore"))
pc(all(constant$data$feature==0))
bad_values<-values;bad_values[1]<-NA
pc(pb(PlotFeature(bad_values,GetCoords(f))))
flipped<-PlotFeature(values,GetCoords(f),reverse_y=TRUE)
pc(!is.null(flipped$scales$get_scales("y")))

pd<-PlotDecay(fp,"L1")
pc(nrow(pd$data)==501L&&pd$data$weight[1]==1&&tail(pd$data$weight,1)==0)
pc(max(pd$data$distance)==500)
pc(grepl("empirical P",pd$labels$subtitle,fixed=TRUE))
pc(!grepl("empirical P",PlotDecay(f,"L1")$labels$subtitle,fixed=TRUE))
pc(grepl("input coordinate units",pd$labels$x,fixed=TRUE))
coupling<-PlotCoupling(fp,"L1")
pc(grepl("empirical P",coupling$labels$subtitle,fixed=TRUE))
subset_plot<-PlotCoupling(fp,"L1",cells=selected_cells)
pc(!grepl("empirical P",subset_plot$labels$subtitle,fixed=TRUE))
pc(grepl("Selected receivers",subset_plot$labels$subtitle,fixed=TRUE))
pc(pb(PlotCoupling(fp,"L1",cells=cells[1:2])))
stale<-fp;stale@Diffusion@ligand_act[,1]<-rev(stale@Diffusion@ligand_act[,1])
pc(pb(PlotCoupling(stale,"L1")))

folder<-tempfile("lidar_plot_tests_");dir.create(folder)
# Rendering tests force ggplot/patchwork to evaluate their layers and scales.
ggplot2::ggsave(file.path(folder,"feature.png"),p,width=6,height=4,dpi=90)
ggplot2::ggsave(file.path(folder,"decay.png"),pd,width=7,height=4,dpi=90)
ggplot2::ggsave(file.path(folder,"coupling.png"),coupling,width=7,height=4,dpi=90)
pc(all(file.info(file.path(folder,c("feature.png","decay.png","coupling.png")))$size>1000))
exports<-ExportPlots(f,features=c("LigandExpr","LigandExprImp"),out_dir=file.path(folder,"shared"),ratio=1.2,zscore=FALSE,format="png",dpi=60)
pc(length(exports$files)==2L)
pc(isTRUE(all.equal(exports$color_limits$LigandExpr,range(lig))))
pc(all(vapply(exports$files,function(a)file.exists(a$path),logical(1))))
independent<-ExportPlots(f,features="LigandExpr",out_dir=file.path(folder,"independent"),shared_scales=FALSE,format="png",dpi=60)
pc(length(independent$color_limits)==0L)
only_decay<-ExportPlots(f,ligand_list="L1",features="Decay",out_dir=file.path(folder,"decayonly"),format="png",dpi=60)
pc(identical(only_decay$files$L1$features,"Decay"))
only_cor<-ExportPlots(f,ligand_list="L1",features="Coupling",out_dir=file.path(folder,"coronly"),format="png",dpi=60)
pc(identical(only_cor$files$L1$features,"Coupling"))
repeat_export<-suppressWarnings(ExportPlots(f,ligand_list="L1",features="Decay",out_dir=file.path(folder,"decayonly"),format="png",dpi=60))
pc(length(repeat_export$files)==0L&&any(grepl("File exists",repeat_export$log$message)))
missing<-suppressWarnings(ExportPlots(f,ligand_list="L2",features="Decay",out_dir=file.path(folder,"missing"),format="png",dpi=60))
pc(length(missing$files)==0L&&any(missing$log$status=="skipped"))
zexport<-ExportPlots(f,features="LigandExpr",scale="zscore",out_dir=file.path(folder,"zscore"),format="png",dpi=60)
pc(length(zexport$files)==2L&&length(zexport$color_limits$LigandExpr)==2L)
pc(pb(ExportPlots(f,ligand_list="missing",out_dir=folder)))
pc(pb(ExportPlots(f,features="unknown",out_dir=folder)))
cat("PLOT PASS:",plot_checks,"checks\n")
cat("Rendered images:",folder,"\n")
preview_dir<-Sys.getenv("LIDAR_PREVIEW_DIR",unset="")
if(nzchar(preview_dir)) {
  dir.create(preview_dir,recursive=TRUE,showWarnings=FALSE)
  file.copy(file.path(folder,c("feature.png","decay.png","coupling.png")),preview_dir,overwrite=TRUE)
}
nodes<-c('n3','n1','n4','n2'); ligands<-c('L2','L1')
g<-igraph::make_ring(4);igraph::V(g)$name<-nodes
s<-matrix(c(1,2,3,4,4,3,2,1),4,dimnames=list(nodes,ligands))
s<-s[c(4,2,1,3),,drop=FALSE]
r<-dnet::dRWR(g,setSeeds=s,restart=1,parallel=FALSE,verbose=FALSE)
a<-.AlignRWR(r,nodes,ligands)
expected<-s[nodes,,drop=FALSE]/10
stopifnot(isTRUE(all.equal(a,expected)))
stopifnot(identical(.AlignRWR(a[c(3,1,4,2),2:1],nodes,ligands),a))
stopifnot(identical(.AlignRWR(unname(a[,1,drop=FALSE]),nodes,ligands[1]),a[,1,drop=FALSE]))
bad<-function(x) inherits(tryCatch(.AlignRWR(x,nodes,ligands),error=function(e)e),'error')
stopifnot(bad(a[-1,,drop=FALSE]))
x<-a;rownames(x)[1]<-'wrong';stopifnot(bad(x))
x<-a;colnames(x)<-c('L1','L1');stopifnot(bad(x))
r1<-.AlignRWR(dnet::dRWR(g,setSeeds=s,restart=.75,parallel=FALSE,verbose=FALSE),nodes,ligands)
g2<-igraph::permute(g,c(3,1,4,2));nodes2<-igraph::V(g2)$name
r2<-.AlignRWR(dnet::dRWR(g2,setSeeds=s[rev(rownames(s)),2:1,drop=FALSE],restart=.75,parallel=FALSE,verbose=FALSE),nodes2,rev(ligands))
stopifnot(isTRUE(all.equal(r1,r2[nodes,ligands],tolerance=1e-10)))
cat('PASS: dnet ordering, single ligand, named reordering, dimension/name rejection, permuted graph and seed invariance\n')
}, envir = new.env(parent = asNamespace("LiDAR")))
