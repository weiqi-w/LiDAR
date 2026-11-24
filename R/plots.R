# Combined LiDAR themes, result accessors and plots.
# Companion to the revised Method_publication.R. No files are sourced implicitly.
# Default: Z-scores, same-feature shared scales, continuous decay,
# unchanged Y direction. No expression-percentage ligand filtering.

.PlotNeed <- function(pkg) {
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Required package is not installed: ",pkg,call.=FALSE)
}
.PlotFlag <- function(x,name) {
  if (!is.logical(x)||length(x)!=1L||is.na(x)) stop(name," must be TRUE or FALSE")
}
.PlotNames <- function(x,name) {
  if (is.null(x)||anyNA(x)||any(!nzchar(x))||anyDuplicated(x)) stop(name," must have unique nonempty names")
}
.PlotPick <- function(request,available,label) {
  .PlotNames(available,label)
  if (is.null(request)) return(available)
  request <- as.character(request)
  .PlotNames(request,paste("Requested",label))
  missing <- setdiff(request,available)
  if (length(missing)) stop("Unknown ",label,": ",paste(missing,collapse=", "))
  request
}
.PlotNumber <- function(x,name,positive=FALSE) {
  if (!is.numeric(x)||length(x)!=1L||!is.finite(x)||
      (positive && x<=0)||(!positive && x<0)) stop("Invalid ",name)
}
.PlotFinite <- function(x,label) {
  if (any(!is.finite(x))) stop(label," contains NA/Inf")
}

theme_Lidar1 <- function(x_text_angle=FALSE,ratio=NULL,axis_text_size=14,
                         axis_title_size=14,legend_text_size=14,
                         legend_position="right",margin_size=35) {
  .PlotNeed("ggplot2"); .PlotFlag(x_text_angle,"x_text_angle")
  for (v in c("axis_text_size","axis_title_size","legend_text_size","margin_size")) .PlotNumber(get(v),v)
  th <- ggplot2::theme_test()+ggplot2::theme(
    panel.border=ggplot2::element_rect(colour="black",fill=NA,linewidth=.6),
    axis.ticks=ggplot2::element_line(colour="black",linewidth=.5),
    axis.text=ggplot2::element_text(size=axis_text_size,colour="black"),
    axis.text.x=ggplot2::element_text(angle=if(x_text_angle)45 else 0,hjust=if(x_text_angle)1 else .5),
    axis.title=ggplot2::element_text(size=axis_title_size,colour="black"),
    legend.title=ggplot2::element_text(size=legend_text_size),
    legend.text=ggplot2::element_text(size=legend_text_size),
    legend.position=legend_position,
    plot.margin=ggplot2::margin(margin_size,margin_size,margin_size,margin_size))
  if (is.null(ratio)) th else { .PlotNumber(ratio,"ratio",TRUE); list(ggplot2::coord_fixed(ratio),th) }
}
theme_Lidar2 <- function(border=FALSE,ratio=NULL,legend_text_size=14,
                         legend_position="right",margin_size=35) {
  .PlotNeed("ggplot2"); .PlotFlag(border,"border")
  .PlotNumber(legend_text_size,"legend_text_size"); .PlotNumber(margin_size,"margin_size")
  th <- ggplot2::theme_test()+ggplot2::theme(
    panel.border=if(border)ggplot2::element_rect(colour="black",fill=NA,linewidth=.6) else ggplot2::element_blank(),
    axis.ticks=ggplot2::element_blank(),axis.text=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),
    legend.title=ggplot2::element_text(size=legend_text_size),legend.text=ggplot2::element_text(size=legend_text_size),
    legend.position=legend_position,plot.margin=ggplot2::margin(margin_size,margin_size,margin_size,margin_size))
  if (is.null(ratio)) th else { .PlotNumber(ratio,"ratio",TRUE); list(ggplot2::coord_fixed(ratio),th) }
}
theme_LiDAR <- function(show_border=FALSE,show_xy_labels=FALSE,show_ticks=FALSE,
                        show_tick_labels=FALSE,show_title=TRUE,show_xy=TRUE,show_legend=TRUE,
                        text_size=14,xlab_size=NULL,xlab_angle=NULL,ylab_size=NULL,ylab_angle=NULL,
                        xtick_size=NULL,xtick_angle=NULL,ytick_size=NULL,ytick_angle=NULL,
                        title_size=NULL,title_angle=NULL,legend_title_size=NULL,legend_text_size=NULL) {
  .PlotNeed("ggplot2")
  for(v in c("show_border","show_xy_labels","show_ticks","show_tick_labels","show_title","show_xy","show_legend")) .PlotFlag(get(v),v)
  if(is.null(text_size))text_size<-14
  .PlotNumber(text_size,"text_size",TRUE)
  for(v in c("xlab_size","ylab_size","xtick_size","ytick_size","title_size","legend_title_size","legend_text_size"))
    if(!is.null(get(v))) .PlotNumber(get(v),v)
  ggplot2::theme_test(base_size=text_size)+ggplot2::theme(
    panel.border=if(show_border)ggplot2::element_rect(fill=NA,colour="black") else ggplot2::element_blank(),
    axis.line=if(show_xy)ggplot2::element_line(colour="black") else ggplot2::element_blank(),
    axis.ticks=if(show_ticks)ggplot2::element_line() else ggplot2::element_blank(),
    axis.title.x=if(show_xy_labels)ggplot2::element_text(size=xlab_size,angle=xlab_angle) else ggplot2::element_blank(),
    axis.title.y=if(show_xy_labels)ggplot2::element_text(size=ylab_size,angle=ylab_angle) else ggplot2::element_blank(),
    axis.text.x=if(show_tick_labels)ggplot2::element_text(size=xtick_size,angle=xtick_angle) else ggplot2::element_blank(),
    axis.text.y=if(show_tick_labels)ggplot2::element_text(size=ytick_size,angle=ytick_angle) else ggplot2::element_blank(),
    plot.title=if(show_title)ggplot2::element_text(size=title_size,angle=title_angle,hjust=.5) else ggplot2::element_blank(),
    legend.position=if(show_legend)"right" else "none",
    legend.title=ggplot2::element_text(size=legend_title_size),legend.text=ggplot2::element_text(size=legend_text_size))
}

select_distinct_colors <- function(n,palette=lidar_colors) {
  if(!is.numeric(n)||length(n)!=1L||!is.finite(n)||n<0||n!=floor(n)) stop("n must be a nonnegative integer")
  if(n>length(palette))stop("n exceeds palette size")
  if(!length(palette)||n==0L)return(character())
  lab <- grDevices::convertColor(t(grDevices::col2rgb(palette)/255),from="sRGB",to="Lab")
  if(n==1L)return(palette[1L])
  selected <- integer(n); selected[1L]<-1L
  for(i in seq.int(2L,n)) {
    previous <- selected[seq_len(i-1L)]
    remaining <- setdiff(seq_len(nrow(lab)),previous)
    distances <- vapply(remaining,function(j) {
      delta <- sweep(lab[previous,,drop=FALSE],2L,lab[j,],"-")
      min(sqrt(rowSums(delta^2)))
    },numeric(1))
    selected[i]<-remaining[which.max(distances)]
  }
  palette[selected]
}
order_colors_by_similarity <- function(colors) {
  if(!length(colors))return(colors)
  lab <- grDevices::convertColor(t(grDevices::col2rgb(colors)/255),from="sRGB",to="Lab")
  n<-length(colors); if(n==1L)return(colors)
  ordering<-integer(n);ordering[1L]<-1L;remaining<-seq.int(2L,n)
  for(i in seq.int(2L,n)) {
    distance<-rowSums(sweep(lab[remaining,,drop=FALSE],2L,lab[ordering[i-1L],],"-")^2)
    ordering[i]<-remaining[which.min(distance)];remaining<-setdiff(remaining,ordering[i])
  }
  colors[ordering]
}

GetCoords <- function(lidardata,cells=NULL) {
  pos<-lidardata@Image@pos
  if(!all(c("x","y")%in%names(pos)))stop("Expected normalized x/y coordinates from CreateLiDAR")
  ids<-.PlotPick(cells,rownames(pos),"coordinate cells")
  pos<-pos[ids,c("x","y"),drop=FALSE]
  if(!all(vapply(pos,is.numeric,logical(1))))stop("Coordinates must be numeric")
  .PlotFinite(as.matrix(pos),"Coordinates")
  data.frame(cell=ids,x=pos$x,y=pos$y,stringsAsFactors=FALSE,row.names=ids)
}
.PlotSubset <- function(mat,features=NULL,cells=NULL,label="matrix") {
  if(length(dim(mat))!=2L)stop(label," must be a matrix")
  ids<-.PlotPick(cells,rownames(mat),paste(label,"cells"))
  cols<-.PlotPick(features,colnames(mat),paste(label,"features"))
  mat[ids,cols,drop=FALSE]
}
GetLigandExpr <- function(lidardata,ligand_list=NULL,cells=NULL,imp=TRUE) {
  .PlotNeed("Matrix");.PlotFlag(imp,"imp")
  mat<-if(imp)lidardata@ImputationData@ligand_exp else lidardata@Data@ligand_exp
  .PlotSubset(t(mat),ligand_list,cells,"Ligand expression")
}
GetReceptorExpr <- function(lidardata,receptor_list=NULL,cells=NULL,imp=TRUE) {
  .PlotNeed("Matrix");.PlotFlag(imp,"imp")
  mat<-if(imp)lidardata@ImputationData@receptor_exp else lidardata@Data@receptor_exp
  .PlotSubset(t(mat),receptor_list,cells,"Receptor expression")
}
GetAvailability <- function(lidardata,ligand_list=NULL,receptor_list=NULL,cells=NULL,imp=TRUE) {
  .PlotNeed("Matrix")
  ligands<-.PlotPick(ligand_list,lidardata@DB@ligandname,"database ligands")
  expr<-GetReceptorExpr(lidardata,receptor_list,cells,imp)
  net<-lidardata@DB@lr_list
  ij<-unique(data.frame(i=match(net$to,colnames(expr)),j=match(net$from,ligands)))
  ij<-ij[stats::complete.cases(ij),,drop=FALSE]
  map<-Matrix::sparseMatrix(i=ij$i,j=ij$j,x=rep(1,nrow(ij)),
    dims=c(ncol(expr),length(ligands)),dimnames=list(colnames(expr),ligands))
  out<-expr%*%map;rownames(out)<-rownames(expr);out
}
GetResponse <- function(lidardata,ligand_list=NULL,cells=NULL)
  .PlotSubset(lidardata@Diffusion@ligand_act,ligand_list,cells,"Response")
GetExposure <- function(lidardata,ligand_list=NULL,cells=NULL)
  .PlotSubset(lidardata@Diffusion@ligand_conc,ligand_list,cells,"Exposure")

GetDecay <- function(lidardata,ligand_list=NULL,coupling=NULL,pval=NULL,n=NULL,ki=NULL,
                     p_column="p.value",r50=NULL,r01=NULL,r=NULL) {
  if (!missing(r)) {
    if (!missing(coupling)) stop("Supply only one of coupling and r")
    coupling <- r
  }
  df<-lidardata@Diffusion@ligand_diff
  if (!"coupling" %in% names(df) && "r" %in% names(df)) {
    names(df)[match("r", names(df))] <- "coupling"
  }
  if(!"Ligand"%in%names(df))stop("Decay results lack Ligand column")
  ids<-.PlotPick(ligand_list,as.character(df$Ligand),"fitted ligands")
  result<-df[match(ids,df$Ligand),,drop=FALSE]
  if(!is.null(pval)) {
    if(!p_column%in%c("p.value","p.empirical"))stop("Invalid p_column")
    .PlotNumber(pval,"pval");if(pval>1)stop("pval must be <=1")
    if(!p_column%in%names(result))stop("Requested P value is unavailable; enable permutation_p during fitting for p.empirical")
    result<-result[is.finite(result[[p_column]])&result[[p_column]]<pval,,drop=FALSE]
  }
  ranges<-list(coupling=coupling,n=n,ki=ki,r50=r50,r01=r01)
  for(field in names(ranges)) {
    limits<-ranges[[field]]
    if(is.null(limits))next
    if(!is.numeric(limits)||length(limits)!=2L||anyNA(limits)||limits[1]>limits[2])stop("Invalid range for ",field)
    if(!field%in%names(result))stop("Missing result column: ",field)
    result<-result[is.finite(result[[field]])&result[[field]]>=limits[1]&result[[field]]<=limits[2],,drop=FALSE]
  }
  attr(result,"fit_config")<-attr(df,"fit_config")
  permutations<-attr(df,"permutation_results")
  if(!is.null(permutations))attr(result,"permutation_results")<-permutations[intersect(result$Ligand,names(permutations))]
  result
}
SelectLigands <- function(lidardata,coupling=0,pval=.01,ligandlist=NULL,p_column="p.value",r=NULL) {
  if (!missing(r)) {
    if (!missing(coupling)) stop("Supply only one of coupling and r")
    coupling <- r
  }
  if(!is.numeric(coupling)||length(coupling)!=1L||!is.finite(coupling)||coupling< -1||coupling>1)stop("coupling must be between -1 and 1")
  df<-GetDecay(lidardata,ligandlist,pval=pval,p_column=p_column)
  as.character(df$Ligand[is.finite(df$coupling)&df$coupling>coupling])
}
# Retained as a standalone descriptive helper, never used for ligand selection.
CalcExprPct <- function(exp) {
  .PlotNeed("Matrix")
  if(length(dim(exp))!=2L||ncol(exp)==0L)stop("Expression must have at least one cell")
  .PlotFinite(exp,"Expression")
  Matrix::rowSums(exp>0)/ncol(exp)
}

.FeatureNames <- c(LigandExpr="Ligand expression",LigandExprImp="Ligand expression (imputed)",
  Response="Ligand response",Exposure="Ligand exposure",Availability="Receptor availability",
  AvailabilityImp="Receptor availability (imputed)",Score="LiDAR score",SenderScore="Sender score",ReceiverScore="Receiver score")
.FeatureAliases <- c(Ligand_exp="LigandExpr",Ligand_exp_imp="LigandExprImp",Ligand_act="Response",Ligand_conc="Exposure",
  Ligand_cores_receptor_exp="Availability",Ligand_cores_receptor_exp_imp="AvailabilityImp",Ligand_diff="Decay",
  Ligand_Diff="Decay",Ligand_conc_act_cor="Coupling")
.CanonicalFeatures <- function(features) {
  features<-as.character(features);match_old<-features%in%names(.FeatureAliases)
  features[match_old]<-unname(.FeatureAliases[features[match_old]])
  bad<-setdiff(features,c(names(.FeatureNames),"Decay","Coupling"))
  if(length(bad))stop("Unknown features: ",paste(bad,collapse=", "))
  unique(features)
}
GetFeatures <- function(lidardata,ligand_list=NULL,receptor_list=NULL,cells=NULL,
  features=c("LigandExpr","LigandExprImp","Response","Exposure","Availability","AvailabilityImp")) {
  features<-.CanonicalFeatures(features)
  if(any(features%in%c("Decay","Coupling")))stop("GetFeatures returns spatial matrices only")
  result<-setNames(vector("list",length(features)),features)
  for(f in features) result[[f]]<-switch(f,
    LigandExpr=GetLigandExpr(lidardata,ligand_list,cells,FALSE),
    LigandExprImp=GetLigandExpr(lidardata,ligand_list,cells,TRUE),
    Response=GetResponse(lidardata,ligand_list,cells),Exposure=GetExposure(lidardata,ligand_list,cells),
    Availability=GetAvailability(lidardata,ligand_list,receptor_list,cells,FALSE),
    AvailabilityImp=GetAvailability(lidardata,ligand_list,receptor_list,cells,TRUE),
    Score=.PlotSubset(lidardata@Lidarscore@aggregated_score,ligand_list,cells,"LiDAR score"),
    SenderScore=.PlotSubset(t(lidardata@Lidarscore@sender_score),ligand_list,cells,"Sender score"),
    ReceiverScore=.PlotSubset(t(lidardata@Lidarscore@receiver_score),ligand_list,cells,"Receiver score"))
  result
}
.PlotValues <- function(values,scale) {
  scale<-match.arg(scale,c("none","zscore"))
  if(!is.numeric(values))stop("Feature values must be numeric")
  .PlotFinite(values,"Feature values")
  result<-as.numeric(values)
  if(scale=="zscore") {
    sdv<-if(length(result)>1L)stats::sd(result) else 0
    # A constant vector has no standardized contrast; display zero and disclose it.
    if(sdv==0) {
      warning("Constant feature: Z-score display set to zero",call.=FALSE)
      result[]<-0
    } else result<-(result-mean(result))/sdv
  }
  names(result)<-names(values);result
}
.SpatialColourScale <- function(feature, col_palette, limits) {
  if (!is.null(col_palette))
    return(ggplot2::scale_color_viridis_c(option = col_palette, limits = limits))
  if (feature %in% c("LigandExpr", "LigandExprImp"))
    return(ggplot2::scale_color_viridis_c(option = "F", direction = -1, end = 0.95, limits = limits))
  if (feature %in% c("Score", "ReceiverScore"))
    return(ggplot2::scale_color_viridis_c(option = "G", direction = -1, begin = 0.3, limits = limits))
  if (feature %in% c("Exposure", "SenderScore", "Response")) {
    colours <- if (feature == "Response")
      c("#355C9A", "#8FB9D4", "#F3F3F3", "#D98C8C", "#A63D40")
    else
      c("#DCEAF4", "#A6CEE3", "#6BAED6", "#3182BD", "#08519C", "#08306B")
    return(ggplot2::scale_color_gradientn(colours = grDevices::colorRampPalette(colours)(100), limits = limits))
  }
  ggplot2::scale_color_viridis_c(option = "magma", limits = limits)
}
PlotFeature <- function(values,cell_positions,show_border=FALSE,scale=c("zscore","none"),
                         title=NULL,labs_params=list(color="Feature"),point_size=.05,
                         col_palette=NULL,reverse_y=FALSE,color_limits=NULL,feature=NULL,zscore=TRUE,ratio=0.8) {
  .PlotNumber(ratio,"ratio",TRUE)
  .PlotFlag(zscore,"zscore")
  if (missing(scale)) scale <- if (zscore) "zscore" else "none"
  .PlotNeed("ggplot2");scale<-match.arg(scale);.PlotFlag(reverse_y,"reverse_y");.PlotFlag(show_border,"show_border")
  .PlotNumber(point_size,"point_size",TRUE)
  if(!all(c("cell","x","y")%in%names(cell_positions)))stop("Coordinates require cell/x/y")
  .PlotNames(cell_positions$cell,"Coordinate cells");.PlotNames(names(values),"Feature cells")
  if(!length(values)||!setequal(names(values),cell_positions$cell))stop("Feature and coordinate cells must match exactly")
  if(!all(vapply(cell_positions[c("x","y")],is.numeric,logical(1))))stop("Coordinates must be numeric")
  .PlotFinite(as.matrix(cell_positions[c("x","y")]),"Coordinates")
  values<-values[match(cell_positions$cell,names(values))]
  values<-.PlotValues(values,scale)
  if (is.null(feature)) feature <- "Feature" else {
    if (length(feature) != 1L) stop("Supply one spatial feature")
    feature <- .CanonicalFeatures(feature)
    if (!feature %in% names(.FeatureNames)) stop("Expected a spatial feature")
  }
  if(!is.null(color_limits)&&(!is.numeric(color_limits)||length(color_limits)!=2L||any(!is.finite(color_limits))||color_limits[1]>color_limits[2]))stop("Invalid color_limits")
  plot_data<-cell_positions;plot_data$feature<-as.numeric(values)
  if(scale=="zscore")labs_params$color<-paste0(if(is.null(labs_params$color))"Feature" else labs_params$color," (Z-score)")
  p<-ggplot2::ggplot(plot_data,ggplot2::aes(x=x,y=y,color=feature))+
    ggplot2::geom_point(size=point_size)+.SpatialColourScale(feature,col_palette,color_limits)+
    ggplot2::coord_fixed(ratio=ratio)+theme_Lidar2(border=show_border)+do.call(ggplot2::labs,labs_params)
  if(!is.null(title))p<-p+ggplot2::ggtitle(title)
  if(reverse_y)p<-p+ggplot2::scale_y_reverse()
  p
}
.DecayConfig <- function(lidardata) {
  cfg<-attr(lidardata@Diffusion@ligand_diff,"fit_config")
  if(is.null(cfg)||is.null(cfg$B)||is.null(cfg$bin)||is.null(cfg$n_bins)||is.null(cfg$imp))
    stop("Missing fit_config: use FitPropagation from the revised publication script")
  .PlotNumber(cfg$B,"fitted B",TRUE);cfg
}
.ContinuousDecay <- function(ki,n,d,B) {
  if(any(!is.finite(c(ki,n,B)))||any(c(ki,n,B)<=0))stop("Invalid fitted decay parameters")
  z<-pmin(d/B,1)
  w<--expm1(n*log(z))*stats::plogis(-n*(log(z)-log(ki)))
  w[z==0]<-1;w[z>=1]<-0;w
}
.PLabel <- function(row) {
  label<-paste0("Exposure-response coupling = ",formatC(row$coupling,digits=3,format="f"),"; Pearson P = ",format.pval(row$p.value,digits=3))
  if("p.empirical"%in%names(row)) {
    empirical<-if(is.finite(row$p.empirical))format.pval(row$p.empirical,digits=3) else "NA (incomplete)"
    label<-paste0(label,"\nempirical P = ",empirical)
  }
  label
}
PlotDecay <- function(lidardata,ligand,show_ranges=TRUE,distance_unit=NULL,n_points=501L) {
  .PlotNeed("ggplot2");.PlotFlag(show_ranges,"show_ranges")
  if(length(ligand)!=1L)stop("Specify one ligand")
  if(length(n_points)!=1L||!is.finite(n_points)||n_points<2||n_points!=floor(n_points))stop("n_points must be an integer >=2")
  row<-GetDecay(lidardata,ligand);cfg<-.DecayConfig(lidardata)
  d<-seq(0,cfg$B,length.out=n_points)
  dat<-data.frame(distance=d,weight=.ContinuousDecay(row$ki,row$n,d,cfg$B))
  label<-if(is.null(distance_unit))"Distance (input coordinate units)" else paste0("Distance (",distance_unit,")")
  p<-ggplot2::ggplot(dat,ggplot2::aes(x=distance,y=weight))+ggplot2::geom_line()+
    ggplot2::labs(title=ligand,subtitle=.PLabel(row),x=label,y="Decay weight",
      caption="Continuous fitted curve.\nExposure uses left-edge distance-bin weights.")+theme_Lidar1()
  if(show_ranges) {
    if(!all(c("r50","r01")%in%names(row)))stop("Fitted results lack r50/r01")
    marks<-data.frame(distance=c(row$r50,row$r01),label=c("r50","r01"))
    .PlotFinite(marks$distance,"Decay ranges")
    p<-p+ggplot2::geom_vline(data=marks,ggplot2::aes(xintercept=distance,linetype=label),colour="grey40")+
      ggplot2::labs(linetype="Continuous range")
  }
  p
}
.CouplingData <- function(lidardata,ligand,cells=NULL) {
  cfg<-.DecayConfig(lidardata)
  GetDecay(lidardata,ligand)
  ids<-.PlotPick(cells,lidardata@cellname,"spatial units")
  # Fit receiver definition always uses all cognate receptors and the fitted imp setting.
  availability<-GetAvailability(lidardata,ligand,cells=lidardata@cellname,imp=cfg$imp)
  .PlotFinite(availability,"Availability")
  full_receivers<-rownames(availability)[as.numeric(availability[,ligand])>0]
  receivers<-ids[ids%in%full_receivers]
  response<-GetResponse(lidardata,ligand,receivers)
  exposure<-GetExposure(lidardata,ligand,receivers)
  dat<-data.frame(cell=receivers,response=as.numeric(response[,ligand]),exposure=as.numeric(exposure[,ligand]))
  .PlotFinite(dat$response,"Response");.PlotFinite(dat$exposure,"Exposure")
  list(data=dat,full=setequal(receivers,full_receivers))
}
PlotCoupling <- function(lidardata,ligand,cells=NULL,point_size=.2) {
  .PlotNeed("ggplot2");.PlotNumber(point_size,"point_size",TRUE)
  if(length(ligand)!=1L)stop("Specify one ligand")
  a<-.CouplingData(lidardata,ligand,cells);dat<-a$data
  if(nrow(dat)<3L||stats::sd(dat$response)==0||stats::sd(dat$exposure)==0)
    stop("Coupling plot requires >=3 receivers with nonconstant response and exposure")
  ct<-stats::cor.test(dat$exposure,dat$response,method="pearson")
  row<-GetDecay(lidardata,ligand)
  if(a$full) {
    if(!isTRUE(all.equal(unname(ct$estimate),as.numeric(row$coupling),tolerance=1e-7)))
      stop("Stored fit and current response/exposure disagree; refit before plotting")
    subtitle<-.PLabel(row)
  } else {
    subtitle<-paste0("Selected receivers: exposure-response coupling = ",formatC(unname(ct$estimate),digits=3,format="f"),
                     "\ndescriptive Pearson P = ",format.pval(ct$p.value,digits=3))
  }
  ggplot2::ggplot(dat,ggplot2::aes(x=exposure,y=response))+
    ggplot2::geom_point(size=point_size,alpha=.6)+
    ggplot2::geom_smooth(method="lm",se=FALSE,colour="navy",formula=y~x)+
    ggplot2::labs(title=ligand,subtitle=subtitle,x="Ligand exposure",y="Ligand response",
      caption=if(a$full)"Potential receivers used for fitting." else "Subset correlation.\nFull-fit empirical P value is not applied to this subset.")+
    theme_Lidar1()
}

ExportPlots <- function(lidardata,cells=NULL,ligand_list=NULL,receptor_list=NULL,
  features=c("LigandExpr","LigandExprImp","Response","Exposure","Availability","AvailabilityImp","Decay","Coupling"),
  out_dir=NULL,scale=c("zscore","none"),show_border=FALSE,col_palette=NULL,reverse_y=FALSE,
  shared_scales=TRUE,distance_unit=NULL,format=c("jpg","png","pdf"),dpi=300,
  point_size=.05,overwrite=FALSE,zscore=TRUE,ratio=0.8) {
  .PlotNeed("ggplot2");.PlotNeed("patchwork")
  .PlotNumber(ratio,"ratio",TRUE)
  .PlotFlag(zscore,"zscore")
  if (missing(scale)) scale <- if (zscore) "zscore" else "none"
  scale<-match.arg(scale);format<-match.arg(format)
  for(v in c("shared_scales","overwrite","show_border","reverse_y")) .PlotFlag(get(v),v)
  .PlotNumber(dpi,"dpi",TRUE);.PlotNumber(point_size,"point_size",TRUE)
  features<-.CanonicalFeatures(features);if(!length(features))stop("Choose at least one feature")
  coords<-GetCoords(lidardata,cells);ids<-coords$cell
  ligands<-.PlotPick(ligand_list,lidardata@DB@ligandname,"database ligands")
  if(is.null(out_dir))out_dir<-paste0("LiDAR_",base::format(Sys.time(),"%Y%m%d_%H%M%S"))
  if(!dir.exists(out_dir)&&!dir.create(out_dir,recursive=TRUE))stop("Cannot create output directory")
  out_dir<-normalizePath(out_dir,winslash="/",mustWork=TRUE)
  manifest<-list();messages<-list();prepared<-list();limits<-list()
  record<-function(l,f,status,message="",path="") {
    messages[[length(messages)+1L]]<<-data.frame(Ligand=l,feature=f,status=status,message=message,path=path,stringsAsFactors=FALSE)
  }
  spatial<-intersect(features,names(.FeatureNames))
  # Load only requested feature matrices. Internal coupling inputs never add panels.
  for(f in spatial) {
    mat<-tryCatch(GetFeatures(lidardata,receptor_list=receptor_list,cells=ids,features=f)[[f]],error=function(e)e)
    if(inherits(mat,"error")) {
      for(l in ligands)record(l,f,"skipped",conditionMessage(mat))
      next
    }
    available<-intersect(ligands,colnames(mat)); vectors<-list()
    for(l in setdiff(ligands,available))record(l,f,"skipped","Feature unavailable for this ligand")
    for(l in available) {
      raw<-setNames(as.numeric(mat[,l]),rownames(mat))
      values<-tryCatch(.PlotValues(raw,scale),error=function(e)e)
      if(inherits(values,"error"))record(l,f,"skipped",conditionMessage(values)) else vectors[[l]]<-values
    }
    prepared[[f]]<-vectors
    if(shared_scales&&length(vectors))limits[[f]]<-range(unlist(vectors,use.names=FALSE))
  }
  safe<-gsub('[<>:"/\\\\|?*[:cntrl:]]',"_",ligands)
  safe<-sub("[. ]+$","",safe);safe[!nzchar(safe)]<-"ligand"
  safe<-make.unique(paste0("ligand_",safe),sep="_")
  for(i in seq_along(ligands)) {
    l<-ligands[i];panels<-list()
    for(f in features) {
      if(f%in%spatial) {
        values<-prepared[[f]][[l]]
        if(is.null(values))next
        label<-.FeatureNames[[f]]
        if(scale=="zscore")label<-paste0(label," (Z-score)")
        p<-tryCatch(PlotFeature(values,coords,show_border=show_border,scale="none",title=paste(l,.FeatureNames[[f]]),
          labs_params=list(color=label),point_size=point_size,col_palette=col_palette,reverse_y=reverse_y,color_limits=limits[[f]],feature=f,ratio=ratio),error=function(e)e)
      } else {
        p<-tryCatch(if(f=="Decay")PlotDecay(lidardata,l,distance_unit=distance_unit) else PlotCoupling(lidardata,l,cells),error=function(e)e)
      }
      if(inherits(p,"error"))record(l,f,"skipped",conditionMessage(p)) else {
        # Compact bottom legends leave enough horizontal room in a 5-inch panel.
        panels[[f]]<-p+ggplot2::theme(legend.position="bottom",
          legend.title=ggplot2::element_text(size=10),legend.text=ggplot2::element_text(size=10),
          plot.title.position="plot",plot.caption.position="plot",plot.margin=ggplot2::margin(14,14,14,14))
      }
    }
    if(!length(panels)) {record(l,"export","skipped","No valid panels");next}
    path<-file.path(out_dir,paste0(safe[i],".",format))
    if(file.exists(path)&&!overwrite) {record(l,"export","skipped","File exists; set overwrite=TRUE",path);next}
    err<-tryCatch({
      # Isolate panel layouts: fixed-aspect spatial plots must not squeeze the
      # decay/coupling axes or impose their legend widths on unrelated panels.
      layout_panels<-lapply(panels,function(p)patchwork::wrap_elements(full=p))
      combined<-patchwork::wrap_plots(layout_panels,ncol=min(2L,length(panels)),
        widths=rep(1,min(2L,length(panels))),heights=rep(1,ceiling(length(panels)/2)))
      ggplot2::ggsave(path,plot=combined,width=5*min(2L,length(panels)),height=4*ceiling(length(panels)/2),units="in",dpi=dpi,bg="white")
      NULL
    },error=function(e)conditionMessage(e))
    if(is.null(err)) {
      record(l,"export","saved",paste(names(panels),collapse=", "),path)
      manifest[[l]]<-list(path=path,features=names(panels))
    } else record(l,"export","failed",err,path)
  }
  log<-if(length(messages))do.call(rbind,messages) else data.frame(Ligand=character(),feature=character(),status=character(),message=character(),path=character())
  if(any(log$status%in%c("skipped","failed")))warning("Some panels/files were skipped or failed; inspect returned $log",call.=FALSE)
  message("Saved ",length(manifest)," ligand plot files to ",out_dir)
  invisible(list(files=manifest,log=log,color_limits=limits,scale=scale,shared_scales=shared_scales,out_dir=out_dir))
}

# Compatible old extractor names; pct is intentionally not an accepted argument.
get_cell_positions<-GetCoords
get_ligand_exp<-GetLigandExpr
get_receptor_exp<-GetReceptorExpr
get_ligand_cores_receptor_exp<-GetAvailability
get_ligand_act<-GetResponse
get_ligand_conc<-GetExposure
get_ligand_diff<-GetDecay
get_succ_ligand<-SelectLigands
calc_pct<-CalcExprPct
make_feature_matrices<-GetFeatures
# Plot signatures changed deliberately: use PlotFeature/PlotDecay/PlotCoupling/ExportPlots.

# Original 120-color palette; retain colors100 as a compatibility alias.
lidar_colors <- c(
  
  ### Blue --
  "#4b6aa8","#537eb7","#3674a2","#2f5f99","#5a7fc6",
  "#3a6fb0","#4c82c4","#3f6e9e","#2f4e88","#6b8fc9",
  
  ### Cyan --
  "#3ca0cf","#4490c4","#63a3b8","#61bada","#76a2be",
  "#6fb9c9","#7fc8d6","#4fb1d3","#5aa9c7","#8bcfe0",
  
  ### Teal --
  "#5aa7a5","#4d9896","#3c8a88","#6fb6b3","#7ec6c2",
  "#58a3a0","#6fbfba","#87d1cc","#4f8f8b","#76c3bf",
  
  ### Green --
  "#408444","#5f9f7f","#64a776","#7ab089","#83ab8e",
  "#6aa96f","#74b07c","#8bc48f","#5b9c60","#9ed3a3",
  
  ### Light green --
  "#a9d3a2","#b6d9a8","#c2e1b5","#9fcd9a","#d0e7c3",
  "#b8d9b4","#c6e4c0","#d8efd2","#a8d0a3","#e2f3dd",
  
  ### Yellow --
  "#ece399","#e6e2a3","#efd9a1","#e9d277","#d5bb72",
  "#e3c986","#f1e29b","#d9c768","#e7d483","#f3e9ad",
  
  ### Orange --
  "#e6b884","#de8b36","#d69a55","#d69971","#e3a472",
  "#e69b63","#d9874b","#e7a979","#f0b98a","#cc7f73",
  
  ### Red --
  "#df5734","#da6f6d","#d25774","#b95055","#c84f4f",
  "#d9625a","#c94343","#e0746e","#b84242","#d86b6b",
  
  ### Pink --
  "#ebb1a4","#e29eaf","#d69ab0","#c376a7","#b85292",
  "#d992b5","#c96aa0","#e3a7c5","#f0bcd2","#b96a9f",
  
  ### Purple --
  "#6c408e","#8d689d","#92699e","#a44e89","#9a70a8",
  "#7e5ca1","#8c6ab3","#6b4f8c","#9d7bb6","#b08ac8",
  
  ### Brown --
  "#c4612f","#b05545","#bc7b63","#a78982","#9f8d89",
  "#b06a4d","#a25b3d","#c08b6d","#8f6d5e","#7a5849",
  
  ### Grey --
  "#696a6c","#8a8c90","#a5a9b0","#bfc2c7","#d4d6db",
  "#9c9ea3","#7f8286","#b5b8bd","#cfd2d6","#e2e5e8"
  
)




colors100 <- lidar_colors

