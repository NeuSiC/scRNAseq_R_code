#Rename column names, making them unique
renameCols <- function(cms, names) {
  if(!length(cms)==length(names)) stop("Names must match number of count matrices.")
  
  mapply(function(c, n) {
    colnames(c) <- lapply(colnames(c), function(cn) paste0(n,cn)) %>% unlist
    return(c)
  }, c=cms, n=names, SIMPLIFY = F)
}

#Mitochondrial fraction
mitoFraction <- function(con, species="human") {
  if(species=="human") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("MT-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else if(species=="mouse") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("mt-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else stop("Species must either be 'human' or 'mouse'.")
}



addEmbeddingP2Web <- function(p2, con, embedding=NULL, name="UMAP") {
  if(is.null(embedding)) embedding <- con$embedding
  
  if(identical(dim(p2$originalP2object$embeddings$PCA[[1]]),dim(embedding))) {
    p2$originalP2object$embeddings$PCA[[name]] <- embedding
    return(p2)
  } else {
    stop("The embedding dimensions of the p2.web object and the input object are not identical.")
  }
}

embedUMAP <- function(con,
                      min.dist=0.01,
                      spread=15,
                      min.prob.lower=1e-7,
                      method=leiden.community,
                      resolution=1,
                      min.group.size=25) {
  message("Creating UMAP embedding...")
  con$embedGraph(method="UMAP", 
                 min.dist=min.dist, 
                 spread=spread,
                 min.prob.lower=min.prob.lower)
  
  message("Estimating clusters...")
  con$findCommunities(method=leiden.community, resolution=resolution, min.group.size=min.group.size)
  
  return(con)
}

buildConosGraph <- function(con,
                            k.conos=15, 
                            k.self=15, 
                            space='PCA', 
                            ncomps=40,
                            n.odgenes=2e3,
                            matching.method='mNN', 
                            metric='angular', 
                            score.component.variance=T,
                            alignment.strength=0,
                            min.dist=0.01, 
                            spread=15,
                            min.prob.lower=1e-3,
                            resolution=1,
                            min.group.size=25) {
  message("Building graph...")
  con$buildGraph(k=k.conos, 
                 k.self=k.self, 
                 space=space, 
                 ncomps=ncomps, 
                 n.odgenes=n.odgenes, 
                 matching.method=matching.method, 
                 metric=metric, 
                 verbose=T, 
                 score.component.variance=score.component.variance,
                 alignment.strength=alignment.strength)
  
  embedUMAP(con=con,
            min.dist=min.dist,
            spread=spread,
            min.prob.lower=min.prob.lower,
            method=leiden.community,
            resolution=resolution,
            min.group.size=min.group.size)
  
  return(con)
}

quickConos <- function(cms, 
                       sample.names,
                       n.cores.p2,
                       n.cores.con,
                       n.odgenes=3e3, 
                       nPcs = 50, 
                       k.p2 = 30, 
                       perplexity = 50, 
                       log.scale = TRUE, 
                       trim = 10, 
                       keep.genes = NULL, 
                       min.cells.per.gene = 3, 
                       min.transcripts.per.cell = 200, 
                       get.largevis = F, 
                       get.tsne = F, 
                       make.geneknn = F,
                       k.conos=15, 
                       k.self=30, 
                       space='PCA', 
                       ncomps=40, 
                       matching.method='mNN', 
                       metric='angular', 
                       score.component.variance=T,
                       alignment.strength=0,
                       min.dist=0.01, 
                       spread=15) {
  if(length(cms)==length(sample.names)) {
    if(any(is.na(sample.names))) stop("Names contains NAs")
    
    message("Performing P2 processing...")
    panel.preprocessed <- lapply(cms, function(x) basicP2proc(x, n.cores = n.cores.p2,
                                                              n.odgenes = n.odgenes, 
                                                              nPcs = nPcs,
                                                              k = k.p2, 
                                                              perplexity = perplexity, 
                                                              log.scale = log.scale, 
                                                              trim = trim, 
                                                              keep.genes = keep.genes, 
                                                              min.cells.per.gene = min.cells.per.gene, 
                                                              min.transcripts.per.cell = min.transcripts.per.cell, 
                                                              get.largevis = get.largevis, 
                                                              get.tsne = get.tsne, 
                                                              make.geneknn = make.geneknn))
    
    names(panel.preprocessed) = sample.names
    con <- Conos$new(panel.preprocessed, n.cores=n.cores.con)
    
    con <- buildConosGraph(con=con,
                           k.conos=k.conos, 
                           k.self=k.self, 
                           space=space, 
                           ncomps=ncomps, 
                           n.odgenes=n.odgenes, 
                           matching.method=matching.method, 
                           metric=metric, 
                           score.component.variance=score.component.variance,
                           alignment.strength=alignment.strength,
                           min.dist=min.dist, 
                           spread=spread)
    
    return(list(con=con, panel.preprocessed=panel.preprocessed))
  } else {
    stop("Sample names must match number of count matrices.")
  }
  
}

collapseAnnotation <- function(anno, label) {
  anno %<>% factor
  idx <- grepl(label,levels(anno))
  cat(paste0("Collapsing ",sum(idx)," labels containing '",label,"' in their name into one label.\n"))
  levels(anno)[idx] <- c(label)
  anno %<>% factor
  return(anno)
}

getConosDepth <- function(con) {
  lapply(con$samples, function(d) d$depth) %>% unlist %>% setNames(.,(strsplit(names(.), ".", T) %>% 
                                                                        sapply(function(d) d[2])))
}

getConosCluster <- function(con, name="leiden") {
  con$clusters[[name]]$groups
}

plotDotMap <- function (markers, 
                        count.matrix, 
                        annotation, 
                        marker.colour="black",
                        cluster.colour="black",
                        text.angle = 45, 
                        gene.order = NULL, 
                        cols = c("blue", "red"),
                        col.min = -2.5,
                        col.max = 2.5,
                        dot.min = 0,
                        dot.scale = 6,
                        scale.by = "radius",
                        scale.min = NA,
                        scale.max = NA,
                        verbose=T) {
  scale.func <- switch(scale.by, 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if(verbose) cat("Plotting .")
  
  if(!is.character(markers)) stop("'markers' must be a character vector.")
  
  missing.markers <- setdiff(markers, colnames(count.matrix))
  if(length(missing.markers)>0) {
    cat("Not all markers are in 'count.matrix'. The following are missing:\n",paste(missing.markers, collapse=" "),"\n")
    stop("Please update 'markers'.")
  }
  
  # From CellAnnotatoR:::plotExpressionViolinMap, should be exchanged with generic function
  p.df <- lapply(markers, function(g) data.frame(Expr = count.matrix[names(annotation), g], Type = annotation, Gene = g)) %>% Reduce(rbind, .)
  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  } else {
    gene.order <- NULL
  }
  
  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene = factor(as.character(Gene), 
                                          levels = gene.order))
  }
  
  # Adapted from Seurat:::DotPlot
  if(verbose) cat(".")
  data.plot <- levels(annotation) %>% lapply(function(t) {
    markers %>% lapply(function(g) {
      df <- p.df %>% filter(Type==t, Gene==g)
      pct.exp <- sum(df$Expr>0)/dim(df)[1]*100
      avg.exp <- mean(df$Expr[df$Expr>0])
      res <- data.frame(gene=g,
                        pct.exp=pct.exp,
                        avg.exp=avg.exp)
      return(res)
    }) %>% Reduce(rbind, .)
  }) %>% 
    setNames(levels(annotation)) %>%
    bind_rows(., .id="cluster")
  
  data.plot$cluster %<>% factor(., levels=rev(unique(.)))
  
  data.plot %<>% arrange(gene)
  
  data.plot$avg.exp.scaled <- data.plot$gene %>% unique %>% sapply(function(g) {
    data.plot %>% .[.$gene == g, 'avg.exp'] %>% 
      scale %>% 
      MinMax(min = col.min, max = col.max)
  }) %>% unlist %>% as.numeric
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  
  cluster.colour %<>% rev
  
  plot <- ggplot(data.plot, aes_string("gene", "cluster")) +
    geom_point(aes_string(size = "pct.exp", color = "avg.exp.scaled")) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.text.x = element_text(angle=text.angle, hjust = 1, colour=marker.colour),
          axis.text.y = element_text(colour=cluster.colour),
          panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(size = guide_legend(title = 'Percent expressed'), color = guide_colorbar(title = 'Average expression')) +
    labs(x = 'Marker', y = 'Cluster') +
    scale_color_gradient(low = cols[1], high = cols[2])
  if(verbose) cat(" done!")
  return(plot)
}

renameAnnotation <- function(annotation, old, new) {
  if(!is.factor(annotation)) stop("Annotation must be a factor.")
  
  levels(annotation)[levels(annotation) %in% old] <- new
  
  return(annotation)
}

dotSize <- function(size, alpha=1) {
  ggplot2::guides(colour = guide_legend(override.aes = list(size=size,
                                                            alpha=alpha)))
}

checkDims <- function(cm, con) {
  cat("Dimensions of cm : ",paste((dim(cm)), collapse=" "),"\n")
  cat("Dimensions of con: ",paste((dim(con$embedding)), collapse=" "),"\n")
  
  if(dim(cm)[2]!=dim(con$embedding)[1])
    stop("Dimensions don't match.")
  
  message("All OK!")
}


#dotSize
dotSize <- function(size, alpha=1) {
  ggplot2::guides(colour = guide_legend(override.aes = list(size=size,
                                                            alpha=alpha)))}

