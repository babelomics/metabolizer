############################################################ 
#                     FUNCTIONS FOR                        #
#         To calculate module activities using             #
#       the  gene expression data (the Metabolizer)        #
#         http://metabolizer.babelomics.org/               #
#                                                          #
# if you use the Metabolizer tool (R code otr webserver)   #
# please cite it as; Please, use the following             #
# convention to cite Metabolizer tool:                     #
#                                                          #
# Cubuk, C., Hidalgo, M., Amadoz, A., Rian, K.,            #
# Salavert, F., Pujana, M., Mateo, F., Herranz, C.,        #
# Carbonell-Caballero, J. and Dopazo, J. (2019).           #
# Differential metabolic activity and discovery of         #
# therapeutic targets using summarized metabolic           # 
# pathway models.                                          #
# npj Systems Biology and Applications, 5(1),              #
# doi: 10.1038/s41540-019-0087-2.                          #
# https://www.nature.com/articles/s41540-019-0087-2        #
#                                                          #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#                    2016-2019                             #
#          @ Sevilla and Valencia, Spain                   #
#                                                          #
############################################################ 

metabolizer <- function(hsa_module_data, combat.vals, output_folder=NULL, saveName=NULL, onesample=F, expbased=T, fluxbased=F, default_value=.5, moduleinfo=NULL){
  
  # has_module_data module data object. This object is named as "hsa_module_data" for all the organisms that are modeled in the Metabolizer.
  # combat.vals gene expression matrix
  # results <- metabolizer(hsa_module_data, combat.vals, output_folder=NULL, saveName=NULL, onesample=F, expbased=T, fluxbased=F, default_value=0.5, moduleinfo="~/metabolizer/files/moduleinfo.RData")
  
  load(moduleinfo)
  
  all_metabolites <- get.metabolite.list(hsa_module_data)
  all_module_genes_vec <- get.All.module.genes(hsa_module_data)
  all_module_rxn_vec <- get.All.module.rxns(hsa_module_data) # buna bak
  rxn_gene_mat <- get.rxn.gene.matrix(hsa_module_data)
  gene_exp_mat <- get.gene.exp.of.module.genes(combat.vals, all_module_genes_vec, min.exp=0.0001, onesample=onesample)
  rxn_vals_mat <- get.RXNvals.from.genes(all_module_rxn_vec, gene_exp_mat, rxn_gene_mat)
  metabolite_matrix <- mat.or.vec(nr = length(all_metabolites), nc = ncol(rxn_vals_mat))
  rownames(metabolite_matrix) <- all_metabolites
  metabolite_matrix[,] <- 1
  moduleActs <- methyways(hsa_module_data,rxn_vals_mat=rxn_vals_mat,expbased=expbased,fluxbased=fluxbased,verbose = F,default_value=default_value, metabolitematrix = metabolite_matrix)
  
  moduleActsList <- sapply(moduleActs, function(x){ x[[1]][1] })
  moduleActsMatrix <- do.call("rbind",moduleActsList)
  if(!is.null(moduleinfo)){
    # rownames(moduleActsMatrix)<- paste0(rownames(moduleActsMatrix),": ",moduleinfo$ModuleDescription[match(rownames(moduleActsMatrix) ,moduleinfo$Module.ID)])
    moduleActsMatrix <- cbind(moduleActsMatrix, moduleinfo)
  }
  
  nodevalslist <- sapply(moduleActs, function(x){ x[[1]][3] })
  nodevalsmatrix <- do.call("rbind",nodevalslist)
  allnodes <- unique(rownames(nodevalsmatrix))
  allnodes <- allnodes[grep("R",allnodes)]
  nodevalsmatrix <- nodevalsmatrix[match(allnodes, rownames(nodevalsmatrix)),]
  
  if(!is.null(output_folder)){  
    write.table(x=moduleActsMatrix,file=paste0(output_folder,"/",saveName,"_paths_vals.txt"),quote=F,sep="\t")
  }
  
  message("The citation for the Metabolizer tool: \n-> https://www.nature.com/articles/s41540-019-0087-2 \n-> doi: 10.1038/s41540-019-0087-2")
  
  return(list("moduleActivities"=moduleActsMatrix,
              "ReactionNodeValues"=nodevalsmatrix,
              "ModuleGenes"=all_module_genes_vec,
              "GeneExp"=gene_exp_mat))
}


methyways <- function(hsa_module_data,rxn_vals_mat,expbased=T,fluxbased=F,flux_dist, verbose=T, default_value=0.5, metabolitematrix=NULL){
  
  library(igraph)
  res<-lapply(hsa_module_data,function(x){ 
    
    graphobject<-x$graphobject
    KEGG_met_path_node<-x$KEGG_met_path_node
    
    if(expbased || !fluxbased ){
      nodes.vals<-node_vals(graphobject,KEGG_met_path_node=KEGG_met_path_node,rxn_vals_mat,default_value)
    }else{
      # calculates based on flux value of each reaction
      message("fluxbased=T; needs node_vals_for_flux from source(\"~/functions.unnecessary.1403209.R\")")
      nodes.vals<-node_vals_for_flux(graphobject,KEGG_met_path_node,flux_dist)
    }
    
    if(!is.null(graphobject$reduced_igraph)){
      subgraph<-graphobject$reduced_igraph
      ininodes<-graphobject$reduced_initial_node
      endnode<-graphobject$reduced_end_node
    }else{
      subgraph<-graphobject$igraph
      ininodes<-graphobject$initial_node
      endnode<-graphobject$end_node 
    }
    
    if(!is.null(metabolitematrix)){
      nodes.vals <- rbind(nodes.vals, metabolitematrix)
    }
    
    path.val<-path.value(nodes.vals=nodes.vals, subgraph ,ininodes,endnode)
    if(verbose){ cat(graphobject$module_name,"...DONE \n") }
    return(list(path.val))
  })
  
  return(res)
}

path.value<-function( nodes.vals, subgraph, ininodes, endnode, method="pond", maxnum = 100, tol = 0.000001, divide=F, response.tol = 0 ){
  # Initialize lists
  ready <- ininodes
  processed <- list()
  
  # Initialize node values
  node.signal <- matrix(NA, ncol=ncol(nodes.vals), nrow = length(V(subgraph)), dimnames = list(V(subgraph)$name, colnames(nodes.vals)))
  endnode.signal.dif <- 10
  
  num <- 0
  reached_last <- F
  while( length(ready) > 0 && num <= maxnum){
    num <- num + 1
    actnode <- ready[[1]]
    old.signal <- node.signal[actnode,]
    
    # Compute node signal
    if(divide && actnode != endnode){
      nfol <- length(incident(subgraph, actnode, mode="out"))
    }else{
      nfol <- 1
    }
    # print(nfol)
    node.signal[actnode,] <- compute.node.signal2(actnode, nodes.vals[actnode,], node.signal, subgraph, method, response.tol) / nfol
    
    # Transmit signal
    nextnodes <- get.edgelist(subgraph)[incident(subgraph, actnode, mode="out"),2]
    dif <- old.signal - node.signal[actnode,]
    
    if(actnode==endnode){
      reached_last <- T
      if(!all(is.na(dif)))
        endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2)))
      #num <- num+1
    }
    if(all(is.na(old.signal)) || endnode.signal.dif[length(endnode.signal.dif)] > tol )
      ready <- unique(c(ready, nextnodes))
    ready <- ready[-1]
  }
  if(reached_last==F){
    endnode.signal.dif <- NA
  }
  return(list(node.signal[endnode,], endnode.signal.dif,nodes.vals))
}


compute.node.signal2<-function(actnode, node.val, node.signal, subgraph, method="pond", response.tol = 0){
  
  incis <- incident(subgraph, actnode, mode="in")
  
  if(length(incis)==0){
    signal <- rep(1, length(node.val))
    
  } else {
    
    # get activators and inhibitors signal
    prevs <- get.edgelist(subgraph)[incis,1]
    input_signals <- node.signal[prevs,,drop=F]
    nas <- is.na(input_signals[,1])
    prevs <- prevs[!nas]
    incis <- incis[!nas]
    input_signals <- input_signals[!nas,,drop=F]
    typeincis <- E(subgraph)$relation[incis]
    activators <- typeincis==1
    nactivators <- sum(activators)
    inhibitors <- typeincis==-1
    ninhibitors <- sum(inhibitors)
    activator_signals <- input_signals[activators,,drop=F]
    inhibitor_signals <- input_signals[inhibitors,,drop=F]
    
    if( method == "sum"){
      s1 <- prettyifelse(nactivators>0, colSums(activator_signals), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, colSums(inhibitor_signals), rep(0,length(node.val)))
      signal <- s1-s2
    }
    else if( method == "pond"){
      s1 <- prettyifelse(nactivators>0, apply(1- activator_signals, 2, prod), rep(0,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, apply(1- inhibitor_signals, 2, prod), rep(1,length(node.val)))
      signal <- (1-s1)*s2
    }
    else if( method == "min"){
      s1 <- prettyifelse(nactivators>0, apply(activator_signals,2,min), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, 1-apply(inhibitor_signals,2,max), rep(1,length(node.val)))
      signal <- s1*s2
    }
    else {
      stop("Unknown propagation rule")
    }
    
    # If signal too low, signal do not propagate
    if(sum(nas) == 0 && signal < response.tol)
      signal <- rep(0,length(node.val))
    
  }
  
  signal[signal>1] <- 1
  signal[signal<0] <- 0
  signal <- signal*node.val
  
  return(signal)
}

prettyifelse<-function(test,una,olaotra){
  if(test){
    return(una)
  } else {
    return(olaotra)
  }
}


# this function checks and ignores missing data (gene exp or node val)
# does not require prior functions such as add missing data
node_vals<-function(graphobject,KEGG_met_path_node,rxn_vals_mat,default_value){
  
  # graphobject is output of moduleSIFandGraph
  # KEGG_met_path_node is output of KEGG_Met_path_parse
  # rxn_vals_mat is nrow=n rxns ncol= m samples nxm matrix
  
  if(!is.null(graphobject$reduced_SIF)){SIF<-graphobject$reduced_SIF}else{SIF<-graphobject$SIF}
  
  entries<-unique(c(SIF[,1],SIF[,3]))
  nodes.vals.total<-data.frame(name=seq(1:length(entries)))
  
  idxleft<-grep("R",SIF[,1])
  SIFleft<-SIF[,1][idxleft]
  idxright<-grep("R",SIF[,3])
  SIFright<-SIF[,3][idxright]
  reactions<-unique(c(SIFleft,SIFright))
  
  y<-sapply(reactions,function(x) {strsplit(x,split=",")[[1]]})
  if(all(reactions == y)){y<-as.list(y)}
  # fake data
  # y[1]<-list(c("R1+R2","R3","R4")) #c
  # y[2]<-list(c("R5+R6")) # a
  # y[3]<-list(c("R7","R8")) #h
  # y[4]<-list(c("R9")) #e
  # y[5]<-list(c("R19+R20+R21","R22+R23")) #b
  # y[6]<-list(c("R24+R25+R26","R27+R28","R29","R30")) #d
  
  y1 <- sapply(y,function(node) {
    
    #complex_enzyme_rxn
    complex_enzyme_rxn <-node[grep("[+]",node)]
    
    if(length(complex_enzyme_rxn)>0){
      node_comp <- sapply(complex_enzyme_rxn,function(comp_node) {
        
        r_comp <- strsplit(comp_node,"[+]")[[1]]
        r_comp <- r_comp[r_comp %in% rownames(rxn_vals_mat)]
        if(length(r_comp)>1){
          calc_comp <- apply(rxn_vals_mat[which(rownames(rxn_vals_mat) %in% r_comp),],2,min)
        }else if(length(r_comp)==1){
          calc_comp <- rxn_vals_mat[which(rownames(rxn_vals_mat) %in% r_comp),]
        }else{ calc_comp <- NULL }
        return(calc_comp)    
      })
      
      # if(class(node_comp)=="list"){
      if(inherits(node_comp, "list")){
        node_comp <- node_comp[sapply(node_comp,function(x) !is.null(x))]
        node_comp  <- do.call("rbind",node_comp)
      }else{ node_comp <- t(node_comp) }
      
      #isoenzymes_rxn
      isoenzymes_rxn <- node[grep(invert = T,"[+]",node)]
      
      if(length(isoenzymes_rxn)>0){
        isoenzymes_rxn <- isoenzymes_rxn[isoenzymes_rxn %in% rownames(rxn_vals_mat)]
        node_iso <- rxn_vals_mat[which(rownames(rxn_vals_mat) %in% isoenzymes_rxn),]
        # if(class(node_iso)=="numeric") { node_iso <- t(as.matrix(node_iso)) }
        if(inherits(node_iso, "numeric")) { node_iso <- t(as.matrix(node_iso)) }
        if(nrow(node_iso)==0){ node_iso<-NULL }
      }else{node_iso<-NULL}
      
      comp_iso_rxn<-rbind(node_comp,node_iso)
      if(!is.null(comp_iso_rxn)){
        calculated_node_val <- apply(comp_iso_rxn,2,max)
      }else{calculated_node_val<-default_value}
      
    }else{
      node <- node[which(node %in% rownames(rxn_vals_mat))]
      idx_iso_node<-which(rownames(rxn_vals_mat) %in% node)
      
      if(length(idx_iso_node)>1){
        node_iso <- rxn_vals_mat[idx_iso_node,]
        calculated_node_val <- apply(node_iso,2,max) 
      }else if(length(idx_iso_node)==1){
        calculated_node_val <- rxn_vals_mat[idx_iso_node,]
      }else{ calculated_node_val<-default_value }
    }
    return(calculated_node_val) 
  } )
  
  ## class() gives following error in R.4.0 and later versions: Error in if (class(y1) == "matrix") { : the condition has length > 1
  # if(class(y1)=="matrix"){
  #   node.vals <- t(y1)
  # }else if(class(y1)=="list"){ 
  if(inherits(y1, "matrix")){
    node.vals <- t(y1)
  }else if(inherits(y1,"list")){
    node.vals <- do.call(what="rbind",y1)
  }else{ node.vals <- mat.or.vec(nr = length(y1) ,nc = ncol(rxn_vals_mat))
  colnames(node.vals) <- colnames(rxn_vals_mat)
  rownames(node.vals) <- names(y1)
  node.vals[,] <- default_value
  message(sprintf("Warning for %s: all the reaction node values are %s",graphobject$module_name,default_value)) 
  cat("\n\n")
  }
  
  if(length(is.na(node.vals[,1]))>0){
    node.vals[is.na(node.vals[,1]),] <- default_value
  }
  return(node.vals)
}


get.All.module.genes <- function(hsa_module_data){
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  all_module_genes <-sapply(hsa_module_data, function(x) {
    module_elements<-unique(as.vector(x$graphobject$SIF[,c(1,3)]))
    rxns<-as.vector(unlist(sapply(module_elements[grep("R",module_elements)],function(x) {strsplit(x,split = ",")})))
    rxns<-as.vector(unlist(sapply(rxns,function(x) {strsplit(x,split = "[+]")})))
    Reax2Entrez <- x$KEGG_met_path_node$Reax2Entrez
    to_split <- grep(" R",Reax2Entrez$ReaxID)
    
    if(length(to_split)>0){
      for(ts in to_split){
        idx_push <- nrow(Reax2Entrez) + 1
        rxn_to_split <- Reax2Entrez$ReaxID[ts]
        item_to_push <- strsplit(rxn_to_split," ")[[1]]
        Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),1] <- rep(Reax2Entrez$EntrezIDs[ts],length(item_to_push)) 
        Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),2] <- item_to_push
      }
      Reax2Entrez <- Reax2Entrez[-to_split,]
    }
    module_genes<-Reax2Entrez$EntrezID[which(Reax2Entrez$ReaxID %in% rxns)]
    return(module_genes)
  }
  )
  all_module_genes<-unique(unlist(all_module_genes))
  all_module_genes_vec<-c()
  for(amg in all_module_genes){  all_module_genes_vec<-c(all_module_genes_vec, strsplit(amg,";")[[1]])}
  all_module_genes_vec<-unique(all_module_genes_vec)
  return(all_module_genes_vec)
}


get.All.module.rxns <- function(hsa_module_data){
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  all_module_rxn_vec <- sapply(hsa_module_data, function(x) {
    module_elements<-unique(as.vector(x$graphobject$SIF[,c(1,3)]))
    rxns<-as.vector(unlist(sapply(module_elements[grep("R",module_elements)],function(x) {strsplit(x,split = ",")})))
    rxns<-as.vector(unlist(sapply(rxns,function(x) {strsplit(x,split = "[+]")})))
    return(rxns)
  }
  )
  all_module_rxn_vec<-unique(unlist(all_module_rxn_vec))
  return(all_module_rxn_vec)
}


get.rxn.gene.matrix <- function(hsa_module_data){
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  rxn_gene_mat<-mat.or.vec(nc = 2, nr = 1)
  colnames(rxn_gene_mat) <- c("EntrezIDs","ReaxID")
  for(mat in 1:length(hsa_module_data)){
    rxn_gene_mat <- rbind(rxn_gene_mat,hsa_module_data[[mat]]$KEGG_met_path_node$Reax2Entrez)
  }
  rxn_gene_mat<-rxn_gene_mat[-1,]
  rxn_gene_mat<-unique(rxn_gene_mat)
  Reax2Entrez <- rxn_gene_mat
  to_split <- grep(" R",Reax2Entrez$ReaxID)
  
  if(length(to_split)>0){
    for(ts in to_split){
      idx_push <- nrow(Reax2Entrez) + 1
      rxn_to_split <- Reax2Entrez$ReaxID[ts]
      item_to_push <- strsplit(rxn_to_split," ")[[1]]
      Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),1] <- rep(Reax2Entrez$EntrezIDs[ts],length(item_to_push)) 
      Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),2] <- item_to_push
    }
    Reax2Entrez <- Reax2Entrez[-to_split,]
  }
  return(Reax2Entrez)  
}

get.gene.exp.of.module.genes <- function(combat.vals, all_module_genes_vec, min.exp=0.001, onesample=F){
  # combat.vals is normalized and batch corrected RNASeq (gene expression) data 
  # all_module_genes_vec: output of "get.All.module.genes"
  # min.exp lets us propagate flux/signal. This is in case RXN = 1 gene = 0.001
  idx_sellect_genes <- match(all_module_genes_vec,rownames(combat.vals))
  idx_NA_genes <- which(is.na(idx_sellect_genes))
  if(length(idx_NA_genes)>0){
    idx_sellect_genes <- idx_sellect_genes[-idx_NA_genes]
  }
  gene_exp_mat<-as.matrix(combat.vals[idx_sellect_genes,])
  # combat returns negative values which can be replaced by 0
  gene_exp_mat[gene_exp_mat<0]<-0
  if(onesample){
    gene_exp_mat <- (gene_exp_mat-min(gene_exp_mat,na.rm=T))/(max(gene_exp_mat,na.rm=T)-min(gene_exp_mat,na.rm=T)) 
  }
  else{
    gene_exp_mat <- t(apply(gene_exp_mat, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))   
  }
  gene_exp_mat[is.na(gene_exp_mat)]<-0
  gene_exp_mat[gene_exp_mat==0] <- min.exp
  return(gene_exp_mat)
}


get.RXNvals.from.genes <- function(all_module_rxn_vec,gene_exp_mat,rxn_gene_mat){
  # all_module_rxn_vec: output of "get.All.module.rxns"
  # gene_exp_mat: output of "get.gene.exp.of.module.genes"
  rxn_vals_mat <- as.matrix(mat.or.vec(nr = length(all_module_rxn_vec) ,nc = ncol(gene_exp_mat)))
  colnames(rxn_vals_mat) <- colnames(gene_exp_mat)
  rownames(rxn_vals_mat) <- all_module_rxn_vec
  
  for(r in 1:length(all_module_rxn_vec)){
    # use grep to deal with M00015: R01251, hsa00330: R01251 R01248
    #r_genes_1 <- rxn_gene_mat$EntrezIDs[which(rxn_gene_mat$ReaxID %in% all_module_rxn_vec[r])]
    r_genes_1 <- rxn_gene_mat$EntrezIDs[grep(all_module_rxn_vec[r],rxn_gene_mat$ReaxID)]
    
    # Added because there are some module reactions which do not have gene info e.g. R02189 in M00001 catalyzed by EC 2.7.1.63 is not exist for homosapiens
    if(length(r_genes_1)>0){
      
      r_genes_2 <- unique(unlist(sapply(r_genes_1,strsplit,";")))
      r_genes_3 <- r_genes_2[r_genes_2 %in% rownames(gene_exp_mat)]
      r_genes_3_vals <- gene_exp_mat[which(rownames(gene_exp_mat) %in% r_genes_3),]
      
      if(is.vector(r_genes_3_vals)){
        r_genes_3_vals <- r_genes_3_vals
      }else if(is.matrix(r_genes_3_vals)){
        r_genes_3_vals <- t(apply(r_genes_3_vals,2, function(x){quantile(x,0.9)}))
      }else{ r_genes_3_vals<-NA }
    }else{
      r_genes_3_vals<-NA
    }
    
    if(1==0){ 
      if(any(is.na(r_genes_3_vals))){cat("Rxn: ",all_module_rxn_vec[r]," has NA \n")}
    }
    rxn_vals_mat[r,] <- r_genes_3_vals
    
    # if((r %% 50)==0){
    #   cat(r,"...DONE\n")
    # }
  }
  cat("Module activity calculations are done successfully \n")    
  rxns_missing_val <- which(apply(rxn_vals_mat,1,function(x) all(is.na(x))))
  # there is no gene expression info for these reactions
  cat(paste0("number of the reactions with missing values: ", length(rxns_missing_val),"\n"))
  rxn_vals_mat <- rxn_vals_mat[-rxns_missing_val,]
  return(rxn_vals_mat)
}


get.metabolite.list <- function(hsa_module_data){
  all_metabolites <- sort(unique(unlist(sapply(hsa_module_data, function(x) {
    nodes <- c(x$graphobject$SIF[,1],x$graphobject$SIF[,3])
    metabolites <- nodes[c(grep("C",nodes),grep("G",nodes))]
    return(metabolites)
  }))))
  return(all_metabolites)
}



