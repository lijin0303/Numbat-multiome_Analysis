dirC <- function(p){if(!dir.exists(p)){dir.create(p,recursive = T)}}
unilen <- function(x){
  return(length(unique(x)))
}
cfct <- function(x,rev=T){
  clusterno <- unilen(x)
  if(rev){
    fx = factor(x,levels = rev(paste0("C",seq_len(clusterno))))
  }else{
    fx = factor(x,levels = paste0("C",seq_len(clusterno))) 
  }
  return(fx)
}
winsor <- function(x,q=0.95){
  cutoff <- quantile(x,q)
  x[x>cutoff] <- cutoff
  return(x)
}
minmax <- function(s){
  normed <- (s-min(s))/(max(s)-min(s))
  return(normed)
}
gridlist <- function(gridL,f){
  pdf(f)
  for(i in seq_along(gridL)){
    grid::grid.draw(gridL[[i]])
    if(i != length(gridL)){grid::grid.newpage()}}
  dev.off()
}
paletteVec <- function(v,p="Zissou1"){
  ele <- levels(v)
  hexcode <- as.character(wes_palette(p))
  psub <- hexcode[1:length(ele)]
  names(psub) <- ele
  return(psub)
}
list2flat <- function(l){
  ll = unlist(lapply(seq_along(l),function(i) rep(names(l)[i],
                                                  each=lengths(l)[i])))
  l2df <-  data.frame(values = unlist(l),
                      names = ll)
  rownames(l2df) <- NULL
  return(l2df)
} 
matrixFill <- function(mat,rowOrder,colOrder,fillin){
  newmat <- matrix(fillin,nrow =length(rowOrder),ncol=length(colOrder))
  colnames(newmat) <- colOrder; rownames(newmat) <- rowOrder
  newmat[rownames(mat),colnames(mat)] <- mat
  return(newmat)
}
matmap <- function(mat,from,to){
  omat <- mat
  if(class(mat)=="data.frame"){omat = as.matrix(omat)}
  
  for(i in seq_along(from)){
    omat[omat==from[i]] <- to[i]
  }
  if (class(to)=="numeric"){
    omat <- matrix(as.numeric(omat), ncol = ncol(omat))
    rownames(omat) = rownames(mat);colnames(omat) = colnames(mat)
  }
  return(omat)
}
nature_supp <- function(id,fno,Opath){
  id2 = rev(rev(strsplit(id,"-")[[1]])[-1])
  id2[2] <- paste0("2",id2[2])
  id2[3] <- gsub("^0+","",id2[3])
  id2 = paste0(id2,collapse = "_")
  Fpath <- glue("https://static-content.springer.com/esm/",
       "art%3A10.1038%2",
       "Fs{id}/MediaObjects/",
       "{id2}_MOESM{fno}_ESM.xlsx")
  message(glue("wget {Fpath} -O {Opath}"))
  return(glue("wget {Fpath} -O {Opath}"))
}
oneCol <- function(obj,f){
  write.table(as.vector(obj),f,col.names = F,row.names = F,quote=F)
}
mgrep <- function(pVec,str){
  return(any(map_lgl(pVec,\(x)grepl(x,str))))
}
nmap <- function(l){
  return(lapply(seq_along(l),\(i) l[[i]] |> mutate(cluster = names(l)[i])))
}
tablem <- function(s,min=1){
  tab <- table(s)
  return(as.vector(sort(tab[tab>min],decreasing = T)))
}
displayCol <- function(hexCol){
  colDF <- as.data.frame(hexCol) |> 
    rownames_to_column("description") |> 
    mutate(id = 1:n())
  colg <- ggplot(colDF, aes(x = 1, y = max(id) - id, fill = description)) +
    geom_tile() +
    geom_text(aes(label = description), size = rel(5)) +
    scale_fill_manual(values = hexCol) +
    theme_void() +
    theme(legend.position = "none")+
    coord_equal(ratio = 0.15)
  return(colg)
}
list2upset <- function(l){
  combvec <-  Reduce("union",l)
  indMat <- bind_cols(map(l,\(x) combvec %in% x))
  indMat[is.na(indMat)] <- F
  indDF <- as.data.frame(1*indMat)
  return(indDF)
}
TODsave <- function(F){
  TOD <- Sys.Date()
  F <- gsub("\\.",paste0("_",TOD,"."),F)
  return(F)
}
matNorm <- function(cM,nby="row"){
  cM <- as.matrix(cM)
  if(nby=="row"){
    cM <-  sweep(cM,1,rowSums(cM),`/`)
  }else{
    cM <-  sweep(cM,2,colSums(cM),`/`)
  }
  return(cM)
}
grpmean <- function(grp,mat){
  if(ncol(mat)==length(grp)){
    mat <- t(mat)
  }else{
    stopifnot(nrow(mat)!=length(grp))
  }
  mat = mat[names(grp),]
  grpsum = rowsum(mat,grp,na.rm=T)
  sorted_grpsize = table(grp)[rownames(grpsum)]
  meangrp = grpsum/as.numeric(sorted_grpsize)
  return(t(meangrp))
}