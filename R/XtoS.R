XtoS <- function(X) {
  id.rows <- apply(X, 1, paste, sep = "", collapse = ".")
  id.rows
  Xtab <- table(id.rows)
  Xtab
  
  Xcat <- names(Xtab)
  Xcat
  
  S <- matrix(as.numeric(unlist(
    strsplit(Xcat, ".", fixed = TRUE))), 
    length(Xtab), ncol(X), byrow = TRUE)
  
  counts <- Xtab
  names(counts) <- NULL
  counts <- as.vector(counts)
  
  if(is.null(rownames(X))){
    rownames(X) <- paste("id", 1:nrow(X), sep = "_")
  }
  
  rownames(S) <- NULL
  
  out <- list(S = S, counts = counts, 
    idX = id.rows, idS = Xcat, Xnames = rownames(X))
  out
}
