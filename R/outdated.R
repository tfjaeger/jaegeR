myCenter <- function(x) {
  if (is.numeric(x)) { return(x - mean(x, na.rm=T)) }
  if (is.factor(x)) {
    x <- as.numeric(x)
    return(x - mean(x, na.rm=T))
  }
  if (is.data.frame(x) || is.matrix(x)) {
    m <- matrix(nrow=nrow(x), ncol=ncol(x))
    colnames(m) <- paste("c", colnames(x), sep="")
    for (i in 1:ncol(x)) {
      m[,i] <- x[,i] - mean(x[,i], na.rm=T)
    }
    return(as.data.frame(m))
  }
}

myScale <- function(x) {
  if (is.numeric(x)) { return(as.numeric(scale(x))) }
  if (is.factor(x)) {
    x <- as.numeric(as.character(x))
    return(as.numeric(scale(x)))
  }
  if (is.data.frame(x) || is.matrix(x)) {
    m <- matrix(nrow=nrow(x), ncol=ncol(x))
    colnames(m) <- paste("s", colnames(x), sep="")
    for (i in 1:ncol(x)) {
      m[,i] <- as.numeric(scale(x[,i]))
    }
    return(as.data.frame(m))
  }
}

myFactorCleanup <- function(x) {
  if (is.factor(x)) {
    if (is.ordered(x)) { return(x) }
    else { return(as.factor(as.character(x))) }
  }
  if (is.data.frame(x) || is.matrix(x)) {
    for (i in 1:ncol(x)) {
      if(is.factor(x[,i]) & !is.ordered(x[,i])) {
        x[,i] <- as.factor(as.character(x[,i]))
      }
    }
    return(as.data.frame(x))
  }
}
