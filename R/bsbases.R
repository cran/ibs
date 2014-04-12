bsbases <-
function(x,knots,ord){
  ord <- as.integer(ord)
  if(length(knots)<=ord)stop("length of knots <= ord")
  if(length(as.double(x))==0){
    return(matrix(nrow=0,ncol=length(knots)-ord))
  }
  tmp <- .Call("bspbases",as.double(x),as.double(knots),as.integer(ord))
  matrix(tmp,nrow=length(x),byrow=TRUE)
}
