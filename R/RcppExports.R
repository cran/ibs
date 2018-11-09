bsbases <-
function(x,knots,ord){
  ord <- as.integer(ord)
  if(length(knots)<=ord)stop("length of knots <= ord!\n")
  if(length(as.double(x))==0){
    return(matrix(nrow=0,ncol=length(knots)-ord))
  }
  tmp <- .Call("_ibs_bsbasesCpp",as.double(x),as.double(knots),as.integer(ord),PACKAGE="ibs")
  matrix(tmp,nrow=length(x),byrow=TRUE)
}
bspline <-
function(x,knots,ord=4,coef=rep(1,length(knots)-ord)){
  ord <- as.integer(ord)
  if(length(coef)!=length(knots)-ord)stop("length(knots)-ord!=length(coef)!")
  if(length(as.double(x))==0)return(numeric(0))
  .Call("_ibs_bsplineCpp",as.double(x),as.integer(ord),as.double(sort(knots)),as.double(coef),PACKAGE="ibs")  
}
ibs <-
function(x,knots,ord=4,coef=rep(1,length(knots)-ord)){
 if(length(coef)!=length(knots)-ord)stop("length(knots)-ord!=length(coef)!")
 if(length(as.double(x))==0)return(numeric(0));
 knots <- sort(knots);
 if(any(x<knots[1] | x>knots[length(knots)-ord+1]))
     stop("Some x value(s) are out of the range from the smallest to the ord-th largest knots!\n")
  .Call("_ibs_ibsCpp",as.double(x),as.integer(ord),
        as.double(knots),
        as.double(coef),PACKAGE="ibs")  
}
