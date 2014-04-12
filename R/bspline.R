bspline <-
function(x,knots,ord,coef=rep(1,length(knots)-ord)){
  ord <- as.integer(ord)
  if(length(coef)!=length(knots)-ord)stop("length(knots)-ord!=length(coef)!")
  if(length(as.double(x))==0)return(numeric(0))
  .Call("bsp1",as.double(x),as.integer(ord),as.double(sort(knots)),as.double(coef))  
}
