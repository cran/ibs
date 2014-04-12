ibs <-
function(x,knots,ord,coef=rep(1,length(knots)-ord)){
 if(length(coef)!=length(knots)-ord)stop("length(knots)-ord!=length(coef)!")
 if(length(as.double(x))==0)return(numeric(0));
 x[x>max(knots)] <- max(knots)
  .Call("Ibsp1",as.double(x),as.integer(ord),
        as.double(c(sort(knots),rep(max(knots),ord))),
        as.double(coef))  
}
