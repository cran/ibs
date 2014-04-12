/* based on the file from the folder seimonex */
#include<R.h>
#include<Rinternals.h>
#include "bsp.h"
double gdiv(double a,double b){if(a==0.0&&b==0.0)return(0.0);else return(a/b);}
double bsp(int i, int ord, double x, int nk, double kns[]){
  if(i<0||i>nk-ord-1){warning("illegal i value: i=%d;nk-ord=%d-%d=%d\n",i,nk,ord,nk-ord);return R_NaN;}
  if(x<kns[i]||x>kns[i+ord])return(0.0);
  int k=nk-1; while(kns[k]==kns[k-1])k--; k--; /* printf("n=%d\n",k); */
  if(ord==1){
    if(i!=k){
      /* if(i==5)printf("kns[%d]=%f,kns[%d]=%f,x=%f\n",i,kns[i],i+1,kns[i+1],x); */
      return((kns[i]<=x && x<kns[i+1])? 1.0 : 0.0);
    }else{
      if(i==k){
	return((kns[i]<=x && x<=kns[i+1])? 1.0 : 0.0);
      }else return R_NaN;
    }
  }else
    return(gdiv((x-kns[i])*bsp(i,ord-1,x,nk,kns), kns[i+ord-1]-kns[i])+
	   gdiv((kns[i+ord]-x)*bsp(i+1,ord-1,x,nk,kns),kns[i+ord]-kns[i+1])
	   );
}
SEXP bsp0(SEXP i, SEXP ord, SEXP x, SEXP nk, SEXP kns){
  SEXP ans;
  PROTECT(ans=allocVector(REALSXP,1));
  REAL(ans)[0]=bsp(INTEGER(i)[0],INTEGER(ord)[0],
		   REAL(x)[0],INTEGER(nk)[0],REAL(kns));
  UNPROTECT(1);
  return ans;
}
SEXP bspbases(SEXP xs, SEXP kns,SEXP order){
  SEXP ans;
  int nx=length(xs), nk=length(kns),i,j;
  int ord=INTEGER(order)[0], nb=nk-ord;
  PROTECT(ans=allocVector(REALSXP,nx*nb));
  for(i=0;i<nx;i++){
    for(j=0;j<nb;j++){
      REAL(ans)[j+i*nb]=bsp(j,ord,REAL(xs)[i],nk,REAL(kns));
      /* printf("%1.2f\t",REAL(ans)[j+i*nb]); */
    }
    /* printf("\n"); */
  }
  UNPROTECT(1);
  /* printf("%1.2f\n",bsp(10,1,1.0,12,REAL(kns))); */
  return(ans);
}
SEXP bsp1(SEXP xs, SEXP order, SEXP kns, SEXP coef){
  SEXP ans;
  int nx=length(xs), nk=length(kns), i, j;
  int ord=INTEGER(order)[0];
  PROTECT(ans=allocVector(REALSXP,nx));
  for(i=0;i<nx;i++)REAL(ans)[i]=0.0;
  for(i=0;i<nx;i++)
    for(j=0; j<nk-ord; j++)
      REAL(ans)[i] += REAL(coef)[j]*bsp(j,ord,REAL(xs)[i],nk,REAL(kns));
  UNPROTECT(1);
  return(ans);
}
/* integral of a B-spline with coefficient vector coef and knot vector
   nks on the interval [min(kns),x] with x in the range of the knots,
   using the formula (33) on page 128 of de Boor (2001) A pratical
   guide to splines (revised edition). We do this with the special
   case of the smallest and largest knots have equal multiplicity, but
   the formula might be true when there is other forms of multiplicity
   or no multiplicity among the knots.
*/
/* It was later discovered that the relationship is true only on the
   interval from the first to second last knot. So it is
   necessary to pad the original knots sequence by an (arbitrary)
   extra value, e.g. the last knot.*/  
SEXP Ibsp1(SEXP xs, SEXP order, SEXP kns, SEXP coef){
  SEXP ans,kns_ext,coef_ext;
  int nx=length(xs), nk=length(kns), i, j, k;
  int ord=INTEGER(order)[0];
  double tmp;
  PROTECT(ans=allocVector(REALSXP,nx));
  PROTECT(kns_ext=allocVector(REALSXP,nk+1));
  PROTECT(coef_ext=allocVector(REALSXP,nk+1));
  for(i=0;i<nk;i++){
    REAL(kns_ext)[i]=REAL(kns)[i];
    REAL(coef_ext)[i]=REAL(coef)[i];
  }
  REAL(kns_ext)[nk]=REAL(kns)[nk-1];
  REAL(coef_ext)[nk]=0.0;
  
  for(k=0;k<nx;k++){
    REAL(ans)[k]=0.0;
    for(i=0;i<nk-ord;i++){
      tmp=0.0;
      for(j=0; j<=i; j++){
	tmp+=REAL(coef_ext)[j]*(REAL(kns_ext)[j+ord]-REAL(kns_ext)[j])/ord;
      }
      REAL(ans)[k]+=tmp* bsp(i,ord+1,REAL(xs)[k],nk+1,REAL(kns_ext));
    }
  }
  UNPROTECT(3);
  return ans;	
}
