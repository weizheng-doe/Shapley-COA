###############################################################
######################code for Algorithm 1#####################
###############################################################

# Note:
# This code only describes the framework of Algorithm 1, 
# some of which need to be predefined, 
# such as the value function (val), 
# and the primitive polynomial (f_d) when the COA method is used with a prime power d.

####=================1. Code for the LS method==========================
##--------predefined function: value function------------
# Input: a set; Output: value
# val<-function(set){}

##--------Subfunction: generate an LS-----------
# Input: d; Output: an LS
onels<-function(d){
  sq <- matrix(nrow=d, ncol=d) 
  sq[1,]<-c(1,sample(d-1)+1);
  for (r in 2:d) {
    sq[r,]<-c(sq[r-1,-1],sq[r-1,1])
  }
  pcol<-sample(d);prow<-c(1,sample(c(2:d)));
  sq<-sq[,pcol]; sq<-sq[prow,];
  return(sq)
}

##--------Main function-----------
# Input: 
# d: the number of players; n: the sample size; val: the predefined value function
# Output:
# sh: a vector including Shapley values of all players.
est.shls<-function(d,n,val){
  sh<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i])-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}



####=================2. Code for the COA method with a prime d==========================
##--------predefined function: value function------------
# Input: a set; Output: value
# val<-function(set){}

##--------Subfunction: generate a COA with a prime d-----------
# Input: a prime d; Output: a COA
onecoa<-function(d){
  firstline<-c(0,sample(d-1)); cz<-matrix(0,d-1,d); D<-matrix(0,d*(d-1),d)
  for(i in 1:(d-1)){
    cz[i,]<-(firstline*i)%%d
  }
  for(j in 0:(d-1)){
    D[(j*(d-1)+1):((j+1)*(d-1)),]<-(cz+j)%%d
  }
  D<-D+1
  return(D)
}

##--------Main function-----------
# Input: 
# d: the number of players; n: the sample size; val: the predefined value function
# Output:
# sh: a vector including Shapley values of all players.
est.shcoa<-function(d,n,val){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i])-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}



####=================3. Code for the COA method with a prime power d==========================
library(gtools) #permutations
##--------predefined function: value function------------
# Input: a set; Output: value
# val<-function(set){}

##--------predefined values: primitive polynomial(f_d)-----------
# For example:
# f_4<-c(1,1,1)
# f_8<-c(1,0,1,1)
# f_9<-c(1,1,2)
# f_16<-c(1,0,0,1,1)
# f_25<-c(1,1,2)
# f_27<-c(1,0,2,1)
# f_32<-c(1,0,0,1,0,1)
# f_49<-c(1,1,3)
# f_64<-c(1,0,0,0,0,1,1)
# f_81<-c(1,0,0,1,2)
# f_125<-c(1,0,3,2)
# f_128<-c(1,0,0,0,0,0,1,1)
# f_243<-c(1,0,0,0,2,1)
# f_256<-c(1,0,0,0,1,1,1,0,1)
# f_625<-c(1,0,1,2,2)

##--------Subfunction1: Polynomial Division-----------
# Input:
# F1 is the dividend, and f2 is the divisor
#For example, if the dividend is x^5+2x^3+1 and the divisor is x^4+2, then f1=c(1,0,2,0,0,1),f2=c(1,0,0,0,2).

poly.div<-function(f1,f2){ 
  d1<-length(f1)-1; d2<-length(f2)-1; 
  if(d1<d2){f1<-c(rep(0,d2-d1),f1); d1<-d2}
  dd<-d1-d2;
  divisend<-f1; divisor<-f2;
  b<-rep(0,dd+1) 
  for(i in 1:(dd+1)){
    b[i]<-divisend[1]/divisor[1]
    for(j in 1:(d2+1)){
      divisend[j]<-divisend[j]-b[i]*divisor[j]
    }
    divisend<-divisend[-1]
  }
  r<-divisend
  return(r)
}

gfpoly.div<-function(f1,f2,s){
  r<-poly.div(f1,f2)
  gfr<-r%%s
  return(gfr)
}

#----------------Subfunction2: polynomial mutiplication of GF(s) with a prime s-------------
gfpoly.multi<-function(f1,f2,s){ #f1,f2 are multipliers
  d1<-length(f1);d2<-length(f2);hd<-d1+d2-1;
  t<-matrix(0,d1,hd)
  for(i in 1:d1){ # each coefficient of f1 multiply d_2
    t[i,i:(i+d2-1)]<-f1[i]*f2 #the second coefficient is one place behind the first
  }
  multi<-apply(t,2,sum)
  gft<-multi%%s #GF(s)
  return(gft)
}

#---------------Subfunction3: polynomial additive of GF(s) with a prime s-------------------
gfpoly.add<-function(f1,f2,s){
  l1<-length(f1); l2<-length(f2);
  if(l1<l2){f1<-c(rep(0,l2-l1),f1);}
  else if(l1>l2){f2<-c(rep(0,l1-l2),f2);}
  f<-(f1+f2)%%s
  return(f)
}


##--------Subfunction4: generate a COA with a prime d-----------
# Input:
# d: number of component
# p: a prime satisfying that d is a power of p
# f_d: primative polynomial on GF(d)
onecoa<-function(d,p,f_d){
  firstr<-c(0,sample(d-1));
  r<-round(log(d,p))
  #multiplication table and addition table on GF(d)
  if(r==1){
    M_d<-matrix(0,p,p); A_d<-matrix(0,p,p);
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        M_d[i,j]<-((i-1)*(j-1))%%p; M_d[j,i]<-M_d[i,j];
        A_d[i,j]<-(i+j-2)%%p; A_d[j,i]<-A_d[i,j];
      }
      M_d[i,i]<-(i-1)^2%%p; A_d[i,i]<-(2*(i-1))%%p;
    }
    M_d[p,p]<-(p-1)^2%%p; A_d[p,p]<-(2*(p-1))%%p;
  }
  else{
    entry<-permutations(p,r,0:(p-1),repeats=TRUE); cases<-permutations(d,2,repeats= TRUE);
    M1<-matrix(nrow=d,ncol=d*r)  
    for(i in 1:(d^2)){ 
      M1[cases[i,1],((cases[i,2]-1)*r+1):(cases[i,2]*r)]<-gfpoly.div(gfpoly.multi(entry[cases[i,1],],entry[cases[i,2],],p),f_d,p)
    }
    A1<-matrix(nrow=d,ncol=d*r)
    for(i in 1:(d^2)){ 
      A1[cases[i,1],((cases[i,2]-1)*r+1):(cases[i,2]*r)]<-gfpoly.add(entry[cases[i,1],],entry[cases[i,2],],p)
    }
    M_d<-matrix(nrow=d,ncol=d); A_d<-matrix(nrow=d,ncol=d)
    for(i in 1:d){
      for(j in 1:d){
        I<-0
        for(l in 1:d){
          if(round(sum(abs(M1[i,((j-1)*r+1):(j*r)]-entry[l,])))==0) {M_d[i,j]<-l-1; I<-I+1;}
          if(round(sum(abs(A1[i,((j-1)*r+1):(j*r)]-entry[l,])))==0) {A_d[i,j]<-l-1; I<-I+1;}
          if(I==2) break;
        }
      }
    }
  }
  #construction
  nr<-d*(d-1); coa<-matrix(nrow=nr,ncol=d);
  for(i in 1:d){
    for(j in 1:(d-1)){
      for(k in 1:d)
        coa[(i-1)*(d-1)+j,k]<-A_d[i,M_d[j+1,firstr[k]+1]+1]
    }
  }
  return(coa+1)
}


##--------Main function-----------
# Input: 
# d: the number of players; n: the sample size; val: the predefined value function
# p: a prime satisfying that d is a power of p
# f_d: primative polynomial on GF(d)

# Output:
# sh: a vector including Shapley values of all players.

est.shcoa<-function(d,p,f_d,n,val){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d,p,f_d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i])-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}