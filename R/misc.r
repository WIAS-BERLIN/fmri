Varcor<-function(lkern,h,d=1){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h)+1,trunc(h)+1)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
2*sum(penl)^2/sum(diff(penl)^2)
}
Varcor.gauss<-function(h,interv = 1){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
h<-h/2.3548*interv
ih<-trunc(4*h)+1
dx<-2*ih+1
d<-length(h)
penl <- dnorm(((-ih[1]):ih[1])/h[1])
if(d==2) penl <- outer(penl,dnorm(((-ih[2]):ih[2])/h[2]),"*")
if(d==3) penl <- outer(penl,outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
2*sum(penl)^2/sum(diff(penl,interv)^2)/interv^(d)
}

SpatialCorr<-function(lkern,h,d=1){
#
#   Calculates the correlation of 
#
#   colored noise that was produced by smoothing with lkern and bandwidth h
#
#  !!! the result is not monotone in h for   lkern="Triangle" (all d) and lkern="Uniform" (d>1)
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h+1),trunc(h)+1)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-as.vector(switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x)))
dim(penl)<-rep(dx,d)
if(d==1) z<-sum(penl[-1]*penl[-dx])/sum(penl^2)
if(d==2) z<-sum(penl[-1,]*penl[-dx,])/sum(penl^2)
if(d==3) z<-sum(penl[-1,,]*penl[-dx,,])/sum(penl^2)
z
}

SpatialCorr.gauss<-function(h,interv=1){
#
#   Calculates the correlation of 
#
#   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
#
#   Result does not depend on d for "Gaussian" kernel !!
#
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
h<-h/2.3548*interv
ih<-trunc(4*h+1)
dx<-2*ih+1
penl<-dnorm(((-ih):ih)/h)
sum(penl[-(1:interv)] * penl[-((dx-interv+1):dx)])/sum(penl^2)
}

Spatialvar<-function(lkern,lkern0,h,h0,d){
#
#   Calculates the factor of variance reduction obtained at bandwidth h in 
#
#   case of colored noise that was produced by smoothing with lkern0 and bandwidth h0
#
#   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h),trunc(h))
ih<-max(1,ih)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
dim(penl)<-rep(dx,d)
if(lkern0=="Gaussian") h0<-h0/2.3548
ih<-switch(lkern0,Gaussian=trunc(4*h0),trunc(h0))
ih<-max(1,ih)
dx0<-2*ih+1
x<- ((-ih):ih)/h0
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl0<-switch(lkern0,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
dim(penl0)<-rep(dx0,d)
penl0<-penl0/sum(penl0)
dz<-dx+dx0-1
z<-array(0,rep(dz,d))
if(d==1){
for(i1 in 1:dx0) {
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
z[-ind1]<-z[-ind1]+penl*penl0[i1]
}
} else if(d==2){
for(i1 in 1:dx0) for(i2 in 1:dx0){
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
ind2<-c(0:(i2-1),(dz-dx0+i2):dz+1)
ind2<-ind2[ind2<=dz][-1]
z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
}
} else if(d==3){
for(i1 in 1:dx0) for(i2 in 1:dx0) for(i3 in 1:dx0){
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
ind2<-c(0:(i2-1),(dz-dx0+i2):dz+1)
ind2<-ind2[ind2<=dz][-1]
ind3<-c(0:(i3-1),(dz-dx0+i3):dz+1)
ind3<-ind3[ind3<=dz][-1]
z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
}
}
sum(z^2)/sum(z)^2
}

Spatialvar.gauss<-function(h,h0,d,interv=1){
#
#   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
#
#   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
#
#   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
  h0 <- pmax(h0,1e-5)
  h<-h/2.3548*interv
if(length(h)==1) h<-rep(h,d)
ih<-trunc(4*h)
ih<-pmax(1,ih)
dx<-2*ih+1
penl<-dnorm(((-ih[1]):ih[1])/h[1])
if(d==2) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
if(d==3) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
dim(penl)<-dx
h0<-h0/2.3548*interv
if(length(h0)==1) h0<-rep(h0,d)
ih<-trunc(4*h0)
ih<-pmax(1,ih)
dx0<-2*ih+1
x<- ((-ih[1]):ih[1])/h0[1]
penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
dim(penl0)<-dx0
penl0<-penl0/sum(penl0)
dz<-dx+dx0-1
z<-array(0,dz)
if(d==1){
for(i1 in 1:dx0) {
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
z[-ind1]<-z[-ind1]+penl*penl0[i1]
}
} else if(d==2){
for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
ind1<-ind1[ind1<=dz[1]][-1]
ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
ind2<-ind2[ind2<=dz[2]][-1]
z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
}
} else if(d==3){
for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]) for(i3 in 1:dx0[3]){
ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
ind1<-ind1[ind1<=dz[1]][-1]
ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
ind2<-ind2[ind2<=dz[2]][-1]
ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
ind3<-ind3[ind3<=dz[3]][-1]
z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
}
}
sum(z^2)/sum(z)^2*interv^d
}

geth.gauss<-function(corr,step=1.002,interv=1){
#   get the   bandwidth for lkern corresponding to a given correlation
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
  if (corr < 0.1) {
    h <- 0
  } else { 
    h <- .5
    z <- 0
    while (z<corr) {
      h <- h*step
      z <- get.corr.gauss(h,interv)
    }
    h <- h/step
  }
  h
}

get3Dh.gauss<-function(vred,h0,vwghts,step=1.002,interv=1){
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
  h0 <- pmax(h0,1e-5)
  n<-length(vred)
vred1<-vred
h<-.5/vwghts
fixed<-rep(FALSE,length(vred))
while(any(!fixed)){
ind<-(1:n)[!fixed][vred[!fixed]>=Spatialvar.gauss(h,1e-5,3,interv)]
vred1[ind]<-Spatialvar.gauss(h,h0,3,interv)
fixed[ind]<-TRUE
h<-h*step
}
hvred<-matrix(0,3,n)
hh<-.01/vwghts
h<-h0
fixed<-rep(FALSE,length(vred))
while(any(!fixed)){
ind<-(1:n)[!fixed][vred1[!fixed]>=Spatialvar.gauss(h,1e-5,3,interv)]
hvred[,ind]<-h
fixed[ind]<-TRUE
hh<-hh*step
h<-sqrt(h0^2+hh^2)
}
t(hvred)
}

