CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   extract element of 3D array
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function getlwght(lwght,dw1,dw2,dw3,j1,j2,j3)
      integer dw1,dw2,dw3,j1,j2,j3
      real*8 lwght(dw1,dw2,dw3)
      getlwght=lwght(j1,j2,j3)
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Plateau
C          Kern=2     Epanechnicov
C          Kern=3     Gaussian
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         lkern=exp(-xsq*8.d0)
      ELSE IF (kern.eq.4) THEN
         lkern=exp(-xsq*18.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute aws-weights  w_{ij}
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awswghts(n1,n2,n3,j1,j2,j3,dv0,thi,theta,
     1                    vwghts,skern,spf,spmin,spmax,bii,wj)
      implicit logical (a-z)
      integer n1,n2,n3,j1,j2,j3,dv0,skern
      real*8 thi(dv0),theta(n1,n2,n3,dv0),vwghts(dv0),spf,spmin,spmax,
     1       bii,wj,wjin
      integer k
      real*8 sij,z
      wjin=wj
      sij=0.d0
C  compute distance in sij
      DO k=1,dv0
         z=thi(k)-theta(j1,j2,j3,k)
         sij=sij+z*z*vwghts(k)
      END DO
      sij=bii*sij
      IF (sij.gt.spmax) THEN
         wj=0.d0
      ELSE IF (skern.eq.1) THEN
C  skern == "Plateau"
         wj=wj*min(1.d0,1.d0-spf*(sij-spmin))
      ELSE IF (skern.eq.2) THEN
C  skern == "Triangle"
         wj=wj*(1.d0-sij)
      ELSE
C  skern == "Exp"
         IF (sij.gt.spmin) wj=wj*exp(-spf*(sij-spmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute aws-weights  w_{ij}
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awswght3(n,jind,dv0,thi,theta,vwghts,skern,
     1                    spf,spmin,spmax,bii,wj)
      implicit logical (a-z)
      integer n,jind,dv0,skern
      real*8 thi(dv0),theta(dv0,n),vwghts(dv0),spf,spmin,spmax,
     1       bii,wj,wjin
      integer k
      real*8 sij,z
      wjin=wj
      sij=0.d0
C  compute distance in sij
      DO k=1,dv0
         z=thi(k)-theta(k,jind)
         sij=sij+z*z*vwghts(k)
      END DO
      sij=bii*sij
      IF (sij.gt.spmax) THEN
         wj=0.d0
      ELSE IF (skern.eq.1) THEN
C  skern == "Plateau"
         wj=wj*min(1.d0,1.d0-spf*(sij-spmin))
      ELSE IF (skern.eq.2) THEN
C  skern == "Triangle"
         wj=wj*(1.d0-sij)
      ELSE
C  skern == "Exp"
         IF (sij.gt.spmin) wj=wj*exp(-spf*(sij-spmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws2(y,si2,mask,wlse,n1,n2,n3,dv,dv0,hakt,lambda,
     1                  theta,ncores,bi,thn,kern,skern,spmin,spmax,
     2                  lwght,wght,vwghts,swjy,thi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,kern,skern,dv,dv0,ncores
      logical aws,wlse,mask(*)
      real*8 y(dv,*),theta(dv0,*),bi(*),thn(dv,*),lambda,spmax,
     1       wght(2),si2(*),hakt,lwght(*),spmin,vwghts(dv0),
     2       thi(dv0,ncores),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,iind,jind,thrednr
      real*8 bii,swj,swjy(dv,ncores),wj,hakt2,spf,si2j,si2i,swjv
      external getlwght
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,n3,kern,skern,dv,dv0,aws,wlse,mask,y,theta,bi,thn,
C$OMP& lambda,spmax,wght,si2,hakt,lwght,spmin,vwghts,thi,spf,
C$OMP& ncores,ih1,ih2,ih3,clw1,clw2,clw3,dlw1,dlw2,dlw3,n,hakt2,swjy)
C$OMP& PRIVATE(iind,i1,i2,i3,k,si2i,bii,swjv,swj,
C$OMP& thrednr,j1,j2,j3,jw1,jw2,jw3,jind,wj,si2j)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
!$         thrednr = omp_get_thread_num()+1
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(mask(iind)) THEN
            DO k=1,dv
               thn(k,iind)=0.d0
            END DO
            CYCLE
         END IF
         si2i=si2(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swjv=0.d0
         DO k=1,dv
            swjy(k,thrednr)=0.d0
         END DO
         DO k=1,dv0
            thi(k,thrednr)=theta(k,iind)
         END DO
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               DO jw1=1,dlw1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+n1*(j2-1)+n1*n2*(j3-1)
                  IF(mask(jind)) CYCLE
                  wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                  if(wj.le.0.d0) CYCLE
                  si2j=si2(jind)
                  IF (aws) THEN
                     call awswght3(n,jind,dv0,thi(1,thrednr),
     1                    theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                  END IF
                  if(wlse) THEN 
                     wj=wj*si2j
                  ELSE
                     swjv=swjv+wj/si2j
                  END IF
                  swj=swj+wj
                  call daxpy(dv,wj,y(1,jind),1,swjy(1,thrednr),1)
               END DO
            END DO
         END DO
         DO k=1,dv
            thn(k,iind)=swjy(k,thrednr)/swj
         END DO
         IF(wlse) THEN
            bi(iind)=swj
         ELSE
            bi(iind)=swj*swj/swjv
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,bi)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsv(y,res,si2,mask,wlse,n1,n2,n3,n4,dv,dv0,hakt,
     1                  lambda,theta,ncores,bi,resnew,thn,kern,skern,
     2                  spmin,spmax,lwght,wght,vwghts,swjy,thi,resi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,n4,kern,skern,dv,dv0,ncores
      logical aws,wlse,mask(*)
      real*8 res(n4,*),y(dv,*),theta(dv0,*),
     1       bi(*),thn(dv,*),lambda,spmax,wght(2),
     1       si2(*),hakt,lwght(*),spmin,vwghts(dv0),thi(*),
     1       resi(*),getlwght,resnew(n4,*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1       clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,thrednr,iind,jind,
     2       sthrednr,rthrednr
      real*8 bii,swj,swjy(*),wj,hakt2,spf,si2j,si2i,swjv,
     1       sresisq,resik
      external getlwght
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(kern,skern,ncores,mask,res,
C$OMP& y,theta,bi,thn,spmax,wght,si2,hakt,lwght,spmin,vwghts,
C$OMP& thi,resi,resnew,ih1,ih2,ih3,
C$OMP& swjy,hakt2,spf)
C$OMP& FIRSTPRIVATE(n,n1,n2,n3,n4,dv,dv0,lambda,clw1,clw2,clw3,
C$OMP& dlw1,dlw2,dlw3,aws,wlse)
C$OMP& PRIVATE(iind,jind,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,k,si2i,
C$OMP& bii,swj,swjv,wj,si2j,sresisq,thrednr,sthrednr,rthrednr,resik)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
!$         thrednr = omp_get_thread_num()+1
         sthrednr = (thrednr-1)*dv
         rthrednr = (thrednr-1)*n4
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(mask(iind)) THEN
            DO k=1,dv
               thn(k,iind)=0.d0
            END DO
            DO k=1,n4
               resnew(k,iind)=0.d0
            END DO
            bi(iind)=1.d0            
            CYCLE
         END IF
         si2i=si2(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swjv=0.d0
         DO k=1,dv
            swjy(k+sthrednr)=0.d0
         END DO
         DO k=1,n4
            resi(k+rthrednr)=0.d0
         END DO
         DO k=1,dv0
            thi(k+sthrednr)=theta(k,iind)
         END DO
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               DO jw1=1,dlw1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+n1*(j2-1)+n1*n2*(j3-1)
                  IF(mask(jind)) CYCLE
                  wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                  if(wj.le.0.d0) CYCLE
                  si2j=si2(jind)
                  IF (aws) THEN
                     call awswght3(n,jind,dv0,thi(1+sthrednr),
     1                  theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                  END IF
                  if(wlse) THEN 
                     wj=wj*si2j
                  ELSE
                     swjv=swjv+wj/si2j
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k+sthrednr)=swjy(k+sthrednr)+wj*y(k,jind)
                  END DO
                  DO k=1,n4
                     resi(k+rthrednr)=resi(k+rthrednr)+wj*res(k,jind)                     
                  END DO
C                  call daxpy(dv,wj,y(1,jind),1,swjy(1,thrednr),1)
C                  call daxpy(n4,wj,res(1,jind),1,resi(1,thrednr),1)
               END DO
            END DO
         END DO
         DO k=1,dv
            thn(k,iind)=swjy(k+sthrednr)/swj
         END DO
         sresisq=0.d0
         DO k=1,n4
            resik=resi(k+rthrednr) 
            resnew(k,iind)=resik/swj
            sresisq=sresisq+resik*resik
         END DO
         bi(iind)=swj*swj/sresisq*(n4-1.d0)
C  thats the inverse variance
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn,bi,resnew)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ihaws2(y,si2,mask,wlse,n1,n2,n3,dv,dv0,hakt,lambda,
     1                  theta,ncores,bi,thn,kern,skern,spmin,spmax,
     2                  lwght,wght,vwghts,swjy,thi)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi       (input)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer dv,dv0,n1,n2,n3,kern,skern,ncores
      logical aws,wlse,mask(*)
      real*8 theta(dv0,*),bi(*),y(dv,*),
     1       lambda,spmax,wght(2),si2(*),thn(dv,*),
     1       hakt,lwght(*),spmin,vwghts(dv0),thi(dv0,ncores),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,thrednr,iind,jind
      real*8 bii,swj,swjy(dv,ncores),wj,hakt2,spf,si2j,si2i
      external getlwght
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      thrednr = 1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(dv,dv0,n1,n2,n3,kern,skern,ncores,aws,wlse,mask,theta,
C$OMP& bi,y,lambda,spmax,wght,si2,thn,hakt,lwght,spmin,vwghts,thi,
C$OMP& ih1,ih2,ih3,clw1,clw2,clw3,dlw1,dlw2,dlw3,n,swjy,hakt2,spf)
C$OMP& PRIVATE(iind,jind,i1,i2,i3,k,si2i,bii,swj,j1,j2,j3,
C$OMP& jw1,jw2,jw3,wj,si2j,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n
!$         thrednr = omp_get_thread_num()+1
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         if(mask(iind)) THEN
            DO k=1,dv
               thn(k,iind)=0.d0
            END DO
            CYCLE
         END IF
         si2i=si2(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,dv
            swjy(k,thrednr)=0.d0
         END DO
         DO k=1,dv0
            thi(k,thrednr)=theta(k,iind)
         END DO
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               DO jw1=1,dlw1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+n1*(j2-1)+n1*n2*(j3-1)
                  IF(mask(jind)) CYCLE
                  wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                  if(wj.le.0.d0) CYCLE
                  si2j=si2(jind)
                  IF (aws) THEN
                     call awswght3(n,jind,dv0,thi(1,thrednr),
     1                    theta,vwghts,skern,spf,spmin,spmax,bii,wj)
                  END IF
                  if(wlse) THEN 
                     wj=wj*si2j
                  END IF
                  swj=swj+wj
                  call daxpy(dv,wj,y(1,jind),1,swjy(1,thrednr),1)
               END DO
            END DO
         END DO
         DO k=1,dv
            thn(k,iind)=swjy(k,thrednr)/swj
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thn)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smooth3d(y,si2,mask,wlse,n1,n2,n3,dv,hakt,
     1                    thn,kern,lwght,wght,swjy)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      integer n1,n2,n3,kern,dv
      logical wlse,mask(n1,n2,n3)
      real*8 y(n1,n2,n3,dv),thn(n1,n2,n3,dv),wght(2),
     1       si2(n1,n2,n3),hakt,lwght(*),getlwght
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 swj,swjy(dv),wj,hakt2
      external getlwght
      hakt2=hakt*hakt
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  DO k=1,dv
                     thn(i1,i2,i3,k)=0.d0
                  END DO
                  CYCLE
               END IF
C   scaling of sij outside the loop
               swj=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        IF(mask(j1,j2,j3)) CYCLE
                        wj=getlwght(lwght,dlw1,dlw2,dlw3,jw1,jw2,jw3)
                        if(wj.le.0.d0) CYCLE
                        if(wlse) THEN 
                           wj=wj*si2(j1,j2,j3)
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(j1,j2,j3,k)
                        END DO
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(i1,i2,i3,k)=swjy(k)/swj
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
C
C   calculate location weights in lwght
C
      subroutine locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      implicit logical (a-z)
      integer dlw1,dlw2,dlw3,kern
      real*8 wght(2),hakt2,lwght(dlw1,dlw2,dlw3),lkern
      external lkern
      real*8 z1,z2,z3
      integer j1,j2,j3,clw1,clw2,clw3,ih1,ih2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      DO j3=1,dlw3
         Do j2=1,dlw2
            DO j1=1,dlw1
               lwght(j1,j2,j3)=0.d0
            END DO
         END DO
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=FLOOR(sqrt(hakt2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            IF(j2.lt.1.or.j2.gt.dlw2) CYCLE
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO j1=clw1-ih1,clw1+ih1
               IF(j1.lt.1.or.j1.gt.dlw1) CYCLE
               z1=clw1-j1
               lwght(j1,j2,j3)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      RETURN
      END
C
C
C
