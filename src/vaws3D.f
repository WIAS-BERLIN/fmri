CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         lkern=1.d0
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=dexp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          vector-valued  Gaussian   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldist(thi,thj,dv,vwghts)
      implicit logical (a-z)
      integer dv,i
      real*8 thi,thj,z,sz,vwghts(dv)
C        Gaussian
      sz=0.d0
      DO i=1,dv
         z=thi-thj
	 sz=sz+z*z*vwghts(i)
      END DO
      kldist=sz
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws(y,fix,n1,n2,n3,dv,dv0,hakt,lambda,theta,bi,bi2,
     1       bi0,ai,kern,spmin,spmax,lwght,wght,vwghts,swjy,thi,thj)
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
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern,dv,dv0
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,lwght(1),spmin,vwghts(dv0),thi(dv0),thj(dv0)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy(dv),z1,z2,z3,wj,hakt2,
     1        bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      n=n1*n2*n3
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght(2)
            z3=z3*z3
            ih2=dsqrt(hakt2-z3)/wght(1)
            jind3=(j3-1)*dlw1*dlw2
	 ELSE
	    jind3=0
	 END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght(1)
               z2=z3+z2*z2
               ih1=dsqrt(hakt2-z2)
               jind2=jind3+(j2-1)*dlw1
	    ELSE
	       jind2=0
	    END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
	       DO k=1,dv
                  swjy(k)=0.d0
	       END DO
	       DO k=1,dv0
	          thi(k)=theta(iind+(k-1)*n)
	       END DO
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
	          jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=dsqrt(hakt2-z3)/wght(1)
                  DO jw2=clw2-ih2,clw2+ih2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
	             jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        wj=lwght(jw1+jwind2)
                        swj0=swj0+wj
                        IF (aws) THEN
	                DO k=1,dv0
	                   thj(k)=theta(jind+(k-1)*n)
	                END DO
                        sij=bii*kldist(thi,thj,dv0,vwghts)
                           IF (sij.gt.spmax) CYCLE
			   IF (sij.gt.spmin) wj=wj*exp(-sij+spmin)
C   if sij <= spmin  this just keeps the location penalty
C    spmin = 0 corresponds to old choice of K_s 
C   new kernel is flat in [0,spmin] and then decays exponentially
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
			DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(jind+(k-1)*n)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,dv
                  ai(iind+(k-1)*n)=swjy(k)
	       END DO
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws(y,fix,si2,n1,n2,n3,dv,dv0,hakt,lambda,theta,bi,
     1    bi2,bi0,vred,ai,kern,spmin,spmax,lwght,wght,vwghts,swjy,
     2    thi,thj)
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
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern,dv,dv0
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,lwght(1),spmin,vwghts(dv0),thi(dv0),thj(dv0)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy(dv),z1,z2,z3,wj,hakt2,
     1        bii0,si2(1),vred(1),sv1,sv2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      n=n1*n2*n3
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght(2)
            z3=z3*z3
            ih2=dsqrt(hakt2-z3)/wght(1)
            jind3=(j3-1)*dlw1*dlw2
	 ELSE
	    jind3=0
	 END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght(1)
               z2=z3+z2*z2
               ih1=dsqrt(hakt2-z2)
               jind2=jind3+(j2-1)*dlw1
	    ELSE
	       jind2=0
	    END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
	       sv1=0.d0
	       sv2=0.d0
	       DO k=1,dv
                  swjy(k)=0.d0
	       END DO
	       DO k=1,dv0
	          thi(k)=theta(iind+(k-1)*n)
	       END DO
               DO jw3=1,dlw3
	          j3=jw3-clw3+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
		  jwind3=(jw3-1)*dlw1*dlw2
	          jind3=(j3-1)*n1*n2
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=dsqrt(hakt2-z3)/wght(1)
                  DO jw2=clw2-ih2,clw2+ih2
	             j2=jw2-clw2+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
		     jwind2=jwind3+(jw2-1)*dlw1
	             jind2=(j2-1)*n1+jind3
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
	                j1=jw1-clw1+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jind=j1+jind2
                        wj=lwght(jw1+jwind2)
                        swj0=swj0+wj
                        IF (aws) THEN
	                DO k=1,dv0
	                   thj(k)=theta(jind+(k-1)*n)
	                END DO
                        sij=bii*kldist(thi,thj,dv0,vwghts)
                           IF (sij.gt.spmax) CYCLE
			   IF (sij.gt.spmin) wj=wj*exp(-sij+spmin)
C   if sij <= spmin  this just keeps the location penalty
C    spmin = 0 corresponds to old choice of K_s 
C   new kernel is flat in [0,spmin] and then decays exponentially
                        END IF
			sv1=sv1+wj
			sv2=sv2+wj*wj
                        swj=swj+wj*si2(jind)
                        swj2=swj2+wj*wj*si2(jind)
			DO k=1,dv
                          swjy(k)=swjy(k)+wj*y(jind+(k-1)*n)*si2(jind)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,dv
                  ai(iind+(k-1)*n)=swjy(k)
	       END DO
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
	       vred(iind)=sv2/sv1/sv1
            END DO
         END DO
      END DO
      RETURN
      END
