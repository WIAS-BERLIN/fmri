CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segm3d(y,fix,res,si2,mask,wlse,n1,n2,n3,nt,hakt,
     1                  lambda,theta,bi,thn,lwght,wght,swres,pval,
     3                  segm,beta,delta,thresh,step,fov,vq,vest0i,
     4                  varest)
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
      integer n1,n2,n3,nt,kern,segm(n1,n2,n3),step
      logical aws,wlse,mask(n1,n2,n3),fix(n1,n2,n3)
      real*8 y(n1,n2,n3),theta(n1,n2,n3),bi(n1,n2,n3),delta,thresh,
     1      thn(n1,n2,n3),lambda,wght(2),si2(n1,n2,n3),pval(n1,n2,n3),
     1      hakt,lwght(1),thi,getlwght,swres(nt),fov,beta,vq(n1,n2,n3),
     1      varest(n1,n2,n3),res(nt,n1,n2,n3),vest0i(n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n
      real*8 bii,swj,swjy,wj,hakt2,spf,si2j,si2i,swjv,cofh,lsi,s,vqi,
     1       varesti,fpchisq,ti,thij,sij,z,si,swr,z1,lfov,linc,sm1
      external getlwght,fpchisq
      kern=1
      hakt2=hakt*hakt
      lfov=log(fov)
      linc=log(1.25d0)
      s=step
      sm1=s-1
      spf=4.d0/3.d0
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
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
      IF(hakt.gt.1.25) THEN
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  thn(i1,i2,i3)=0.d0
                  CYCLE
               END IF
               vqi=vq(i1,i2,i3)
               thi=theta(i1,i2,i3)
               si2i=vest0i(i1,i2,i3)
               varesti=varest(i1,i2,i3)
               cofh=sqrt(beta*log(varesti*si2i*fov))-0.17d0*sm1
C               lsi = min(sm1,-log(varesti*si2i)/linc)
C               cofh = sqrt(beta*log(varesti*si2i*fov))-0.17d0*lsi
C   this should be more conservative using actual variance reduction instead of theoretical
               ti=max(0.d0,abs(thi)-delta)
               IF(ti/sqrt(varesti/vqi)-cofh.gt.thresh) THEN
                  z=ti/sqrt(varesti/vqi)-cofh-thresh
                  pval(i1,i2,i3)=exp(0.25d0*z*z)
C                  pval(i1,i2,i3)=1.d0-fpchisq(z,1.d0,1,0)
C                  pval(i1,i2,i3)=1.d0
               ELSE
                   pval(i1,i2,i3)=1.d0
C                   pval(i1,i2,i3)=fpchisq(ti*si2i,1.d0,1,0)
               END IF
            END DO
         END DO
      END DO
      END IF 
C   scaling of sij outside the loop
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
                  thn(i1,i2,i3)=0.d0
                  CYCLE
               END IF
               IF (fix(i1,i2,i3)) CYCLE
               vqi=vq(i1,i2,i3)
               si2i=vest0i(i1,i2,i3)
               varesti=varest(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
               swj=0.d0
               swjv=0.d0
               swjy=0.d0
               DO k=1,nt
                  swres(k)=0.d0
               END DO
               thi=theta(i1,i2,i3)
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
                        si2j=si2(j1,j2,j3)
                        IF (aws) THEN
                           thij=thi-theta(j1,j2,j3)
                           sij=thij*thij*bii
                           if(abs(segm(i1,i2,i3)).eq.1) THEN
                           if(segm(i1,i2,i3)*segm(j1,j2,j3).gt.0) THEN
C
C   allow for nonadaptive smoothing if values are in the same segment
C
                                 sij=pval(i1,i2,i3)*sij
                              ELSE
                                 CYCLE
                              END IF
                           END IF
                           IF(sij.gt.1.d0) CYCLE
                        IF(sij.gt.0.25d0) wj=wj*(1.d0-spf*(sij-0.25d0))
                        END IF
                        if(wlse) THEN 
                           wj=wj*si2j
                        ELSE
                           swjv=swjv+wj/si2j
                        END IF
                        swj=swj+wj
                        swjy=swjy+wj*y(j1,j2,j3)
C  weighted sum of residuals
                        call daxpy(nt,wj,res(1,j1,j2,j3),1,swres,1)
                     END DO
                  END DO
               END DO
               z=0.d0
               z1=0.d0
               DO k=1,nt
                  swr=swres(k)/swj
                  z1=z1+swr
                  z=z+swr*swr
               END DO
               thi=swjy/swj
               thn(i1,i2,i3)=thi
               IF(wlse) THEN
                  bi(i1,i2,i3)=swj
               ELSE
                  bi(i1,i2,i3)=swj*swj/swjv
               END IF
               z1=z1/nt
               si = (z/nt - z1*z1)
               si = si/nt
               varest(i1,i2,i3)=si
               cofh = sqrt(beta*log(si*si2i*fov))-0.17d0*step
C               lsi = min(s,-log(si*si2i)/linc)
C               cofh = sqrt(beta*log(si*si2i*fov))-0.17d0*lsi
C   this should be more conservative using actual variance reduction instead of theoretical
               si=sqrt(si/vqi)
               if((thi+delta)/si+cofh.lt.-thresh) THEN
                  segm(i1,i2,i3)=-1
               ELSE IF ((thi-delta)/si-cofh.gt.thresh) THEN
                  segm(i1,i2,i3)=1
C               ELSE
C                  segm(i1,i2,i3)=0
               END IF
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
