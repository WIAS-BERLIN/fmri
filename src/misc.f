      subroutine corr2(res,mask,n1,n2,n3,nv,scorr)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 scorr(n1,n2,n3,3),res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      real*8 z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr
      integer i1,i2,i3,i4
      zk=nv
C  correlation in x
      do i1=1,n1-1
         do i2=1,n2
            do i3=1,n3
               scorr(i1,i2,i3,1)=0.d0
               if (.not.(mask(i1,i2,i3).and.mask(i1+1,i2,i3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+1,i2,i3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) scorr(i1,i2,i3,1)=zcorr/zk/dsqrt(vrm)
            enddo
         enddo
      enddo
C  correlation in y
      do i1=1,n1
         do i2=1,n2-1
            do i3=1,n3
               scorr(i1,i2,i3,2)=0.d0
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2+1,i3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2+1,i3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) scorr(i1,i2,i3,2)=zcorr/zk/dsqrt(vrm)
            enddo
         enddo
      enddo
C  correlation in z
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3-1
               scorr(i1,i2,i3,3)=0.d0
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2,i3+1))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2,i3+1,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) scorr(i1,i2,i3,3)=zcorr/zk/dsqrt(vrm)
            enddo
         enddo
      enddo
      return
      end
      subroutine icorr(res,mask,n1,n2,n3,nv,scorr)

      implicit logical(a-z)
      integer n1,n2,n3,nv,res(n1,n2,n3,nv)
      real*8 scorr(3)
      logical mask(n1,n2,n3)
      real*8 z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,mscorr
      integer i1,i2,i3,i4,k
      zk=nv
C  correlation in x
      k=0
      mscorr=0.d0
      do i1=1,n1-1
         do i2=1,n2
            do i3=1,n3
               if (.not.(mask(i1,i2,i3).and.mask(i1+1,i2,i3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+1,i2,i3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  mscorr=mscorr+zcorr/zk/dsqrt(vrm)
                  k=k+1
               END IF
            enddo
         enddo
      enddo
      scorr(1)=mscorr/k
C  correlation in y
      k=0
      mscorr=0.d0
      do i1=1,n1
         do i2=1,n2-1
            do i3=1,n3
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2+1,i3))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2+1,i3,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  mscorr=mscorr+zcorr/zk/dsqrt(vrm)
                  k=k+1
               END IF
            enddo
         enddo
      enddo
      scorr(2)=mscorr/k
C  correlation in z
      k=0
      mscorr=0.d0
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3-1
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2,i3+1))) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2,i3+1,i4)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/(zk-1.d0)
               vrmp1=y2/(zk-1.d0)
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  mscorr=mscorr+zcorr/zk/dsqrt(vrm)
                  k=k+1
               END IF
            enddo
         enddo
      enddo
      scorr(3)=mscorr/k
      return
      end
      subroutine ivar(res,resscale,mask,n1,n2,n3,nv,var)

      implicit logical(a-z)
      integer n1,n2,n3,nv,res(n1,n2,n3,nv)
      real*8 resscale,var(n1,n2,n3)
      logical mask(n1,n2,n3)
      real*8 z2,zk,resi,ressc2
      integer i1,i2,i3,i4
      zk=nv
      ressc2=resscale*resscale
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3
               var(i1,i2,i3)=1.d20
               if (.not.mask(i1,i2,i3)) CYCLE
               z2=0.d0
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  z2=z2+resi*resi
               enddo
               var(i1,i2,i3)=z2/(zk-1.d0)*ressc2
            enddo
         enddo
      enddo
      return
      end
      subroutine corr(res,mask,n1,n2,n3,nv,scorr)

      implicit logical(a-z)
      integer n1,n2,n3,nv
      real*8 scorr(3),res(n1,n2,n3,nv)
      logical mask(n1,n2,n3)
      real*8 z,z2,resi,resip1,vrm,zk,zcorr
      integer i1,i2,i3,i4,k
      
      scorr(1)=0.d0
      scorr(2)=0.d0
      scorr(3)=0.d0
C  correlation in x
      z=0.d0
      z2=0.d0
      zcorr=0.d0
      k=0
      do i1=1,n1-1
         do i2=1,n2
            do i3=1,n3
               if (.not.(mask(i1,i2,i3).and.mask(i1+1,i2,i3))) CYCLE
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1+1,i2,i3,i4)
                  if (resi.eq.0.d0.or.resip1.eq.0.d0) CYCLE
                  z=z+resi
                  z2=z2+resi*resi
                  zcorr=zcorr+resi*resip1
                  k=k+1
               enddo
            enddo
         enddo
      enddo
      if (k.gt.0) THEN
         zk=k
         z=z/zk
         vrm=(z2-zk*(z*z))/(zk-1.d0)
         scorr(1)=zcorr/zk/vrm
      ELSE
         scorr(1)=0.d0
         call dblepr("All zero residuals x",20,scorr(1),1)
      END IF
      if(scorr(1).gt.1.d0-1.d-8) THEN
         call dblepr("scorr(1) was",12,scorr(1),1)
         scorr(1) = 1.d0-1.d-8
         call dblepr("scorr(1) reset to",17,scorr(1),1)
      ELSE IF (scorr(1).lt.-1.d0+1.d-8) THEN
         call dblepr("scorr(1) was",12,scorr(1),1)
         scorr(1) = 1.d0-1.d-8
         call dblepr("scorr(1) reset to",17,scorr(1),1)
      END IF
C  correlation in y
      z=0.d0
      z2=0.d0
      zcorr=0.d0
      k=0
      do i1=1,n1
         do i2=1,n2-1
            do i3=1,n3
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2+1,i3))) CYCLE
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2+1,i3,i4)
                  if (resi.eq.0.d0.or.resip1.eq.0.d0) CYCLE
                  z=z+resi
                  z2=z2+resi*resi
                  zcorr=zcorr+resi*resip1
                  k=k+1
               enddo
            enddo
         enddo
      enddo
      if (k.gt.0) THEN
         zk=k
         z=z/zk
         vrm=(z2-zk*(z*z))/(zk-1.d0)
         scorr(2)=zcorr/zk/vrm
      ELSE
         scorr(2)=0.d0
         call dblepr("All zero residuals y",20,scorr(2),1)
      END IF
      if(scorr(2).gt.1.d0-1.d-8) THEN
         call dblepr("scorr(2) was",12,scorr(2),1)
         scorr(2) = 1.d0-1.d-8
         call dblepr("scorr(2) reset to",17,scorr(2),1)
      ELSE IF (scorr(2).lt.-1.d0+1.d-8) THEN
         call dblepr("scorr(2) was",12,scorr(2),1)
         scorr(1) = 1.d0-1.d-8
         call dblepr("scorr(2) reset to",17,scorr(2),1)
      END IF
C  correlation in z
      z=0.d0
      z2=0.d0
      zcorr=0.d0
      k=0
      do i1=1,n1
         do i2=1,n2
            do i3=1,n3-1
               if (.not.(mask(i1,i2,i3).and.mask(i1,i2,i3+1))) CYCLE
               do i4=1,nv
                  resi=res(i1,i2,i3,i4)
                  resip1=res(i1,i2,i3+1,i4)
                  if (resi.eq.0.d0.or.resip1.eq.0.d0) CYCLE
                  z=z+resi
                  z2=z2+resi*resi
                  zcorr=zcorr+resi*resip1
                  k=k+1
               enddo
            enddo
         enddo
      enddo
      if (k.gt.0) THEN
         zk=k
         z=z/zk
         vrm=(z2-zk*(z*z))/(zk-1.d0)
         scorr(3)=zcorr/zk/vrm
      ELSE
         scorr(3)=0.d0
         call dblepr("All zero residuals x",20,scorr(3),1)
      END IF
      if(scorr(3).gt.1.d0-1.d-8) THEN
         call dblepr("scorr(3) was",12,scorr(3),1)
         scorr(3) = 1.d0-1.d-8
         call dblepr("scorr(3) reset to",17,scorr(3),1)
      ELSE IF (scorr(3).lt.-1.d0+1.d-8) THEN
         call dblepr("scorr(3) was",12,scorr(3),1)
         scorr(3) = 1.d0-1.d-8
         call dblepr("scorr(3) reset to",17,scorr(3),1)
      END IF
      return
      end
