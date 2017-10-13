      function sinc(x)
        implicit none
        real*8 :: sinc
        real*8 :: x, xpi
        real*8, parameter :: pi = acos(-1d0)
        xpi = x*pi
        sinc = 1.d0
        if(x /= 0) sinc = sin(xpi)/xpi
      end function sinc

      subroutine sincfilter(t,nt,x,nx,ft)
        implicit none
        integer :: nt, nx
        real*8 :: t(nt), x(nx), ft(nt)
        real*8 :: sinc
        integer :: i, j
        real*8 :: y
        external sinc
        DO i = 1, nt
           y=0.d0
           DO j = 1, nx
             y = y + x(j)*sinc(t(i)-j)
           END DO
           ft(i) = y
        END DO
        return
      END

      subroutine slicetim(x,nt,n1,n2,n3,y,t,sliceord)
        implicit none
        integer :: n1,n2,n3,nt,sliceord(n3)
        real*8 :: x(nt,n1,n2,n3), y(nt,n1,n2,n3), t(nt)
        integer :: i1,i2,i3,it
        real*8 :: rn3,dt
        rn3 = n3
        DO i3=1,n3
           dt=sliceord(i3)-1
           DO it=1,nt
              t(it)=it-dt/rn3
           END DO
           DO i2=1,n2
              DO i1=1,n1
                 call sincfilter(t,nt,x(1,i1,i2,i3),nt,y(1,i1,i2,i3))
              END DO
           END DO
        END DO
        return
      END
