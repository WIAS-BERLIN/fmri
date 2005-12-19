      subroutine eigen3(ddiff, diff, evectors, evalues, evaonly, ierr)
      integer ddiff, ierr(ddiff)
      real*8 diff(ddiff,6), evectors(9,ddiff), evalues(3,ddiff)
      logical evaonly
      real*8 xmat(9),work1(3),work2(3)
      integer i,n
      n=3
      DO i=1,ddiff
         xmat(1)=diff(i,1)
         xmat(2)=diff(i,4)
         xmat(3)=diff(i,5)
         xmat(4)=diff(i,4)
         xmat(5)=diff(i,2)
         xmat(6)=diff(i,6)
         xmat(7)=diff(i,5)
         xmat(8)=diff(i,6)
         xmat(9)=diff(i,3)
	 call rs(n, n, xmat, evalues(1,i), evaonly, evectors(1,i), 
     1            work1, work2, ierr(i)) 
      END DO 
      RETURN
      END