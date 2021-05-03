      subroutine const(ninten,intdl,timedl,deltfq,iplane,mxplan)
c
c     + + + PURPOSE + + +
c     Generates A step function to represent delta-t and delta--
c     intensity for cases where intensity is constant.  (Updates
c     INTDL(20) & TIMEDL(20). )
c
c     Note: Either CONST is called from DISAG, or DBLEX is, but
c           NOT BOTH!
c
c     Called from DISAG.
c     Author(s): Lane, Lopez, Stone
c     Reference in User Guide:
c
c     Changes:
c           1) FQX added to store the value of FQ, to reduce
c              flow of data to module.
c           2) FQ is simply initialized to zero in DISAG.  It
c              was removed from this module, and a zero was
c              substituted.
c
c     Recoded by: Charles R. Meyer.
c
c     + + + KEYWORDS + + +
c
c     + + + PARAMETERS + + +
      
	implicit none
	integer mxplan,ninten(mxplan),iplane
	real intdl(21),timedl(21),deltfq
c
c     + + + LOCAL VARIABLES + + +
      real fqx
      integer i
c	real deltfq
c
c     + + + LOCAL DEFINITIONS + + +
c     fqx    - temporarily holds DELTFQ*I.
c
c     + + + END SPECIFICATIONS + + +
c
      fqx = 0.0
c
      do 10 i = 1, ninten(iplane) - 1
        fqx = fqx + deltfq
        timedl(i+1) = fqx
        intdl(i) = 1.0
   10 continue
c
      return
      end
