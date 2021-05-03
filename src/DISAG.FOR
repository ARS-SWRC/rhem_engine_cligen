	 subroutine disag(dur,intsty,timem,ninten,intdl,timedl,
     1                  deltfq,fq,p,ip,timep,rain,stmdur,iplane,mxplan)
c
c     + + + PURPOSE + + +
c     Disaggregates storms into a double exponential intensity
c     pattern with relative time to peak intensity, Tp, and
c     relative maximum intensity, Ip = Max INT/AVE intensity,
c     satisfying 0 < Tp < 1 and Ip >= 1.
c     Also reads precipitation amount P, and storm duration DUR.
c
c     Called from IDAT.
c     Author(s): Lane, Lopez, Stone, FERRIS
c     Updated by: Dennis C. Flanagan and Reza Savabi
c
c     + + + KEYWORDS + + +
c
c     + + + ARGUMENT DECLARATIONS + + +
      
	implicit none
c
c     + + + ARGUMENT DEFINITIONS + + +
	integer mxplan,ninten(mxplan),iplane
	real intsty(100),timem(100),dur
	real intdl(21),timedl(21)
	real deltfq,fq,p,ip,timep
	real rain(mxplan),stmdur

c     + + + LOCAL VARIABLES + + +
      real frstim, sectim
      integer loopfg, i, im
c
c     + + + LOCAL DEFINITIONS + + +
c     frstim - real value for previous time when checking if
c              disaggregated time step is less than 5 minutes
c     sectim - real value for current time when checking if
c              disaggregated time step is less than 5 minutes
c     loopfg - flag.  1=exit loop; 2=exit loop and go to 20;
c                     0=next iteration of loop.
c
c     + + + SUBROUTINES CALLED + + +
c     const
c     dblex
c
c     + + + END SPECIFICATIONS + + +
c
c     initialize variables
c
C      open(725,file='rainfall_int_24Jul2013.txt',status='unknown')
C      ninten(iplane) = 11 original statement
      ninten(iplane) = 21
      deltfq = 1.0 / float(ninten(iplane)-1)
      fq = 0.0
      timedl(1) = 0.0
      intdl(ninten(iplane)) = 0.0
      p = 0.
c
c     ****************************
c     Truth Table of       *
c     Cases Covered:       *
c     *
c     | IR | SM | RN | Case  *
c     ---------------------  *
c     |  n |  N |  Y |  1    *
c     |  n |  Y |  Y |  2    *
c     |  n |  Y |  N |  2    *
c     |  y |  N |  n |  3    *
c     |  y |  N |  y |  4    *
c     *
c*****************************
c      STMDUR = DURATION OF RAINFALL IF RAINFALL OCCURS (STMGET.FOR)

c    RAINFALL ONLY - no irrigation OR SNOWMELT
c
          p = rain(iplane)
          dur = stmdur

c
      if(ip.lt.1.0)ip = 1.0
      if (timep.gt.1.0.or.ip.eq.1.0) then
        timep = 1.0
      else if (timep.le.0.0)then
        timep = 0.01
      else
        continue
      end if
c
   20 continue
      loopfg = 0
c
c     Call CONST if intensity is constant,
c     call DBLEX if it is not.
c
      if (timep.ge.1.0.and.ip.le.1.0) then
        call const(ninten,intdl,timedl,deltfq,iplane,mxplan)
      else
        call dblex(ninten,intdl,timedl,ip,timep,deltfq,iplane,mxplan)
      end if
c
!      frstim = timedl(1) * dur
!c
!      i = 1
!   30 continue
!      i = i + 1
!      sectim = timedl(i) * dur
!      if (sectim-frstim.lt.300.) then
!c
!c       Time step is less than 5 minutes.  Decrease the number
!c       of dimensionless steps and try again.
!c
!        ninten(iplane) = ninten(iplane) - 1
!c
!        if (ninten(iplane).le.2) then
!c
!c         Disaggregated rainfall distribution is set to
!c         constant intensity with 2 time steps.
!c
!          timedl(2) = 1.
!          intdl(1) = 1.
!          intdl(2) = 1.
!          ninten(iplane) = 2
!        else
!c
!c         Re-initialize for decreased step.
!c
!          deltfq = 1.0 / float(ninten(iplane)-1)
!          fq = 0.0
!          timedl(1) = 0.0
!          intdl(ninten(iplane)) = 0.0
!          loopfg = 2
!        end if
!      else
!        frstim = sectim
!      end if
!      if (i.lt.ninten(iplane).and.loopfg.eq.0) go to 30
!      if (loopfg.eq.2) go to 20
!c
!c     Make sure last intensity value is 0.0
!c
      intdl(ninten(iplane)) = 0.0

c     Calculate actual time and intensity
c
      do 60 i = 1, ninten(iplane)
        timem(i) = timedl(i) * dur
c
c       changed p and dur to meters and second above so changed here.
c
        intsty(i) = intdl(i) * p / dur
C-------------------------------------
C--------- Mariano Hernandez 22May2012
C	  write(725,1725)i, timem(i)/60., intsty(i)*1000.*3600.
C1725    format(i5,f20.10,f20.10)
C-------------------------------------
   60 continue
c
      return
      end
