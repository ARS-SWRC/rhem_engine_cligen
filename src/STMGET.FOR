      subroutine stmget(event,day,mon,year,ibrkpt,mxint,dur,
     1 intsty,timem,ninten,ip,timep,p,norain,prcp,rain,stmdur,
     1 nplane,mxplan, filen)
c******************************************************************
c                                                                 *
c  This subroutine is called from MAIN to read in storm data      *
c  on a storm-by-storm basis. It calls SR SUMRNF.                 *
c                                                                 *
c******************************************************************
	implicit none
	integer mxplan

c******************************************************************
      integer  ibrkpt, ninten(75), norain(75), nplane
      integer event, day, mon, year, filen
      real dur, intsty(100), timem(100), ip, timep, p
      real rain(1000), prcp, stmdur
c                                                                 *
c    ibrkpt - flag for breakpoint rainfall inputs 0 - no brkpts   *
c    mxint - maximum rainfall intensity based on average          *
c            intensity value and Ip value                         *
c                                                                 *
c******************************************************************
c
c     real mxint,totrn
      integer nbrkpt, i, jkjk, jkl, kjl
      real xrain, xstmdr, xtimep, xip,mxint,avrint
c
      save
c
      mxint = 0.
	nbrkpt = 0
c
      do 10 kjl =1,nplane
        norain(kjl) = 1
 10   continue

c
	call eatcom (filen)

      if (ibrkpt.eq.0) then

      read (filen,*) event, day, mon, year, xrain, xstmdr, xtimep, xip
c
        do 30 jkjk=1,nplane
          rain(jkjk) = xrain
 30     continue
        prcp = xrain
        stmdur = xstmdr
c
        if (stmdur .gt. 23.999) stmdur = 23.999
c
        timep = xtimep
        ip = xip
c
        if (stmdur.gt.0.0) then
c
          do 40 jkjk=1,nplane
            rain(jkjk) = rain(jkjk) * 0.001
 40       continue
c
          prcp = prcp * 0.001
          stmdur = stmdur * 3600.0
          avrint = rain(1) / stmdur
          mxint = ip * avrint
        else
          if(prcp.gt.0.0)then
            write (6,1300)
            prcp = 0.0
            do 45 jkjk=1,nplane
 45           rain(jkjk) = 0.0
          end if
        end if
c
      else
c       READ IN BREAKPOINT DATA

C	   read (filen,*)  event, day, mon, year, nbrkpt

C	  call eatcom (filen)
c
c       NINT IS THE NUMBER OF BREAKPOINTS OCCURRING ON THE GIVEN DAY.
c       READ IN BREAKPOINT DATA IF RAINFALL OCCURS
c
        if (nbrkpt.gt.0) then

C          call brkpt(nbrkpt,mxint,dur,intsty,timem,p,prcp,stmdur)
c
          do 50 jkjk=1,nplane
            rain(jkjk) = prcp
            ninten(jkjk) = nbrkpt
 50       continue
c
        else
c
          do 60 jkjk=1,nplane
            rain(jkjk) = 0.0
            ninten(jkjk) = 0
 60       continue
c
          prcp = 0.0
        end if
      end if
c
      if (prcp.le.0.0) then
        do 80 jkl=1,nplane
          norain(jkl) = 0
 80     continue
        return
      end if
c
      return

 1300 format ('********************************************************'
     1    /'*                       WARNING!!                      *'/
     1    '* Input precipitation event detected which has a       *'/
     1    '* positive value for precipitation depth and a ZERO    *'/
     1    '* value for precipitation duration.  it has been changed  *'/
     1    '* to 0.0.            *'/
     1    '*                                                      *'/
     1    '* PLEASE CHECK YOUR STORM INPUT FILE                   *'/
     1    ,' *** WARNING ***',/
     1    '*                                                      *'/
     1    '********************************************************')
      end
