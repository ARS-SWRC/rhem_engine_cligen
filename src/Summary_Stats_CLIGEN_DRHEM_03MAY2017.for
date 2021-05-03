C     Code type: FORTRAN subroutine
      
      subroutine summary_stats(nevent, iEvent, iDay, iMonth,
     & iYear, Area, Rainfall, Infil, Runoff, PeakR,
     & SoilLoss, Sloss, depo, sy, susp, file_unit, file_name )
      
C     Compiler: Microsoft Visual Studio 2010

C     Programmed by: Mariano Hernandez

C     Date: 18 JAN 2013

C     Description: 

C     Computes summary statistics, and the 2-yr, 5-yr, 10-yr, 25-yr, 50-yr,
C     and 100-yr return periods for maximum annual soil loss events.
C------------------------------------------------------------------------------
C     Argument Declarations:
      
      integer iEvent, iDay, iMonth, iYear, nevent, files
      real Area, Rainfall, Infil, Runoff, PeakR, SoilLoss, SLoss, depo,
     &     sy, susp
      character*150 fields
C------------------------------------------------------------------------------

C     Argument Definitions:
      
C         iEvent      -   
C         iDay        -
C         iMonth      -
C         iYear       -
C         Area        -
C         Rainfall    -
C         Infil       -
C         Runoff      -
C         PeakR       -
C         SoilLoss    -
C         Depo        - deposition
C         SLoss       - soil loss from integration of depnet curve
C         sy          - sediment yield
C         susp        - suspended soil
C------------------------------------------------------------------------------

C     Local Variables:

      integer i,n,nYr, indBeg, indFin, nRP, File_Error, file_unit, AR1
      integer nSoilLossTolerance, kount
      real AnnualRainfall, AnnualInfil, AnnualRunoff, AnnualSoilLoss
      real AnnualRainfallAverage, AnnualDepo, AnnualSloss, AnnualSY
      real AnnualDepoAverage, AnnualSLossAverage, AnnualSYAverage
      real MaxRainfallEventYr, MaxRunoffEventYr, MaxSoilLossEventYr
      real MaxSLossEventYr, MaxSYEventYr, MaxDepoEventYr
      real arrayRP, prob, ReturnPeriod, SortAscReturnPeriod,
     &     yRainfall, yRunoff, ySoilLoss, yDepo, ySloss, ySY,
     &     outRainfall,outRainfallAnnual, outRunoff, outRunoffAnnual,
     &     outSoilLoss, outDepo, outSLoss,outSLossAnnual, outSY,
     &     outSYAnnual, y2Rainfall,y2Runoff, y2SoilLoss,y2Depo,y2SLoss,
     &     y2SY, SoilLossTolerance, outAnnualSoilLoss,
     &     y2RainfallAnnual, y2RunoffAnnual, y2SLossAnnual,y2SYAnnual,
     &     prob_up205, prob_05_10, prob_10_30, prob_over3,
     &     AnnualSoilLoss4Risk, y2AnnualSoilLoss4Risk,
     &     outAnnualSoilLoss4Risk, yAnnualSoilLoss4Risk
      character(len=150) :: file_name

C     Subroutines Called:
C         sort
C         spline
C         splint


      dimension iEvent(110000), iDay(110000), iMonth(110000),
     &        Rainfall(110000), Infil(110000), iYear(110000),
     &        Runoff(110000), PeakR(110000), SoilLoss(110000),
     &        SLoss(110000), SY(110000), susp(110000),
     &        Depo(110000), kountP(110000),
     &        n(300), AnnualRainfall(300), indBeg(300), indFin(300),
     &        AnnualRunoff(300), AnnualSoilLoss(300), AnnualDepo(300),
     &        AnnualInfil(300), MaxRunoffEventYr(300),AnnualSLoss(300),
     &        AnnualSY(300), AnnualSusp(300), AnnualSoilLoss4Risk(300), 
     &        MaxSoilLossEventYr(300), prob(300), returnPeriod(300),
     &        SortAscReturPeriod(300), arrayRP(6), yRainfall(6),
     &        yRainfallAnnual(6), yRunoffAnnual(6),ySLossAnnual(6), 
     &        yRunoff(6), ySoilLoss(6), yDepo(6), ySLoss(6), ySY(6),
     &        ySYAnnual(6), ySusp(6), MaxSuspEventYr(300),
     &        MaxRainfallEventYr(300), MaxDepoEventYr(300),
     &        MaxSLossEventYr(300), MaxSYEventYr(300),
     &        y2Rainfall(300),y2RainfallAnnual(300),
     &        y2Runoff(300),y2RunoffAnnual(300), y2SoilLoss(300),
     &        y2Depo(300), y2SLoss(300), y2SLossAnnual(300),
     &        y2SY(300), y2SYAnnual(300), y2Susp(300),
     &        arraySLT(3), R(300), rank(300),
     &        AR1(300), y2AnnualSoilLoss(300), SoilLossTolerance(3),
     &        AnnualNew(300), jj(300),
     &        y2AnnualSoilLoss4Risk(300),yAnnualSoilLoss4Risk(3) 
      
      logical :: jj      
C------------------------------------------------------------------------------      
C     Subroutines Called:
C         sort
C         spline
C         splint
      
C     + + + OPEN FILE +  + +
C      open(725, file='CligenOutPut.dat', status='unknown')
      open(file_unit, file = file_name ,status='unknown')
C     + + + DATA INITIALIZATION + + +

      yp1 = 1.E30
      ypn = 1.E30
      arrayRP(1) = 2.
      arrayRP(2) = 5.
      arrayRP(3) = 10.
      arrayRP(4) = 25.
      arrayRP(5) = 50.
      arrayRP(6) = 100.
      nRP = 6
C --- BEGIN Mariano Hernandez 11Nov2013 -----------      
C --- SLT = soil loss tolerance      
      nSoilLossTolerance = 3
C      arraySLT(1) = 0.5
C      arraySLT(2) = 1.0
C      arraySLT(3) = 3.0
C --- soil loss tolerance changes on 05Dec 2013      
      arraySLT(1) = 0.25
      arraySLT(2) = 0.50
      arraySLT(3) = 1.50
C --- END -----------------------------------------
            
C      p = 0
      
C      do
C          p = p + 1  
C          read(725, 10, iostat = file_error) iEvent(p), iDay(p),
C     &     iMonth(p), iYear(p), Area(p), Rainfall(p), Infil(p),
C     &     Runoff(p), PeakR(p), SoilLoss(p)
C          
C          if(file_error .ne. 0) exit
          
C      end do

C     + + +  Find last year of simulation + + +

      nYr =iYear(nevent)
C ------------------------------------------------------------------
C --- 03June2013 Mariano Hernandez ---------------------------------

      do i = 1, nYr
          kountP(i) = 0
      end do
      
      do j = 1, nevent
          
          do i = 1, nYr
              
              if ( iYear(j) == i ) then
                  kountP( i ) = kountP( i ) + 1
              else
                  
              endif
              
          end do
          
      end do
      
      numeroYr = COUNT(kountP .NE. 0, DIM=1)
      
      do i = 1, nYr
          
          idummy = kountP( i )

          if ( idummy .NE. 0 ) then
          
              n( i ) = 0
             
              do j = 1, nevent
            
                 if ( iyear( j ) == i ) then
              
                  n( i ) = n( i ) + 1
              
C                 indFin ( i ) = ievent( j )  
                  indFin ( i ) = j
                 else
              
                 endif
          
              end do
          
          indBeg ( i ) = indFin ( i ) - n( i ) + 1
          else
          end if
          
      end do 
      
C     + + + Annual Estimate Values + + +
C     + + Rainfall, Runoff, Soil Loss + +

      do  k = 1, nevent
          
          if( SoilLoss(k) .gt. SLoss(k) )then
              SoilLoss(k) = SLoss(k)
          endif
          
      end do
      
      do i = 1, nYr
          idummy = kountP(i)
          if ( idummy .NE. 0 )then
      
          AnnualRainfall(i) = sum( Rainfall( indBeg(i) : indFin(i) ) )
          AnnualInfil(i) = sum( Infil( indBeg(i) : indFin(i) ) )
          AnnualRunoff(i) = sum( Runoff( indBeg(i) : indFin(i) ) )
          AnnualSoilLoss(i) = sum( SoilLoss( indBeg(i) : indFin(i) ) )
          AnnualDepo(i) = sum( Depo( indBeg(i) : indFin(i) ) )
C --- SLoss (kg/ha) = soil loss   ---------------------------------------------          
          AnnualSLoss(i) = sum( abs( SLoss( indBeg(i) : indFin(i) ) ) )
C ---------------------------------------------------------------------
C --- SY (kg/ha) = sediment yield   -------------------------------------------           
          AnnualSY(i) = sum( SY( indBeg(i) : indFin(i) ) )
C ---------------------------------------------------------------------          
          AnnualSusp(i) = sum( susp( indBeg(i) : indFin(i) ) )

      MaxRainfallEventYr(i) = maxval( Rainfall( indBeg(i) : indFin(i) ))
      MaxRunoffEventYr(i) = maxval( Runoff( indBeg(i) : indFin(i) ) )
      MaxSoilLossEventYr(i) = maxval( SoilLoss( indBeg(i) : indFin(i) ))
      MaxDepoEventYr(i) = maxval( Depo( indBeg(i) : indFin(i) ))
      MaxSLossEventYr(i) = maxval( SLoss( indBeg(i) : indFin(i) ) ) ! (kg/ha)
      MaxSYEventYr(i) = maxval( SY( indBeg(i) : indFin(i) ) )
      MaxSuspEventYr(i) = maxval( Susp( indBeg(i) : indFin(i) ) )
          else
          end if
          
      end do
      
725   format(1x, i5, 3F20.10)      
C     + + + Annual Average Values + + +

      AnnualRainfallAverage = ( (sum( AnnualRainfall )/numeroYr)/Area )*
     &                                                             1000.
      
      AnnualRunoffAverage = ( (sum( AnnualRunoff )/numeroYr )/Area )   *
     &                                                             1000.
       
      AnnualSoilLossAverage = sum( AnnualSoilLoss )/numeroYr
      
      AnnualDepoAverage = ( sum( AnnualDepo )/numeroYr )
      
      AnnualSLossAverage = ( sum( AnnualSLoss )/numeroYr ) ! Average Soil Loss (kg/ha) in Summary Output (ton/ha) 

      AnnualSYAverage = ( sum( AnnualSY )/numeroYr ) ! Average Sediment Yield (kg/ha) in Summary Output (ton/ha) 
      
      AnnualSuspAverage = ( sum( AnnualSusp )/numeroYr )
C---BEGIN Mariano Hernandez 03May2017
C---Sort Annual Totals
      call sort( numeroYr, AnnualRainfall )
      call sort( numeroYr, AnnualRunoff )
      call sort( numeroYr, AnnualSLoss )
      call sort( numeroYr, AnnualSY )      
C---END 03May2017
      
      call sort( numeroYr, MaxRainfallEventYr )
      call sort( numeroYr, MaxRunoffEventYr )
      call sort( numeroYr, MaxSoilLossEventYr )
      call sort( numeroYr, MaxDepoEventYr )
      call sort( numeroYr, MaxSLossEventYr )
      call sort( numeroYr, MaxSYEventYr )
C --- BEGIN Mariano Hernnadez 31Oct2013 ---------------      
C --- Risk Assessment
      AnnualSoilLoss4Risk=AnnualSLoss/1000.    ! Annual Soil Loss (kg/ha); AnnualSoilLoss4Risk (ton/ha)
C     AnnualSoilLoss4Risk=AnnualSoilLoss/1000. ! Annual Sediment Yield
      call sort( numeroYr, AnnualSoilLoss4Risk )
C ----------------------------------------------------
      AR1 = MINLOC(AnnualSoilLoss4Risk,DIM=1,
     &             MASK= AnnualSoilLoss4Risk .GT. 0.)
      
      minValue = minval(AR1)
      
      kount = 0
      
      do i = minValue, numeroYr
          kount = kount + 1
          AnnualNew(kount) = AnnualSoilLoss4Risk(i) ! AnnualNew (ton/ha); AnnualSoilLoss4Risk (ton/ha)
      end do
      
C --- Estimation of the Kaplan-Meier Estimator 
      do k = 1, kount
          rank(k) = kount - (k-1)
      end do
      
      R(1) = 1.
      
      do i = 2, kount
          R(i) = ( ( rank(i) - 1. )/rank(i) ) * R(i-1)
      end do
      
      do k = 2, kount
          
          SoilLossDiff = abs( AnnualNew(k-1) - AnnualNew(k) )
          if ( SoilLossDiff .LT. 0.0001) then
              jj(k) = .TRUE.
          else
              jj(k) = .FALSE.
          end if
          
      end do
      
      do k = 2, kount
          
          if ( jj(k) .EQ. .TRUE.) then
              AnnualNew(k-1) = AnnualNew(k) * 0.999
          else
              AnnualNew(k)   = AnnualNew(k) * 1.00
          end if
          
      end do
      
C --- BEGIN Mariano Hernandez 22Nov2013
C --- Check if AnnualNew MAX is < than 0.5 ton/ha (FIRST SOIL LOSS TOLERANCE)
C --- Set Prob = 100
      
       if ( AnnualNew(1) .GT. arraySLT(3) ) then
           
            yAnnualSoilLoss4Risk(1)= 1.
            yAnnualSoilLoss4Risk(2)= 1.
            yAnnualSoilLoss4Risk(3)= 1.
      
       else if ( AnnualNew(kount) .LE. arraySLT(1) ) then
          
            yAnnualSoilLoss4Risk(1)= 0.
            yAnnualSoilLoss4Risk(2)= 0.
            yAnnualSoilLoss4Risk(3)= 0.
            
       else if ( AnnualNew(1) .GT. arraySLT(2) .AND.
     &           AnnualNew(1) .LE. arraySLT(3) ) then
                 
            yAnnualSoilLoss4Risk(1) = 1.
            yAnnualSoilLoss4Risk(2) = 1.
            
       do m = 3, nSoilLossTolerance        
          call spline( AnnualNew(1:kount), R(1:kount),
     &            kount, yp1, ypn, y2AnnualSoilLoss4Risk )
          
          call splint( AnnualNew(1:kount), R(1:kount),
     &     y2AnnualSoilLoss4Risk, kount, arraySLT(m),
     &     outAnnualSoilLoss4Risk)     
          yAnnualSoilLoss4Risk(m) = outAnnualSoilLoss4Risk         
       end do
       
      else if ( AnnualNew(1) .GT. arraySLT(1) .AND.
     &          AnnualNew(1) .LE. arraySLT(2) ) then
      
           yAnnualSoilLoss4Risk(1)= 1.
      
      do m = 2, nSoilLossTolerance        
          call spline( AnnualNew(1:kount), R(1:kount),
     &            kount, yp1, ypn, y2AnnualSoilLoss4Risk )
          
          call splint( AnnualNew(1:kount), R(1:kount),
     &     y2AnnualSoilLoss4Risk, kount, arraySLT(m),
     &     outAnnualSoilLoss4Risk)     
           yAnnualSoilLoss4Risk(m) = outAnnualSoilLoss4Risk        
      end do
      
      else if ( AnnualNew(1) .GT. arraySLT(1) .AND.
     &          AnnualNew(1) .LE. arraySLT(2) ) then
          
      do m = 1, nSoilLossTolerance       
          call spline( AnnualNew(1:kount), R(1:kount),
     &            kount, yp1, ypn, y2AnnualSoilLoss4Risk )
          
          call splint( AnnualNew(1:kount), R(1:kount),
     &     y2AnnualSoilLoss4Risk, kount, arraySLT(m),
     &     outAnnualSoilLoss4Risk)     
           yAnnualSoilLoss4Risk(m) = outAnnualSoilLoss4Risk     
      end do
      
      
      else if ( AnnualNew(1) .LE. arraySLT(1) ) then
          
      do m = 1, nSoilLossTolerance       
          call spline( AnnualNew(1:kount), R(1:kount),
     &            kount, yp1, ypn, y2AnnualSoilLoss4Risk )
          
          call splint( AnnualNew(1:kount), R(1:kount),
     &     y2AnnualSoilLoss4Risk, kount, arraySLT(m),
     &     outAnnualSoilLoss4Risk)     
           yAnnualSoilLoss4Risk(m) = outAnnualSoilLoss4Risk     
      end do
          
      end if 
      
C --- BEGIN Mariano Hernandez 15Nov2013 -------------------------------
C --- mantain a monotonically decreasing function to prevent negative 
C --- probability values      
      if (yAnnualSoilLoss4Risk(2) .LT. yAnnualSoilLoss4Risk(3) ) then
          
          diff_yAnnualSoilLoss4Risk =( yAnnualSoilLoss4Risk(3) -
     &                                 yAnnualSoilLoss4Risk(2) )/2.
          yAnnualSoilLoss4Risk(3) = yAnnualSoilLoss4Risk(2) - 
     &                              diff_yAnnualSoilLoss4Risk     
      else
          yAnnualSoilLoss4Risk(3) = yAnnualSoilLoss4Risk(3)
      end if
C ---END --------------------------------------------------------------      
      
      
C --- Estimation of the Probability of Occurance for Risk Assessment ---
C --- First Case when: P(X < 0.5)
C --- BEGIN Mariano Hernandez 10Dec2013
      if ( yAnnualSoilLoss4Risk(1) .LT. 0.0 )
     &                                     yAnnualSoilLoss4Risk(1) = 0.
      if ( yAnnualSoilLoss4Risk(1) .GT. 1.0 )
     &                                     yAnnualSoilLoss4Risk(1) = 1.
C --- END 10Dec2013      
      prob_up205 = 1. - yAnnualSoilLoss4Risk(1)
C --- Second Case when: P( >= 0.5 X <= 1.0)
      if ( yAnnualSoilLoss4Risk(2) .LT. 0.0 )
     &                                     yAnnualSoilLoss4Risk(2) = 0.
      if ( yAnnualSoilLoss4Risk(2) .GT. 1.0 )
     &                                     yAnnualSoilLoss4Risk(2) = 1.
C      prob_05_10 = ( 1. - yAnnualSoilLoss4Risk(2) ) - 
C     &             ( 1. - yAnnualSoilLoss4Risk(1) )
      prob_05_10 = yAnnualSoilLoss4Risk(1) - yAnnualSoilLoss4Risk(2)
C --- Third Case when: P( >= 1.0 X <= 3.0 )
      if ( yAnnualSoilLoss4Risk(2) .LT. 0.0 )
     &                                     yAnnualSoilLoss4Risk(2) = 0.
      if ( yAnnualSoilLoss4Risk(3) .LT. 0.0 .OR.
     &     yAnnualSoilLoss4Risk(3) .GT. 1.0 )
     &                                     yAnnualSoilLoss4Risk(3) = 0.
C      prob_10_30 = ( 1. - yAnnualSoilLoss4Risk(3) ) - 
C     &             ( 1. - yAnnualSoilLoss4Risk(2) )
      prob_10_30 = yAnnualSoilLoss4Risk(2) - yAnnualSoilLoss4Risk(3)
C --- Fourth Case when: P( X > 3.0 )
      if ( yAnnualSoilLoss4Risk(3) .LT. 0.0 .OR.
     &     yAnnualSoilLoss4Risk(3) .GT. 1.0 )
     &                                     yAnnualSoilLoss4Risk(3) = 0.
      prob_over3 = yAnnualSoilLoss4Risk(3)
C --- Check sum probabilities 
      test = prob_up205 + prob_05_10 + prob_10_30 + prob_over3
C --- END Mariano Hernandez 11Nov2013 ----------------------------------      
C      else
C          prob_up205 = 1.
C          prob_05_10 = 0.
C          prob_10_30 = 0.
C          prob_over3 = 0.
C      end if
C --- END Mariano Hernandez 22Nov2013     
      
C     + + + Compute  Probability + + +
C100   continue
      do i = 1, numeroYr
          
         prob( i ) = i/( numeroYr + 1. )
         returnPeriod ( i ) = 1./prob( i )
         
      end do
      
      SortAscReturPeriod( 1 ) = ( numeroYr )
      
      do i = 1, numeroYr-1
          
          SortAscReturPeriod( i+1 ) =  returnPeriod( numeroYr - i )
          
      end do

C---BEGIN Mariano Hernandez  03May2017
C---Yearly Totals Frequency Analysis
      do m = 1, nRP
          
C         call SPLINE to get second derivatives

         call spline( SortAscReturPeriod, AnnualRainfall,
     &                 numeroYr, yp1, ypn, y2RainfallAnnual)
         
         call spline( SortAscReturPeriod, AnnualRunoff,
     &                 numeroYr, yp1, ypn, y2RunoffAnnual)
        
        call spline( SortAscReturPeriod, AnnualSLoss,
     &                 numeroYr, yp1, ypn, y2SLossAnnual)
          
        call spline( SortAscReturPeriod, AnnualSY,
     &                 numeroYr, yp1, ypn, y2SYAnnual)        

C         call SPLINT for interpolations

          call splint( SortAscReturPeriod, AnnualRainfall,
     &    y2RainfallAnnual, numeroYr, arrayRP(m), outRainfallAnnual )

          call splint( SortAscReturPeriod, AnnualRunoff,
     &    y2RunoffAnnual, numeroYr, arrayRP(m), outRunoffAnnual )
          
          call splint( SortAscReturPeriod, AnnualSLoss,
     &    y2SLossAnnual, numeroYr, arrayRP(m), outSLossAnnual )
          
          call splint( SortAscReturPeriod, AnnualSY,
     &    y2SYAnnual, numeroYr, arrayRP(m), outSYAnnual )          
          
          yRainfallAnnual (m) = ( outRainfallAnnual / Area ) * 1000. ! (mm)
          yRunoffAnnual (m) = ( outRunoffAnnual / Area ) * 1000. ! (mm)
          ySLossAnnual (m) = outSLossAnnual ! Soil Loss (Mg/ha)
          ySYAnnual (m) =  outSYAnnual ! Sediment Yield (Mg/ha)
      
      end do

C---END Yearly Totals Frequency Analysis
      
      do m = 1, nRP
          
C         call SPLINE to get second derivatives

         call spline( SortAscReturPeriod, MaxRainfallEventYr,
     &                 numeroYr, yp1, ypn, y2Rainfall)
         
         call spline( SortAscReturPeriod, MaxRunoffEventYr,
     &                 numeroYr, yp1, ypn, y2Runoff)

        call spline( SortAscReturPeriod, MaxSoilLossEventYr,
     &                 numeroYr, yp1, ypn, y2SoilLoss)

       call spline( SortAscReturPeriod, MaxDepoEventYr,
     &                 numeroYr, yp1, ypn, y2Depo)
        
        call spline( SortAscReturPeriod, MaxSLossEventYr,
     &                 numeroYr, yp1, ypn, y2SLoss)
          
        call spline( SortAscReturPeriod, MaxSYEventYr,
     &                 numeroYr, yp1, ypn, y2SY)        

C         call SPLINT for interpolations

          call splint( SortAscReturPeriod, MaxRainfallEventYr,
     &                 y2Rainfall, numeroYr, arrayRP(m), outRainfall )

          call splint( SortAscReturPeriod, MaxRunoffEventYr,
     &                 y2Runoff, numeroYr, arrayRP(m), outRunoff )

          call splint( SortAscReturPeriod, MaxSoilLossEventYr,
     &                  y2SoilLoss, numeroYr, arrayRP(m), outSoilLoss )
          
          call splint( SortAscReturPeriod, MaxDepoEventYr,
     &                  y2Depo, numeroYr, arrayRP(m), outDepo )
          
          call splint( SortAscReturPeriod, MaxSLossEventYr,
     &                  y2SLoss, numeroYr, arrayRP(m), outSLoss )
          
          call splint( SortAscReturPeriod, MaxSYEventYr,
     &                  y2SY, numeroYr, arrayRP(m), outSY )          
          
          yRainfall (m) = ( outRainfall / Area ) * 1000. ! (mm)
          yRunoff (m) = ( outRunoff / Area ) * 1000. ! (mm)
          ySoilLoss (m) = outSoilLoss ! (kg/Ha)
          yDepo (m) = outDepo
          ySLoss (m) = outSLoss
          ySY (m) =  outSY
      
      end do

C     + + + Write Summary Output + + +
C---  BEGIN Mariano Hernandez 03May2017 ---
C---  Change from ton to Mg in Avg-Soil-Loss(ton/ha/year) to (Mg/ha/year)
C---  Change from ton to Mg in Avg-SY(ton/ha/year) to (Mg/ha/year)

      write(file_unit, 12) "-ANNUAL-AVERAGES-"
C     - - -  Annual Averages - - -
      write(file_unit, 15) "Avg. Precipitation (mm/year) =", 
     &                                        AnnualRainfallAverage
      write(file_unit, 15) "Avg-Runoff(mm/year)= ", 
     &                                          AnnualRunoffAverage
      write(file_unit, 15) "Avg-Soil-Loss(Mg/ha/year)=",
     &                                    AnnualSLossAverage/1000.
      
      write(file_unit, 15) "Avg-SY(Mg/ha/year)=", 
     &                                    AnnualSYAverage/1000.
      
C--- BEGIN Mariano Hernandez 03May2017
C--- Write output for RETURN-FREQUENCY-RESULTS-YEARLY TOTALS
      
      write(file_unit,11)
      write(file_unit, 12) "RETURN-FREQUENCY-RESULTS-YEARLY TOTALS"
      write(file_unit, 14) "Variable        2-yr            5-yr
     &            10-yr           25-yr           50-yr       100-yr" 
C     - - - Return Frequency Results - - - 
      write(file_unit, 20) "Rain (mm)",(yRainfallAnnual(i), i = 1, nRP )
      write(file_unit, 20) "Runoff(mm)",( yRunoffAnnual(i), i = 1, nRP )
      write(file_unit, 20) "Soil-Loss(Mg/ha)",
     &                             ( ySLossAnnual(i)/1000., i = 1, nRP )
      write(file_unit, 20) "Sediment-Yield(Mg/ha)",
     &                            ( ySYAnnual(i)/1000., i = 1, nRP )
          write(file_unit, 11)     
C---END Output for RETURN-FREQUENCY-RESULTS-YEARLY TOTALS     
C---BEGIN Mariano Hernandez 03May2017
C---Write output for RETURN_FREQUENCY-RESULTS-YEARLY MAXIMUM DAILY
          
      write(file_unit,11)
      write(file_unit,12)"RETURN-FREQUENCY-RESULTS-YEARLY MAXIMUM DAILY"
      write(file_unit, 14) "Variable        2-yr            5-yr
     &            10-yr           25-yr           50-yr       100-yr" 
C     - - - Return Frequency Results - - - 
          write(file_unit, 20) "Rain (mm)", ( yRainfall(i), i = 1, nRP )
          write(file_unit, 20) "Runoff(mm)", ( yRunoff(i), i = 1, nRP )
          write(file_unit, 20) "Soil-Loss(Mg/ha)",
     &                               ( ySLoss(i)/1000., i = 1, nRP )
          write(file_unit, 20) "Sediment-Yield(Mg/ha)",
     &                               ( ySoilLoss(i)/1000., i = 1, nRP )
          write(file_unit, 11)
C --- Risk Assessment Analysis ------------------------------------------------
C          write(file_unit, 12) "-RISK-ASSESSMENT-ANALYSIS-"
C          write(file_unit, 16) " Severity Description - Total Annual
C     & Soil Loss Tolerance - Probability -"
C          write(file_unit, 17) "(ton/ha)"
C          write(file_unit, 21) "Negligible                   Up to 0.5",
C     &                                                    prob_up205/1.
C          write(file_unit, 21) "Acceptable                   0.5 - 1.0",
C     &                                                    prob_05_10/1.
C          write(file_unit, 21) "Undesirable                  1.0 - 3.0",
C     &                                                    prob_10_30/1.
C          write(file_unit, 21) "Unacceptable                 Over  3.0",
C     &                                                    prob_over3/1.
C          
C          write(file_unit, 21) "Negligible                 Up to 0.25",
C     &                                                    prob_up205/1.
C          write(file_unit, 21) "Acceptable                0.25 - 0.50",
C     &                                                    prob_05_10/1.
C          write(file_unit, 21) "Undesirable               0.50 - 1.50",
C     &                                                    prob_10_30/1.
C          write(file_unit, 21) "Unacceptable               Over  1.50",
C     &                                                    prob_over3/1.          
C -----------------------------------------------------------------------------          
          
C10    format(4i5, 6F15.5)
11    format(/)
12    format(5x,a/)
14    format(2x, a/) 
15    format(2x, a, 2x, F15.5)
16    format(2x, a)
17    format(40x, a)      
20    format(2x, a, 6F15.5)
21    format(10x,a,  F20.5)      
      end
      
      SUBROUTINE sort(n,arr)  
      INTEGER n,M,NSTACK  
      REAL arr(n)  
      PARAMETER (M=7,NSTACK=50)  
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)  
      REAL a,temp  
      jstack=0  
      l=1  
      ir=n  
1     if(ir-l.lt.M)then  
        do 12 j=l+1,ir  
          a=arr(j)  
          do 11 i=j-1,l,-1  
            if(arr(i).le.a)goto 2  
            arr(i+1)=arr(i)  
11        continue  
          i=l-1  
2         arr(i+1)=a  
12      continue  
        if(jstack.eq.0)return  
        ir=istack(jstack)  
        l=istack(jstack-1)  
        jstack=jstack-2  
      else  
        k=(l+ir)/2  
        temp=arr(k)  
        arr(k)=arr(l+1)  
        arr(l+1)=temp  
        if(arr(l).gt.arr(ir))then  
          temp=arr(l)  
          arr(l)=arr(ir)  
          arr(ir)=temp  
        endif  
        if(arr(l+1).gt.arr(ir))then  
          temp=arr(l+1)  
          arr(l+1)=arr(ir)  
          arr(ir)=temp  
        endif  
        if(arr(l).gt.arr(l+1))then  
          temp=arr(l)  
          arr(l)=arr(l+1)  
          arr(l+1)=temp  
        endif  
        i=l+1  
        j=ir  
        a=arr(l+1)  
3       continue  
          i=i+1  
        if(arr(i).lt.a)goto 3  
4       continue  
          j=j-1  
        if(arr(j).gt.a)goto 4  
        if(j.lt.i)goto 5  
        temp=arr(i)  
        arr(i)=arr(j)  
        arr(j)=temp  
        goto 3  
5       arr(l+1)=arr(j)  
        arr(j)=a  
        jstack=jstack+2  
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'  
        if(ir-i+1.ge.j-l)then  
          istack(jstack)=ir  
          istack(jstack-1)=i  
          ir=j-1  
        else  
          istack(jstack)=j-1  
          istack(jstack-1)=l  
          l=i  
        endif  
      endif  
      goto 1  
      END
      
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)  
      INTEGER n,NMAX  
      REAL yp1,ypn,x(n),y(n),y2(n)  
      PARAMETER (NMAX=500)  
      INTEGER i,k  
      REAL p,qn,sig,un,u(NMAX)  
      if (yp1.gt..99e30) then  
        y2(1)=0.  
        u(1)=0.  
      else  
        y2(1)=-0.5  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
      endif  
      do 11 i=2,n-1  
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
        p=sig*y2(i-1)+2.  
        y2(i)=(sig-1.)/p  
        u(i)=(6.*((y(i+1)-y(i))/(x(i+  
     &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*  
     &u(i-1))/p  
11    continue  
      if (ypn.gt..99e30) then  
        qn=0.  
        un=0.  
      else  
        qn=0.5  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
      endif  
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)  
12    continue  
      return  
      END
      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)  
      INTEGER n  
      REAL x,y,xa(n),y2a(n),ya(n)  
      INTEGER k,khi,klo  
      REAL a,b,h  
      klo=1  
      khi=n  
1     if (khi-klo.gt.1) then  
        k=(khi+klo)/2  
        if(xa(k).gt.x)then  
          khi=k  
        else  
          klo=k  
        endif  
      goto 1  
      endif  
      h=xa(khi)-xa(klo)  
      if (h.eq.0.) pause 'bad xa input in splint'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**  
     &2)/6.  
      return  
      END
      