! Front end program unit to run K2 in batch mode from AGWA 1.5
!
! kin.fil file is a series of lines, one line for each simulation, with the following comma-separated values:
!
! parameter-file, rainfall-file, output-file, title, tfin, delt, cour?, sed?, multiplier-file?, table?
!
! The title must be enclosed in double quotes.
!
! Entries after delt are optional, will default to 'N'


  program AGWAshell
  
  use runpars
  use elpars
  use multip
  use kinsed_def

  integer             :: file_error, rfile, file_unit,number_of_arguments
  integer(4)          :: arg_len
  logical             :: char2log
  real                :: x, y, t(1000), z(1000), sat, sumbal(10)
  real                :: outRainfall, outInfil,stmdur
  real                :: outRunoff, outPeakR, outsusp
  real                :: outSoilLoss,outsloss,outdpo, outsy, Area
  real                :: inRainfall(110000), inInfil(110000)
  real                :: inRunoff(110000), inPeakR(110000)
  real                :: inSoilLoss(110000), inDpo(110000)
  real                :: inSLoss(110000), inSY(110000), insusp(110000)
  real                :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9
  real                :: time_begin, time_end, mxint
  integer             :: inEvent(110000), inDay(110000)
  integer             :: inMonth(110000), inYear(110000)
  character(len=2)    :: flag,flag2print
  character(len=20)   :: rmks_char, rmn_char, rmcv_char, rmg_char
  character(len=20)   :: rmin_char, rmcoh_char, rmspl_char
  character(len=50)   :: version
  character(len=1000) :: line
  character(len=1)    :: bell, c
  character(len=10)   :: block_name
  character(len=150)  :: runfile, auxfile, fields(11), file_name

10 format(a) ! for general character i/o
      
  version = 'Compiled 07-24-2012 with Intel Fortran 11.1'
  bell    = char(7)
  call cpu_time( time_begin )
! Default values

  cour = .false.
  sed  = .false.
  tabl = .false.
  suma = .false.
!
!  number_of_arguments = nargs( )
!  call getarg( 1, flag, arg_len )
  if (COMMAND_ARGUMENT_COUNT () /= 0 ) then
       call GET_COMMAND_ARGUMENT ( 1, flag, arg_len )
       if (flag == '-b') then
!      if(number_of_arguments == 3)then
          if (COMMAND_ARGUMENT_COUNT () == 2 ) then
           call GET_COMMAND_ARGUMENT ( 2, runfile, arg_len )
!          call getarg( 2, runfile, arg_len )
       else
          runfile='kin.fil'
          arg_len = 7
         end if
       end if
  end if
!  
!  
  !if(nargs() == 1) then
  !  runfile = 'kin.fil'
  !  arg_len = 7
  !else
  !  call getarg( 1, runfile, arg_len )
  !endif
!
  open(4, file = runfile(1:arg_len), status = 'old', iostat = file_error)
  if(file_error .ne. 0) call errxit('AGWAshell', "Can't open kin.fil")

! Begin simulation loop

   do p = 1, 1

    read(4, 10, iostat = file_error) line

    if(file_error .ne. 0) exit

!   Replace double quotes around title field with spaces

    i = index(line, '"')
    line(i:i) = ' '
    j = index(line, '"')
    line(j:j) = ' '

!   Temporarily replace any commas in title field with bell characters
   
    k = index(line(i:j), ',')
    do while( k > 0)
      k = i + k - 1
      line(k:k) = bell
      k = index(line(i:j), ',')
    end do

!   Identify comma-separated fields
   
    fields = ' '
    i = 1
    j = 1
    k = index(line, ',')
    do while( k > 0)
      line(k:k) = ' '
      fields(i) = line(j:k-1)
      i = i + 1
      j = k + 1
      k = index(line, ',')
    end do
    fields(i) = line(j:)
    flag2print = fields(11) 

!   Restore commas to title field
   
    k = index(fields(4), bell)
    do while( k > 0)
      fields(4)(k:k) = ','
      k = index(fields(4), bell)
    end do

!   Open parameter, rainfall and output files - the rainfall file is treated
!   as optional (the file unit is set to -1 as a flag for subroutine rain to
!   create a single raingage entry with zero depths)

    open(files(1), file = trim(fields(1)), status = 'old', iostat = file_error)

    if(file_error .ne. 0) call errxit('AGWAshell', "Can't open "//trim(fields(1)))

    if(len_trim(fields(2)) == 0) then ! no rainfall file specified     
      rfile  = -1
      fields(2) = 'None'
    else
      open(files(2), file = fields(2), status = 'old', iostat = file_error)
      if(file_error .ne. 0) call errxit('AGWAshell', "Can't open "//fields(2))
      rfile  = files(2)
    end if
  
    open(files(3), file = fields(3), status = 'unknown', iostat = file_error)
    if(file_error .ne. 0) call errxit('AGWAshell', "Can't open "//fields(3))

!   Decode character input into floating and logical values (module runpars)

!    read(fields(5) , *) tfin
    read(fields(6) , *) delt
   
!   If no value for a logical variable is specified, it remains at its default value
   
    if(len_trim(fields(7))  > 0) cour = char2log(fields(7))
    if(len_trim(fields(8))  > 0) sed  = char2log(fields(8))
    if(len_trim(fields(10)) > 0) tabl = char2log(fields(10))
    if(len_trim(fields(11)) > 0) suma = char2log(fields(11))
    
!   Two options for parameter multipliers:

!   1) Use file specified in kin.fil
!   2) No multipliers 

    if(len_trim(fields(9))      == 0 .or. &
         verify(fields(9),'N ') == 0 .or. &
         verify(fields(9),'n ') == 0) then ! no multipliers

      rmks_char  = '1.0'
      rmn_char   = '1.0'
      rmcv_char  = '1.0'
      rmg_char   = '1.0'
      rmin_char  = '1.0'
      rmcoh_char = '1.0'
      rmspl_char = '1.0'
  
      fields(9) = 'N'
  
    else
  
      open(7, file = fields(9), status = 'old', iostat = file_error)
      if(file_error .ne. 0) call errxit('AGWAshell', "Can't open "//fields(9))
  
      read(7, 10) rmks_char
      read(7, 10) rmn_char
      read(7, 10) rmcv_char
      read(7, 10) rmg_char
      read(7, 10) rmin_char
      read(7, 10) rmcoh_char
      read(7, 10) rmspl_char
      close(7)
  
    end if
  
!   Decode multipliers into floating values (module multip)
   
    read(rmks_char,  *) rmks
    read(rmn_char,   *) rmn
    read(rmcv_char,  *) rmcv
    read(rmg_char,   *) rmg
    read(rmin_char,  *) rmin
    read(rmcoh_char, *) rmcoh
    read(rmspl_char, *) rmspl
    
    if( len_trim(fields(11))     == 0 .or. &
        verify(fields(11), 'N ') == 0 .or. &
        verify(fields(11), 'n ') == 0 ) then ! No Summary Large File
        fields(11) = 'N'
   else
       k = index(fields(3), '.')
       open(23, file = fields(3)(1:k)//'out', status='unknown', iostat = file_error)
       if(file_error .ne. 0) call errxit('DRHEM', "Can't open "//fields(3)(1:k)//'out')
       
!   Write info to output file
    write (23, 1957)'Event-Id', 'Day', 'Month', 'Year', 'P-rain', &
      'Max.P.Int', 'P-dur', 'Q-runoff', 'Q-Peak', 'Q-dur',        &
      'Sed-Yield', 'Soil-Loss'
    
    write (23, 1958)'#', '#', '#','#', 'mm', 'mm/hr', 'min', 'mm', 'mm/hr', 'min',    &
      'ton/ha', 'ton/ha'
   end if
   
       
   
!   Write info to output file
!    write (23, 1957)'Event-Id', 'Day', 'Month', 'Year', 'P(rain)', &
!      'Max P Int', 'P-dur', 'Q(runoff)', 'Q-Peak', 'Q-dur',        &
!      'Sed Yield', 'Soil Loss'
    
!    write (23, 1958) 'mm', 'mm/hr', 'min', 'mm', 'mm/hr', 'min',    &
!      'kg/ha', 'kg/ha'
    
1957 format(1x, a8, 3x, a3, 3x, a5, 3x, a4, 4x, a7, 4x, a10,       &
            4x, a5, 4x, a9, 4x, a6, 6x, a5, 7x, a9, 3x, a9)
     
!1958 format(31x, a7, 4x, a10, 6x, a5, 3x, a9, 7x, a6,              &
!             5x, a5, 5x, a9, 3x, a9)    
     
1958 format(1x, a7, 3x, a3, 3x, a5, 3x, a4, 4x, a7, 4x, a10, 6x, a5, 3x, a9, 7x, a6,              &
             5x, a5, 5x, a9, 3x, a9)
!    write(files(3), 30) version, trim(fields(4)), trim(fields(1)), &
!                        trim(fields(2)), (trim(fields(i)), i = 5, 10)

!30  format(' KINEROS2 - ',a//           &
!           ' Title: ',a//                     &
!           ' Parameter File Used....... ',a/  &
!           ' Rainfall File Used........ ',a// &
!           ' Length of Run, minutes.... ',a/  &
!           ' Time Step, minutes........ ',a/  &
!           ' Use Courant criteria?..... ',a/  &
!           ' Simulate Sed. Transport?.. ',a/  &
!           ' Multiplier file (if any).. ',a/  &
!           ' Tabular Summary?.......... ',a)
!
!    write(files(3), 40) rmks_char, rmn_char, rmcv_char, rmg_char, &
!                        rmin_char, rmcoh_char, rmspl_char
!
!40  format(/' Multipliers Used:'//           &
!            ' Saturated Conductivity... ',a/ &
!            ' Manning n................ ',a/ &
!           ' CV of Ksat............... ',a/  &
!           ' Capillary Drive Coeff.... ',a/  &
!           ' Intercepted Depth........ ',a/  &
!           ' Sediment Cohesion Coeff.. ',a/  &
!           ' Sediment Splash Coeff.... ',a/)

!   Initialize global volumes
    
!    do j = 1, 10
!    sumbal(j) = 0.
!    end do
    
!   Compute number of time steps
    
!    limit = int (tfin / delt) + 1
    
!   Convert time increment to seconds
    
    delt = delt * 60.
    
!   Read global parameter block
    
    call reader(files(1), block_name, ierr)
    if(ierr .gt. 0) call errxit('AGWAshell', "Invalid global block")
    
    call getr4('C', 0, clen, ierr)
    if(ierr .ne. 0) call errxit('AGWAshell', "char. length not found")
    
!   System of units
    
    call getstr ('U', 0, c, l, ierr)
    if(ierr .ne. 0) call errxit('AGWAshell', "units not specified")
    
    if(c .eq. 'm' .or. c .eq. 'M') then ! metric
      units  = 1
      dlab   = 'm.'
      dilb   = 'mm'
      conv   = 1000.
      wt     = 1000.  ! kg/m3 water density
      arconv = 10000. ! m^2 per ha.
      wconv  = 1000.  ! kg/tonne
      grav   = 9.81
      bdep   = 656.
    else if(c .eq. 'e' .or. c .eq. 'E') then ! english
      units  = 2
      dlab   = 'ft'
      dilb   = 'in'
      conv   = 12.
      wt     = 62.4    
      arconv = 43560.
      wconv  = 2000.
      grav   = 32.2
      bdep   = 200.
    else
      call errxit('K2shell', "units not recognized")
    end if
    
!    write (*, '(/)')
    
!   Rain on channel?
    
    call getstr( 'CHR', 0, c, l, ierr)
    if(ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then
      chrain = .true.
    else
      chrain = .false. ! default
    end if
    
!   Need temperature (C) to compute kinematic viscosity for laminar flow option in plane
    
    call getr4('T', 0, temp, ierr)
    if(ierr .gt. 0) temp = 0.
    if(units .eq. 2) temp = (temp - 32.) * 5. / 9.
    xnu = .0000017756 / (1.+ 0.03368 * temp + 0.000221 * temp*temp) ! m**2/s
    if(units .eq. 2) xnu = xnu * 10.76496 ! convert to ft**2/s
    
! ----- CLIGEN- Mariano Hernandez 12Nov2012 -----------------------------------
! Obtain the number of events produced by CLIGEN
    call eatcom ( rfile )
    read ( rfile, * ) nevent
    read ( rfile, * ) ibrkpt
! -----------------------------------------------------------------------------
        
    if(sed) call sed00()
    
    call wrt00(files(3))
    
    call reader( files(1), block_name, ierr )
          
    
    if(ierr .gt. 0) then
        if(ierr .eq. 1) exit ! end of file
        call errxit('AGWAshell', "error reading parameter file")                                                                   
    end if
    
    do m = 1, nevent
        
        do j = 1, 10
           sumbal(j) = 0.
        end do
        
        do j = 1, 1000
          time_hyd_min(j) = 0.
          time_hyd_max(j) = 0.
        end do
            
        
        inEvent(m) = 0
        inDay(m) = 0
        inMonth(m) = 0
        inYear(m) = 0
        inRainfall(m) = 0.
        inInfil(m) = 0.
        inRunoff(m) = 0.
        inPeakR(m) = 0.
        inSoilLoss(m) = 0.
        inSLoss(m) = 0.
        indpo(m) = 0.
        inSY(m) = 0.
        insusp(m) = 0.
        
    call rain( rfile, idevent, iday, imon, iyear, tfin, stmdur, prcp, mxint )
    
!   Compute number of time steps

    limit = int (tfin / delt) + 1
    
!   Convert time increment to seconds
    
!    delt = delt * 60.
    
    call clerk(limit)
    
    
!    call rain( rfile, idevent, iday, imon, iyear, tfin )
    
    
    do j = 1, 4
        
      dtm( j ) = delt
      
    end do
    
    nelt = 0
    ltab = 0
    
!   Element processing loop...
    
!    do
!      call reader( files(1), block_name, ierr )
          
    
!      if(ierr .gt. 0) then
!        if(ierr .eq. 1) exit ! end of file
!        call errxit('AGWAshell', "error reading parameter file")                                                                   
!      end if
     
      nprint = 0
      nchan  = 1
!      nelt   = nelt + 1 ! these will not agree when there is a
!      ltab   = ltab + 1 ! compound channel with overbank infiltration
    
!     Print = 3 option
    
      call geti4 ('PR', 0, nprint, ierr)
    
      if(nprint .eq. 3) then
        call getstr ('FI', 0, auxfile, j, ierr)
        if (ierr .gt. 0) then
          call errxit('AGWAshell', "PRINT=3, but no file specified")                                                                   
        else
          open (files(4), file = auxfile, status = 'unknown', iostat = file_error)
          if (file_error .ne. 0) call errxit('K2shell', "can't open "//auxfile)                                                                   
        end if
      end if
    
!     Report progress
    
!      write (*, "('+ Processing ', a8, i8)") block_name, id
    
      if(block_name(1:5) .eq. 'PLANE'  .or.                &
        (block_name(1:7) .eq. 'CHANNEL' .and. chrain) .or. &
         block_name(1:4) .eq. 'POND' .or.                  &
         block_name(1:5) .eq. 'URBAN') then ! Rainfall
      
        call geti4 ('RA', 0, ngage, ierr)
        if(ierr .gt. 0) then
          ngage = 0
          call getr4 ('X', 0, x, ierr)
          if(ierr .gt. 0) x = 0.
          call getr4 ('Y', 0, y, ierr)
          if(ierr .gt. 0) y = 0.
        end if
        call interp (x, y, t, z, nd, sat, ngage)
      end if
    
      if(block_name(1:5) .eq. 'PLANE') then  
        call plane(t, z, nd, sat)
!      else if(block_name(1:7) .eq. 'CHANNEL') then
!        call channl(t, z, nd, sat)
!      else if(block_name(1:4) .eq. 'PIPE') then
!        call pipe()                      
!      else if(block_name(1:6) .eq. 'INJECT') then
!        call inject()
!      else if(block_name(1:4) .eq. 'POND') then
!        call pond(t, z, nd)
!      else if(block_name(1:5) .eq. 'URBAN') then                       
!        call urban(t, z, nd, sat)
!      else if(block_name(1:5) .eq. 'ADDER') then                       
!        call adder()
      else
        call errxit(block_name, 'invalid element')
      end if
      
!     Add element contribution to global volume balance componets
    
      do i = 1, 9
        sumbal(i) = sumbal(i) + qbal(i)
      end do
    
      call writer( files(4), idEvent, iDay, iMon, iYear,   &
      outRainfall, outInfil, outRunoff, outPeakR,          &
      outSoilLoss, outsloss, outdpo, outsusp, stmdur, prcp, mxint )
!   End of element processing loop
    
!    end do
    
!   Outflow from final element
    
    sumbal(10) = qbal(10)
    Area = qbal(1)
    inEvent(m) = idEvent 
    inDay(m) = iDay
    inMonth(m) = iMon
    inYear(m) = iYear
!    inRainfall(m) = outRainfall
    inRainfall(m) = qbal(2)
    inInfil(m) = outInfil
!    inRunoff(m) = outRunoff
    inRunoff(m) = qbal(10) !m^3
    inPeakR(m) = outPeakR
    inSoilLoss(m) = outSoilLoss
    inSLoss(m) = outsloss
    indpo(m) = outdpo
    insy(m) = outSoilLoss
    insusp(m) = outsusp
    
!   Write  event summary
    
!    call event( sumbal, idevent, iday, imon, iyear )

    end do
    
!    do i = 1, 4
!     close( files(i) )
!    end do
!    close( files(1) )
!    close( files(2) )
!    close( files(3) )
!    close( files(4) )

!   Write spreadsheet file for final element - name is same as output file but with 'csv' extension

!    k = index(fields(3), '.')
!    open(files(4), file = fields(3)(1:k)//'csv', status = 'unknown', iostat = file_error)
!    if(file_error .ne. 0) call errxit('AGWAshell', "Can't open "//fields(3)(1:k)//'csv')
!    nprint = 3
!    call writer(files(4), idevent, iday, imon, iyear )
!    close(files(4))

! End of simulation loop
    
   end do
   file_unit = files(3)
   file_name = fields(3)
    call Summary_Stats(nevent, inEvent, inDay, inMonth,    &
    inYear, Area, inRainfall, inInfil, inRunoff,           &
    inPeakR, inSoilLoss, insloss*(-1.), indpo, insy, insusp,     &
    file_unit, file_name )
    
    
     call cpu_time( time_end )
!    write (*,*) tt1, tt6
     WRITE (*,*) 'Time of operation was ', time_end - time_begin, ' seconds'
end program AGWAshell


logical function char2log(char)

! Converts character values Y|N or y|n into true|false
! Returns char as uppercase

  character(len=*) :: char

  if(verify(char, 'Y ') == 0 .or. verify(char, 'y ') == 0) then
    char2log = .true.
    char = 'Y'
  else
    char2log = .false.
    char = 'N'
  end if
  return
end

subroutine errxit (id, message)

! Write error message to stdout & output file

   use runpars

   character(len=*) id, message

!   write(files(3), 10) id, trim(message)
!   write(*, 10)        id, trim(message)
!10 format(//, ' error - ', a, /, 1x, a)
   stop ' '
end
