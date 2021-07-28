
    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ===========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ===========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    30/11/2018    M. Poschmann        Creation
    !
    ! Purpose:
    ! =========
    ! Executable program for input file form of Thermochimica.
    !
    !-------------------------------------------------------------------------------------------------------------

program ThermochimicaInputScriptMode

  USE ModuleThermoIO
  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none
  character(1024) :: cInputFile
  real(8) :: dTempLow, dTempHigh, dDeltaT, dPressLow, dPressHigh, dDeltaP
  integer :: i, nT, j, nP

  ! Read input argument to get filename
  call get_command_argument(1, cInputFile)
  if (len_trim(cInputFile) == 0) then
    print *, 'No input file specified'
    call EXIT(1)
  endif

  ! Call input parser
  call ParseInput(cInputFile,dTempLow,dTempHigh,dDeltaT,dPressLow,dPressHigh,dDeltaP)
  ! Parse the ChemSage data-file:
  call ParseCSDataFile(cThermoFileName)

  if ((dTempHigh == dTempLow) .OR. (dDeltaT == 0)) then
    nT = 0
  else
    nT = CEILING(DABS(dTempHigh-dTempLow)/DABS(dDeltaT))
  end if

  if ((dPressHigh == dPressLow) .OR. (dDeltaP == 0)) then
    nP = 0
  else
    nP = CEILING(DABS(dPressHigh-dPressLow)/DABS(dDeltaP))
  end if

  if (lStepTogether) then
    nP = 0
    dDeltaP = (dPressHigh - dPressLow) / nT
  end if

  do i = 0, nT
    dTemperature = dTempLow + i*dDeltaT
    if ((dTempHigh > dTempLow) .AND. (dTemperature > dTempHigh)) dTemperature = dTempHigh
    if ((dTempHigh < dTempLow) .AND. (dTemperature < dTempHigh)) dTemperature = dTempHigh

    do j = 0, nP
      if (lStepTogether) then
        dPressure = dPressLow + i*dDeltaP
      else
        dPressure = dPressLow + j*dDeltaP
      end if
      if ((dPressHigh > dPressLow) .AND. (dPressure > dPressHigh)) dPressure = dPressHigh
      if ((dPressHigh < dPressLow) .AND. (dPressure < dPressHigh)) dPressure = dPressHigh

      ! Call Thermochimica:
      call Thermochimica

      if (lReinitRequested) call SaveReinitData

      ! Perform post-processing of results:
      if (iPrintResultsMode > 0)  call PrintResults

      ! Reset Thermochimica:
      call ResetThermo

      ! Call the debugger:
      call ThermoDebug
    end do
  end do

  call ResetThermoAll

end program ThermochimicaInputScriptMode
