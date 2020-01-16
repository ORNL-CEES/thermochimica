program msre

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer :: i, nCalc, nFail, iOutFrequency
    real(8)   :: dTotalMass, dmLi, dmBe, dmF, dmTh, dmU, dmPu
    real(8)   :: minTemp, maxTemp, difTemp, rTemp
    real(8)   :: minPres, maxPres, difPres, rPres
    character(120)         :: cPassOutFileName, cFailOutFileName

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'LiBeUPuThF_T1.dat'
    ! cThermoFileName       = DATA_DIRECTORY // 'MSAR_Ver_1.dat'

    ! Specify values:
    minTemp = 800D0
    maxTemp = 1200D0
    difTemp = maxTemp - minTemp
    minPres = 1D0
    maxPres = 1D0
    difPres = maxPres - minPres
    nCalc = 100000
    cPassOutFileName = 'carpetMSREPass.txt'
    cFailOutFileName = 'carpetMSREFail.txt'

    ! Specify output and debug modes:
    iPrintResultsMode     = 0
    lDebugMode            = .FALSE.
    iOutFrequency         = 500

    nFail = 0
    i = 0
    LOOP_carpet: do while (i < nCalc)

      i = i + 1

      rTemp = RAND(0)
      rPres = RAND(0)
      dTemperature          = minTemp + difTemp * rTemp
      dPressure             = minPres + difPres * rPres
      dElementMass          = 0D0

      dElementMass(3)  = 55.5D0 + 20D0*RAND(0)  ! Li nominal 65%
      dElementMass(4)  = 26.2D0 + 6D0*RAND(0)  ! Be nominal 29.2%
      dElementMass(90) = 5.02D0  + RAND(0)  ! Th (Zr) nominal 5.52%
      dElementMass(92) = 0.73D0 + 0.2D0*RAND(0)  ! U nominal 0.83%
      if (RAND(0) > 0.5D0) dElementMass(94) = 0.01D0*RAND(0)  ! Pu nominal 0%
      dElementMass(9) = (dElementMass(3) + 2D0*dElementMass(4) + 4D0*dElementMass(90) + (3.0D0 + RAND(0))*dElementMass(92) &
                         + 3D0*dElementMass(94))     ! F


      ! if (i < 40000) cycle LOOP_carpet
      dTotalMass = dElementMass(3)+dElementMass(4)+dElementMass(9)+dElementMass(90)+dElementMass(92)+dElementMass(94)

      dElementMass(3)  = dElementMass(3)  / dTotalMass
      dElementMass(4)  = dElementMass(4)  / dTotalMass
      dElementMass(9)  = dElementMass(9)  / dTotalMass
      dElementMass(90) = dElementMass(90) / dTotalMass
      dElementMass(92) = dElementMass(92) / dTotalMass
      dElementMass(94) = dElementMass(94) / dTotalMass

      call ParseCSDataFile(cThermoFileName)
      call Thermochimica

      dTotalMass = dElementMass(3)+dElementMass(4)+dElementMass(9)+dElementMass(90)+dElementMass(92)+dElementMass(94)

      dmLi = dElementMass(3)  / dTotalMass
      dmBe = dElementMass(4)  / dTotalMass
      dmF  = dElementMass(9)  / dTotalMass
      dmTh = dElementMass(90) / dTotalMass
      dmU  = dElementMass(92) / dTotalMass
      dmPu = dElementMass(94) / dTotalMass

      OPEN(23, FILE = cPassOutFileName)
      OPEN(24, FILE = cFailOutFileName)
      if (INFOThermo == 99) then
          i = i - 1
      else if (INFOThermo > 0) then
        nFail = nFail + 1
        call ThermoDebug
        WRITE(24,*) 'FAILED: INFO ', INFOThermo
        WRITE(24,*) 'Temperature, pressure, masses: ', dTemperature, dPressure, dmLi, &
                    dmBe, dmF, dmTh, dMu, dmPu
      else
        WRITE(23,*) 'Temperature, pressure, masses: ', dTemperature, dPressure, dmLi, &
                    dmBe, dmF, dmTh, dMu, dmPu
      end if

      if (iPrintResultsMode > 0)  call PrintResults

      ! Destruct everything:
      call ResetThermoAll

      if (MODULO(i,iOutFrequency) == 0) print *, i, "tests complete"

    end do LOOP_carpet

    print *, nFail, ' tests failed out of ', nCalc

end program msre
