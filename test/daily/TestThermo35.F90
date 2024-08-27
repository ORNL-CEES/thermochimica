
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo56.F90
    !> \brief   Spot test - Fe-Cu-C 1400 K.
    !> \author  M. Poschmann
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    06/30/2020    M. Poschmann        Original code
    !    06/05/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including magnetic solution species and excess magnetic terms. Also tests
    !! G-type asymmetric (and symmetric) excess mixing energy implementation for SUBG.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo35

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j, iFirst
    real(8) :: T, B, StructureFactor
    logical :: fccPass, liquidPass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'CuFeC-Kang.dat'

    ! Specify values:
    dTemperature          = 1400D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(6)       = 1.0D0                              ! C
    dElementMass(26)      = 1.0D0                              ! Fe
    dElementMass(29)      = 1.0D0                              ! Cu

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    fccPass = .FALSE.
    liquidPass = .FALSE.
    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dGibbsEnergySys - (-1.73325D05))/((-1.73325D05))) < 1D-3) then
            do i = 1, nSolnPhases
                j = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(j) == 'FCC_A1') then
                    call CompMagneticTemperatureMoment(j,T,B)
                    iFirst = nSpeciesPhase(j-1) + 1
                    StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
                    T = -T * StructureFactor
                    B = -B * StructureFactor
                    if ((DABS((T-59.4D0)/59.4D0) < 1D-3) .AND. (DABS((B-0.62064D0)/0.62064D0) < 1D-3)) then
                        fccPass = .TRUE.
                    end if
                else if (cSolnPhaseName(j) == 'Liquid') then
                    if (DABS((dMolesPhase(nElements + 1 - i)-2.9666D0)/2.9666D0) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j))-4.4975D-2)/4.4975D-2) < 1D-3) then
                            liquidPass = .TRUE.
                        end if
                    end if
                end if
            end do
        end if
    end if

    if (fccPass .AND. liquidPass) then
        ! The test passed:
        print *, 'TestThermo35: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo35: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo35
