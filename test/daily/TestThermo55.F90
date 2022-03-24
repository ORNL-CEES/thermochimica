
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo55.F90
    !> \brief   Spot test - Ni-Cr-Fe-H 300 K.
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
    !    05/21/2020    M. Poschmann         Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including magnetic solution species and excess magnetic terms.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo55

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j, iFirst
    real(8) :: T, B, StructureFactor
    logical :: fccPass, bccPass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC-noSUBI.dat'

    ! Specify values:
    dPressure              = 1000D0
    dTemperature           = 300D0
    dElementMass(23)       = 0.1D0          ! V
    dElementMass(24)       = 1D0            ! Cr
    dElementMass(26)       = 2D0            ! Fe
    dElementMass(28)       = 1D0            ! Ni
    dElementMass(50)       = 0.1D0          ! Sn
    dElementMass(1)        = 1D0            ! H

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    lFuzzyStoich = .FALSE.
    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    fccPass = .FALSE.
    bccPass = .FALSE.
    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dGibbsEnergySys - (-4.96674D04))/((-4.96674D04))) < 1D-3) then
            do i = 1, nSolnPhases
                j = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(j) == 'FCC_A1') then
                    call CompMagneticTemperatureMoment(j,T,B)
                    iFirst = nSpeciesPhase(j-1) + 1
                    StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
                    T = -T * StructureFactor
                    B = -B * StructureFactor
                    if ((DABS((T-97.24D0)/97.24D0) < 1D-3) .AND. (DABS((B-0.13862D0)/0.13862D0) < 1D-3)) then
                        fccPass = .TRUE.
                    end if
                else if (cSolnPhaseName(j) == 'BCC_A2') then
                    call CompMagneticTemperatureMoment(j,T,B)
                    if ((DABS((T-641.38D0)/641.38D0) < 1D-3) .AND. (DABS((B-1.5248D0)/1.5248D0) < 1D-3)) then
                        bccPass = .TRUE.
                    end if
                end if
            end do
        end if
    end if

    if (fccPass .AND. bccPass) then
        ! The test passed:
        print *, 'TestThermo55: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo55: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo55
