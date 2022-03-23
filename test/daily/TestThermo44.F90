
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo44.F90
    !> \brief   Spot test - 1000K with 40% Pd, 60% Tc.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
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
    !    05/14/2013    M.H.A. Piro         Original code
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use Kaye's Pd-Ru-Tc-Mo system
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for the Pd-Ru-Tc-Mo system at 1000K with 40% Pd, 60% Tc.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo44

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: i,j,k
    logical :: s1pass, s2pass, cppass

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000D0
    dElementMass(46)       = 40D0        ! Pd
    dElementMass(43)       = 60D0        ! Tc

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    s1pass = .FALSE.
    s2pass = .FALSE.
    cppass = .FALSE.
    ! Check results:
    if (INFOThermo == 0) then
        if (DABS(dGibbsEnergySys - (-5.04309D6))/(-5.04309D6) < 1D-3) then
            do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'HCPN') then
                    do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                        if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                            if (DABS(dMolFraction(j) - 0.4D0)/0.4D0 < 1D-3) s1pass = .TRUE.
                        else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                            if (DABS(dMolFraction(j) - 0.6D0)/0.6D0 < 1D-3) s2pass = .TRUE.
                        end if
                    end do
                end if
            end do
            if (ABS(dHeatCapacity - 3023.31)/3023.31 < 1D-3) cppass = .TRUE.
        end if
    end if

    if (s1pass .AND. s2pass .AND. cppass) then
        ! The test passed:
        print *, 'TestThermo44: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo44: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo44
