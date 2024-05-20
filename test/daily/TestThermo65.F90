
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo65.F90
    !> \brief   Spot test - 400K with 25% Mo, 25% Ru, 25% Pd, 25% Tc.
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
    !    10/01/2021    M. Poschmann         Original code
    !    05/06/2024    A.E.F. Fitzsimmons   Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for the Pd-Ru-Tc-Mo system at 400K with 25% Mo, 25% Ru, 25% Pd, 25% Tc. The database is
    !! modified such that it contains only the BCC phase, but with ternary miscibility gap possible.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo65

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: i,j,k
    logical :: p1pass, p2pass, p3pass

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ternaryMiscibility-Kaye.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 400D0
    dElementMass(42)       = 1D0        ! Mo
    dElementMass(43)       = 1D0        ! Tc
    dElementMass(44)       = 1D0        ! Ru
    dElementMass(46)       = 1D0        ! Pd

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    p1pass = .FALSE.
    p2pass = .FALSE.
    p3pass = .FALSE.
    ! Check results:
    if (INFOThermo == 0) then
        if (DABS((dGibbsEnergySys - (-4.48928D4))/(-4.48928D4)) < 1D-3) then
            do i = 1, nSolnPhases
                j = nElements + 1 - i
                k = -iAssemblage(j)
                if (cSolnPhaseName(k) == 'BCCN') then
                    if (DABS(dMolesPhase(j) - 1.8231D0)/1.8231D0 < 1D-3) p1pass = .TRUE.
                    if (DABS(dMolesPhase(j) - 1.1145D0)/1.1145D0 < 1D-3) p2pass = .TRUE.
                    if (DABS(dMolesPhase(j) - 1.0625D0)/1.0625D0 < 1D-3) p3pass = .TRUE.
                end if
            end do
        end if
    end if

    if (p1pass .AND. p2pass .AND. p3pass) then
        ! The test passed:
        print *, 'TestThermo65: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo65: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo65
