
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo76.F90
    !> \brief   Spot test - 2500K with 1.0 Cs - 0.5 Te
    !> \author  M.H.A. Piro, B.A.T. Breeden
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
    !    02/24/2022    B.A.T. Breeden      SUBI Test Case
    !    17/04/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the open literature Cs-Te assessment file at 2500K with 1.0 mol of Cs and 0.5 mol of Te.
    !!  It also tests mixing term Case #9 of the SUBI phase.
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: T. N. Pham Thi, J. C. Dumas, V. Bouineau, N. Dupin, C. Gueneau, S. Gosse, 
    !!  P. Benigni, P. Maugis and J. Rogez, "Thermodynamic assessment of the Cs–Te binary system," Calphad,
    !!  vol. 48, pp. 1-12, 2015.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo76

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4
    real(8) :: pcheck1, pcheck2, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'CsTe-2.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(52)        = 0.5D0          ! Te
    dElementMass(55)        = 1.0D0          ! Cs

    ! Liquid #1
    sfcheck1 = 1D0               !Cs
    sfcheck2 = 0.34222D0         !Va
    sfcheck3 = 0.48667D0         !Cs2Te
    sfcheck4 = 0.17111D0         !Te

    pcheck1 = -377921D0        !Cs
    pcheck2 = -381967D0        !Te

    gibbscheck = -568905D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    s1pass = .FALSE.
    s2pass = .FALSE.
    s3pass = .FALSE.
    s4pass = .FALSE.

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS((dGibbsEnergySys - (gibbscheck))/(gibbscheck)) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3)) then
            do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'IONIC_LIQUID') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Cs+') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Cs2Te') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Te') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck4)/sfcheck4 < 1D-3) s4pass = .TRUE.
                            end if
                        end do
                    end do

                end if
            end do
        end if
    end if

    if (s1pass .AND. &
        s2pass .AND. &
        s3pass .AND. &
        s4pass) then
        ! The test passed:
        print *, 'TestThermo76: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo76: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo76
