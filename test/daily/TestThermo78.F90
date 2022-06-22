
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo78.F90
    !> \brief   Spot test - 2000K with  0.5 Cr - 0.8 Sn - 0.6 O.
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
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the Zirc Data file at 2000K with 0.5 mols of Cr, 0.8 mols of Sn, and 0.6 mols of O.
    !!  It also tests mixing term Case #12 of the SUBI phase. Permission was granted from N. Dupin to make
    !!  use of the Zirc DAT file. The data file for this test case has been modified from it's original state.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo78

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4
    real(8) :: pcheck1, pcheck2, pcheck3, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq_mod2.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(24)        = 0.5D0        ! Cr
    dElementMass(50)        = 0.8D0        ! Sn
    dElementMass(8)         = 0.6D0        ! O

    ! Liquid #1
    sfcheck1 = 0.50000D0         !Cr+3
    sfcheck2 = 0.50000D0         !Sn+2
    sfcheck3 = 0.76923D0         !Va
    sfcheck4 = 0.23077D0         !O2SN

    pcheck1 = -723575D0       !Cr
    pcheck2 = -467079D0       !O
    pcheck3 = -780227D0       !Sn

    gibbscheck = -1266220D0

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
        if ((DABS(dGibbsEnergySys - (gibbscheck))/(gibbscheck) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3).AND. &
            (DABS((dElementPotential(3)*dIdealConstant*dTemperature - pcheck3)/pcheck3) < 1D-3)) then
            do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'IONIC_LIQ') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'CR+3') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'SN+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O2SN') then
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
        print *, 'TestThermo78: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo78: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo78
