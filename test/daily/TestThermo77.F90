
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo77.F90
    !> \brief   Spot test - 2500K with 1.0 Nb - 0.7 Sn - 0.3 O.
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
    !!  results for the Zirc Data file at 2500K with 1.0 mol of Nb, 0.7 mols of Sn, and 0.3 of O. It also
    !!  tests mixing term Case #7 of the SUBI phase, with the presence of a miscibility gap. Permission
    !!  was granted from N. Dupin to make use of the Zirc DAT file. The data file for this test case has been
    !!  modified from it's original state.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo77

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4
    real(8) :: sfcheck5, sfcheck6, sfcheck7
    real(8) :: pcheck1, pcheck2, pcheck3, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass, s5pass
    logical :: s6pass, s7pass, s8pass, s9pass, s10pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq_mod1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(41)        = 1D0          ! Nb
    dElementMass(8)         = 0.3D0        ! O
    dElementMass(50)        = 0.7D0        ! Sn

    ! Liquid #1
    sfcheck1 = 0.99586D0          !Nb+2
    sfcheck2 = 0.327163D0         !O-2
    sfcheck3 = 0.563388D0         !Va
    sfcheck4 = 0.109449D0         !NBO5/2

    ! Liquid #2
    sfcheck5 = 0.418758D0         !Nb+2
    sfcheck6 = 0.581242D0         !Sn+2
    sfcheck7 = 0.99976D0          !Va

    pcheck1 = -183696D0       !Nb
    pcheck2 = -593052D0       !O
    pcheck3 = -249101D0       !Sn

    gibbscheck = -535982D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    s1pass = .FALSE.
    s2pass = .FALSE.
    s3pass = .FALSE.
    s4pass = .FALSE.
    s5pass = .FALSE.
    s6pass = .FALSE.
    s7pass = .FALSE.
    s8pass = .FALSE.
    s9pass = .FALSE.
    s10pass = .FALSE.

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS((dGibbsEnergySys - (gibbscheck))/(gibbscheck)) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3).AND. &
            (DABS((dElementPotential(3)*dIdealConstant*dTemperature - pcheck3)/pcheck3) < 1D-3)) then
            loop_checkPhases: do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'IONIC_LIQ') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'NB+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'NBO5/2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck4)/sfcheck4 < 1D-3) s4pass = .TRUE.
                            end if
                        end do
                    end do
                
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'NB+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck5)/sfcheck5 < 1D-3) s5pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'SN+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck6)/sfcheck6 < 1D-3) s6pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck7)/sfcheck7 < 1D-3) s7pass = .TRUE.
                            end if
                        end do
                    end do
                end if
            end do loop_checkPhases
        end if
    end if

    if (s1pass .AND. &
        s2pass .AND. &
        s3pass .AND. &
        s4pass .AND. &
        s5pass .AND. &
        s6pass .AND. &
        s7pass) then
        ! The test passed:
        print *, 'TestThermo77: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo77: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo77
