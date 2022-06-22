
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo70.F90
    !> \brief   Spot test - 1500K with 70% Sn, 30% O.
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
    !    01/10/2022    B.A.T. Breeden      Modification to use Dupin's Zirc Data base with SUBI
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the Zirc Data file at 1500K with 0.7 mols of Sn and 0.3 of O. It also tests mixing term Case #2
    !!  an #4 of the SUBI phase, with the presence of a miscibility gap. Permission was granted from N. Dupin
    !!  to make use of the Zirc DAT file.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo70

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4
    real(8) :: sfcheck5, sfcheck6, sfcheck7, sfcheck8
    real(8) :: pcheck1, pcheck2, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass
    logical :: s5pass, s6pass, s7pass, s8pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1500D0
    dElementMass(8)        = 0.3D0          ! O
    dElementMass(50)       = 0.7D0          ! Sn
    ! Liquid #1
    sfcheck1 = 1D0              !Sn
    sfcheck2 = 1.55834D-2       !O-2
    sfcheck3 = 9.84334D-1       !Va
    sfcheck4 = 8.22893D-5       !O2Sn
    ! Liquid #2
    sfcheck5 = 1D0              !Sn
    sfcheck6 = 9.12945D-01      !O-2
    sfcheck7 = 1.35005D-2       !Va
    sfcheck8 = 7.35549D-2       !O2Sn

    pcheck1 = -310505D0        ! O
    pcheck2 = -125800D0        ! Sn
    gibbscheck = -181212D0

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

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (gibbscheck))/(gibbscheck) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3)) then
            !do i = 1, nSolnPhases
            loop_checkPhases: do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'IONIC_LIQ') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'SN+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O2SN') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck4)/sfcheck4 < 1D-3) s4pass = .TRUE.
                            end if
                        end do
                    end do
                
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'SN+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck5)/sfcheck5 < 1D-3) s5pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck6)/sfcheck6 < 1D-3) s6pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck7)/sfcheck7 < 1D-3) s7pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'O2SN') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck8)/sfcheck8 < 1D-3) s8pass = .TRUE.
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
        s7pass .AND. &
        s8pass) then
        ! The test passed:
        print *, 'TestThermo70: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo70: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo70
