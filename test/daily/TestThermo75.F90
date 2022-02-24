
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo75.F90
    !> \brief   Spot test - 1500K with 0.5 Ca - 0.2 Mn - 0.4 Fe.
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
    !!  results for the open literature Ca-Mn-S assessment file at 1500K with 0.5 mol of Ca, 0.2 mol of Mn, and
    !!  0.4 mol of Fe. It also tests mixing term Case #8 of the SUBI phase.
    !!
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: D. Dilner, "Thermodynamic description of the Fe–Mn–Ca–Mg–S system," Calphad,
    !!  vol. 53, pp. 55-61, 2016.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo75

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
    cThermoFileName        = DATA_DIRECTORY // 'FeMnCaS_mod2.DAT'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1500D0
    dElementMass(26)       = 0.4D0          ! Fe
    dElementMass(25)       = 0.2D0          ! Mn
    dElementMass(20)       = 0.5D0          ! Ca

    ! Liquid #1
    sfcheck1 = 0.45455D0         !Ca+2
    sfcheck2 = 0.18182D0         !Mn+2
    sfcheck3 = 0.36364D0         !Fe+2
    sfcheck4 = 1.0D0             !Va

    pcheck1 = -113841.3D0       !Ca
    pcheck2 = -92132.39D0       !Mn
    pcheck3 = -114102.3D0       !Fe

    gibbscheck = -116594D0

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
                if (cSolnPhaseName(k) == 'IONIC_LIQ#1') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Ca+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Mn+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Fe+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
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
        print *, 'TestThermo75: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo75: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo75
