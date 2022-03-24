
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo74.F90
    !> \brief   Spot test - 1900K with , 2 mol Fe, 4 mol Mn, 3 mol S.
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
    !!  results for the open literature Ca-Mn-S assessment file at 2500K with 1 mol of Ca, 1 mol of Mn, and
    !!  1 mol of S. It also tests mixing term Case #1, # 2, and #3 of the SUBI phase with a miscibility gap
    !!  present.
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: D. Dilner, "Thermodynamic description of the Fe–Mn–Ca–Mg–S system," Calphad,
    !!  vol. 53, pp. 55-61, 2016.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo74

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4
    real(8) :: sfcheck5, sfcheck6, sfcheck7
    real(8) :: pcheck1, pcheck2, pcheck3, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass
    logical :: s5pass, s6pass, s7pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'FeMnCaS_mod1.DAT'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1900D0
    dElementMass(26)       = 2D0          ! Fe
    dElementMass(25)       = 4D0          ! Mn
    dElementMass(16)       = 3D0          ! S

    ! Liquid #1
    sfcheck1 = 3.2331D-1         !Fe+2
    sfcheck2 = 6.7669D-1         !Mn+2
    sfcheck3 = 5.1140D-1         !S-2
    sfcheck4 = 4.8860D-1         !Va

    ! Liquid #2
    sfcheck5 = 7.7271D-1         !Fe+2
    sfcheck6 = 2.27D-1           !Mn+2
    sfcheck7 = 9.9999D-1         !Va

    pcheck1 = -121307D0       !Ca
    pcheck2 = -155496D0       !Mn
    pcheck3 = -360095D0       !S

    gibbscheck = -1944880D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    lFuzzyStoich = .FALSE.
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

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (gibbscheck))/(gibbscheck) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3).AND. &
            (DABS((dElementPotential(3)*dIdealConstant*dTemperature - pcheck3)/pcheck3) < 1D-3)) then
            do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(k) == 'IONIC_LIQ#2') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Fe+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Mn+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'S-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck4)/sfcheck4 < 1D-3) s4pass = .TRUE.
                            end if
                        end do
                    end do
                else if (cSolnPhaseName(k) == 'IONIC_LIQ#1') then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Fe+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck5)/sfcheck5 < 1D-3) s5pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Mn+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck6)/sfcheck6 < 1D-3) s6pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck7)/sfcheck7 < 1D-3) s7pass = .TRUE.
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
        s4pass .AND. &
        s5pass .AND. &
        s6pass .AND. &
        s7pass) then
        ! The test passed:
        print *, 'TestThermo74: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo74: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo74
