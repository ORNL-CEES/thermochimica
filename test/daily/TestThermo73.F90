
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo73.F90
    !> \brief   Spot test - 2500K with 33% Ca, 33% Mn, 33% S.
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
    !    17/04/2024    A.E.F. Fitzsimmons  Naming convention change
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

program TestThermo73

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3, sfcheck4, sfcheck5
    real(8) :: sfcheck6, sfcheck7, sfcheck8, sfcheck9, sfcheck10
    real(8) :: pcheck1, pcheck2, pcheck3, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass, s4pass, s5pass
    logical :: s6pass, s7pass, s8pass, s9pass, s10pass


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'CaMnS.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(20)        = 1D0          ! Ca
    dElementMass(25)        = 1D0          ! Mn
    dElementMass(16)        = 1D0          ! S

    ! Liquid #1
    sfcheck1 = 8.06572D-1         !Ca
    sfcheck2 = 1.93428D-1         !Mn
    sfcheck3 = 8.18084D-1         !S-2
    sfcheck4 = 1.81915D-1         !Va
    sfcheck5 = 7.69522D-7         !S

    ! Liquid #2
    sfcheck6 = 1.96945D-2         !Ca
    sfcheck7 = 9.80306D-1         !Mn
    sfcheck8 = 1.6531D-3          !S-2
    sfcheck9 = 9.98344D-1         !Va
    sfcheck10 = 3.04003D-6        !S

    pcheck1 = -251878D0       !Ca
    pcheck2 = -201710D0       !Mn
    pcheck3 = -457030D0       !S

    gibbscheck = -910619D0

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
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Ca+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Mn+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'S-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck4)/sfcheck4 < 1D-3) s4pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'S') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck5)/sfcheck5 < 1D-3) s5pass = .TRUE.
                            end if
                        end do
                    end do
                
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Ca+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck6)/sfcheck6 < 1D-3) s6pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Mn+2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck7)/sfcheck7 < 1D-3) s7pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'S-2') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck8)/sfcheck8 < 1D-3) s8pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck9)/sfcheck9 < 1D-3) s9pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'S') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck10)/sfcheck10 < 1D-3) s10pass = .TRUE.
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
        s8pass .AND. &
        s9pass .AND. &
        s10pass) then
        ! The test passed:
        print *, 'TestThermo73: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo73: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo73
