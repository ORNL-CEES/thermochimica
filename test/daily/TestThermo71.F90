
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo71.F90
    !> \brief   Spot test - 2000K with 20% Cr, 70% Zr, 10% O.
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
    !\details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the Zirc Data file at 1500K with 0.7 mols of Cn, 0.2 mols of Zr and 0.1 of O. It also
    !!  tests mixing term Case #2 an #3 of the SUBI phase. Permission was granted from N. Dupin to make
    !!  use of the Zirc DAT file.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo71

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: sfcheck1, sfcheck2, sfcheck3

    real(8) :: pcheck1, pcheck2, pcheck3, gibbscheck
    integer :: i,j,k,l
    logical :: s1pass, s2pass, s3pass



    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(8)         = 0.1D0          ! O
    dElementMass(24)        = 0.2D0          ! Cr
    dElementMass(40)        = 0.7D0          ! Zr

    ! Liquid #1
    sfcheck1 = 3.30097D-1      !Cr
    sfcheck2 = 6.69903D-1      !Zr
    sfcheck3 = 9.99997D-1      !Va

    pcheck1 = -126330D0        ! Cr
    pcheck2 = -628950D0        ! O
    pcheck3 = -151241D0        ! Zr

    gibbscheck = -194030D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    s1pass = .FALSE.
    s2pass = .FALSE.
    s3pass = .FALSE.
    !s4pass = .FALSE.

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (gibbscheck))/(gibbscheck) < 1D-3) .AND. &
            (DABS((dElementPotential(1)*dIdealConstant*dTemperature - pcheck1)/pcheck1) < 1D-3).AND. &
            (DABS((dElementPotential(2)*dIdealConstant*dTemperature - pcheck2)/pcheck2) < 1D-3).AND. &
            (DABS((dElementPotential(3)*dIdealConstant*dTemperature - pcheck3)/pcheck3) < 1D-3)) then
            loop_checkPhases: do i = 1, nSolnPhases
                k = -iAssemblage(nElements + 1 - i)
                if ((cSolnPhaseName(k) == 'IONIC_LIQ#1') .OR. (cSolnPhaseName(k) == 'IONIC_LIQ#2')) then
                    do j = 1, 2
                        do l = 1, nConstituentSublattice(i,j)
                            if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'CR+3') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck1)/sfcheck1 < 1D-3) s1pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'ZR+4') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck2)/sfcheck2 < 1D-3) s2pass = .TRUE.
                            else if (TRIM(ADJUSTL(cConstituentNameSUB(i,j,l))) == 'Va') then
                                if (DABS(dSiteFraction(i,j,l) - sfcheck3)/sfcheck3 < 1D-3) s3pass = .TRUE.
                            end if
                        end do
                    end do
                    exit loop_checkPhases
                end if
            end do loop_checkPhases
        end if
    end if

    if (s1pass .AND. &
        s2pass .AND. &
        s3pass) then
        ! The test passed:
        print *, 'TestThermo71: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo71: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo71
