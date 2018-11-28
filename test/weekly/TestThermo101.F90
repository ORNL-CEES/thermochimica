
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo101.F90
    !> \brief   Convergence test - Zr-H.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date           Programmer          Description of change
    !    ----           ----------          ---------------------
    !    07/17/2013     M.H.A. Piro         Original code
    !    05/13/2014     M.H.A. Piro         Write data to a data-file and check to see if the epsilon phase
    !                                        is predicted to be stable (which would be incorrect).
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct results
    !! for a sample problem representing the Zr-H binary system.  Unlike other tests involving CEF phases, this
    !! particular system contains phases where there is only one constituent on the first sublattice and the
    !! mixing terms represent constituents mixing on the second sublattice (i.e., other phases mix constituents
    !! on the first sublattice).  A reference for the thermodynamic database for this test follows:
    !!
    !!   N. Dupin, I. Ansara, C. Servant, C. Toffolon, C. Lemaignan, J.C. Brachet, "A Thermodynamic Database for
    !!   Zirconium Alloys," Journal of Nuclear Materials, 275 (1999) 287-295.
    !!
    !! This test carpet bombs a region of the Zr-H phase diagram to ensure that the epsilon phase does not form.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo101

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer:: i, k, m, b, INFO

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZRHD_MHP.dat'

    ! Specify values:
    dPressure              = 1D0
    iCounter               = 0
    INFO                   = 0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Create a data-file to output results:
    open (unit = 4, file = 'output_TestThermo101.dat')
    write (4,*) 'iCounter | Temperature [K] | Mole Fraction H'

    ! Temperature loop:
    LOOP_T: do i = 300, 900, 1

        ! Composition loop:
        LOOP_B: do b = 1, 100, 1

            ! Define elemental masses:
            dElementMass(1)  = DFLOAT(b) * 1D-3
            dElementMass(40) = 1D0 - dElementMass(1)

            ! Definte temperature:
            dTemperature = DFLOAT(i)

            ! Call Thermochimica:
            call Thermochimica

            ! Exit if an error has been found:
            if (INFOThermo /= 0) exit LOOP_T

            ! Check results:
            if (nSolnPhases == 1) then
                k = nElements
                k = -iAssemblage(k)
                write (4,*) INFOThermo, dTemperature, dElementMass(1) / SUM(dElementMass), cSolnPhaseName(k)

                ! Check to see if the epsilon phase is stable:
                if (cSolnPhaseName(k) == 'ZRH2_EPSILON') then
                    INFO = 101
                    exit LOOP_T
                end if

            elseif (nSolnPhases == 2) then
                k = nElements
                k = -iAssemblage(k)
                m = nElements - 1
                m = -iAssemblage(m)
                write (4,*) INFOThermo, dTemperature, dElementMass(1) / SUM(dElementMass), cSolnPhaseName(k), cSolnPhaseName(m)

                ! Check to see if the epsilon phase is stable:
                if ((cSolnPhaseName(k) == 'ZRH2_EPSILON').OR.(cSolnPhaseName(m) == 'ZRH2_EPSILON')) then
                    INFO = 101
                    exit LOOP_T
                end if

            else
                ! An unexpected error occured.  Record an error code and exit.
                INFO = 10
                exit LOOP_T
            end if

        end do LOOP_B   ! End composition loop
    end do LOOP_T       ! End temperature loop

    if ((INFOThermo == 0).AND.(INFO == 0)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo101: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call PrintResults
        call EXIT(0)
    else
        ! The unit test failed.
        print *, 'TestThermo101: FAIL <--- '
        ! Reset Thermochimica:
        call ResetThermo
        call PrintResults
        call EXIT(1)
    end if

end program TestThermo101
