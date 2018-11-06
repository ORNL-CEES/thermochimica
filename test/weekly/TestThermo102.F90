

    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date           Programmer          Description of change
    !    ----           ----------          ---------------------
    !    07/17/2013     M.H.A. Piro         Original code
    !    05/13/2014     M.H.A. Piro         Write data to a data-file and check to see if the epsilon phase
    !                                        is predicted to be stable (which would be incorrect).
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! sample problem representing the Zr-H binary system.  Unlike other tests involving CEF phases, this particular
    ! system contains phases where there is only one constituent on the first sublattice and the mixing terms
    ! represent constituents mixing on the second sublattice (i.e., other phases mix constituents on the first
    ! sublattice).  A reference for the thermodynamic database for this test follows:
    !
    !   N. Dupin, I. Ansara, C. Servant, C. Toffolon, C. Lemaignan, J.C. Brachet, "A Thermodynamic Database for
    !   Zirconium Alloys," Journal of Nuclear Materials, 275 (1999) 287-295.
    !
    ! This test carpet bombs a region of the Zr-H phase diagram to ensure that the epsilon phase does not form.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo102

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer:: i, j, b, iCount
    real(8):: dTemp


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'UPMZ_MHP.dat'

    ! Specify values:
    dPressure              = 1D0
    iCount                 = 0

    ! Initialize random number generator:
    !call init_random_seed()

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Temperature loop:
    LOOP_T: do i = 300, 1500, 1

        ! Definte temperature:
        dTemperature = DFLOAT(i)

        ! Composition loop:
        LOOP_B: do b = 1, 100, 1

            ! Define elemental masses:
            do j = 1, 118

                call random_number(dTemp)
                dElementMass(j) = dTemp

            end do

            iCount = iCount + 1
            print *, 'iCount = ', iCount

            print *, dTemperature, dElementMass(94), dElementMass(92), dElementMass(42), dElementMass(40)

            ! Call Thermochimica:
            call Thermochimica

            if (INFOThermo /= 0) exit LOOP_T

        end do LOOP_B   ! End composition loop
    end do LOOP_T       ! End temperature loop

    if ((INFOThermo == 0)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo102: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call PrintResults
        call EXIT(0)
    else
        ! The unit test failed.
        print *, 'TestThermo102: FAIL <--- '
        ! Reset Thermochimica:
        call ResetThermo
        call PrintResults
        call EXIT(0)
    end if

end program TestThermo102
