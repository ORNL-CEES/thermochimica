
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo102.F90
    !> \brief   Convergence test - currently broken test with no input.
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
    !> \details 
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
    cThermoFileName        = DATA_DIRECTORY // 'UPMZ_MHP.dat'

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
        call EXIT(1)
    end if

end program TestThermo102
