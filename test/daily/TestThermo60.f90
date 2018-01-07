

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
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    06/10/2016    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to give bogus input to the subroutine GetOutSiteFraction, which
    ! should return an error and exit gracefully.  Specifically, a solution phase is requested to get site
    ! fraction info that is included in the data-file, but is not stable under the conditions specified.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo60

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    integer       :: INFO
    integer       :: iSublatticeOut, iConstituentOut
    character(25) :: cSolnOut
    real(8)       :: dSiteFractionOut


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = '../data/Pu_U_O_CEA.dat'

    ! Initialize variables:
    dTemperature            = 400.0D0
    dPressure               = 1.0D0
    dElementMass            = 0.0D0
    dElementMass(8)         = 7.5D0
    dElementMass(92)        = 3.0D0

    ! Specify output variables to be fetched:
    cSolnOut          = 'O2ZRU_C'
    dSiteFractionOut  = 0D0
    iPrintResultsMode = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Get the site fraction of vacanies on the normal anion sublattice:
    iSublatticeOut  = 2
    iConstituentOut = 2
    call GetOutputSiteFraction(cSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)

    ! This particular test is expecting GetOutputSiteFraction to return an error.
    if (INFO == 1) then
        print *, 'TestThermo60: PASS'
    else
        ! The test failed.
        print *, 'TestThermo60: FAIL <---'
        STATUS  = 1
    end if

    ! iPrintResultsMode = 2
    ! Perform post-processing of results:
    ! if (iPrintResultsMode > 0) call PrintResults

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo60
