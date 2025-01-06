
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
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
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

program TestThermo52

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    ! Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq_mod1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(41)       = 1D0          ! Nb
    dElementMass(8)        = 0.3D0        ! O
    dElementMass(50)       = 0.7D0        ! Sn

    ! Init test values
    dGibbsCheck            = -5.14706D05
    dHeatCapacityCheck     = -1.7023D09
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [6, 14, 20, 21] ! O2Sn, NB+2:O-2, O2Sn, NB+2:O-2
    dMolFractionTest       = [4.5414D-11, 5.1351D-01, 1.0010D-027,  2.1018D-03]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    ! Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo52: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo52: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo52
