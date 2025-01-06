
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

program TestThermo46

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
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(8)        = 0.1D0          ! O
    dElementMass(24)       = 0.2D0          ! Cr
    dElementMass(40)       = 0.7D0          ! Zr

    ! Init test values
    dGibbsCheck            = -1.94030D05
    dHeatCapacityCheck     = 43.8857
    nSpeciesTest           = 5
    iSpeciesIndexTest      = [1, 10, 18, 44, 50] !Cr, Cr+3:O-2, Cr:O, ZR:CR, ZR:O
    dMolFractionTest       = [3.7689D-04, 2.4239D-07, 6.8007D-06, 7.5995D-05, 1.0000D-05]
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
        print *, 'TestThermo46: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo46: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo46
