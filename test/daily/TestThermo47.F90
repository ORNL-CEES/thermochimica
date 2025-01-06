
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo72.F90
    !> \brief   Spot test - 700K with 20% Cs, 80% Te.
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
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !  results for the open literature Cs-Te assessment file at 700K with 0.2 mols of Cs and 0.8 of Te.
    !  It also tests mixing term Case #5 of the SUBI phase.
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: T. N. Pham Thi, J. C. Dumas, V. Bouineau, N. Dupin, C. Gueneau, S. Gosse,
    !!  P. Benigni, P. Maugis and J. Rogez, "Thermodynamic assessment of the Cs–Te binary system," Calphad,
    !!  vol. 48, pp. 1-12, 2015.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo47

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    !Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'CsTe-1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 700D0
    dElementMass(52)       = 1D0          ! Te
    dElementMass(55)       = 1D0          ! Cs

    !Init test values
    dGibbsCheck            = -2.93180D05
    dHeatCapacityCheck     = 9.26111D01
    nSpeciesTest           = 2
    iSpeciesIndexTest      = [16, 17] !Cs2Te, Te
    dMolFractionTest       = [0.41360D0, 0.58640D0]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo47: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo47: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if


end program TestThermo47
