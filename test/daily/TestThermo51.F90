
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo76.F90
    !> \brief   Spot test - 2500K with 1.0 Cs - 0.5 Te
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
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the open literature Cs-Te assessment file at 2500K with 1.0 mol of Cs and 0.5 mol of Te.
    !!  It also tests mixing term Case #9 of the SUBI phase.
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: T. N. Pham Thi, J. C. Dumas, V. Bouineau, N. Dupin, C. Gueneau, S. Gosse, 
    !!  P. Benigni, P. Maugis and J. Rogez, "Thermodynamic assessment of the Cs–Te binary system," Calphad,
    !!  vol. 48, pp. 1-12, 2015.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo51

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
    cThermoFileName        = DATA_DIRECTORY // 'CsTe-2.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(52)       = 0.5D0          ! Te
    dElementMass(55)       = 1.0D0          ! Cs

    ! Init test values
    dGibbsCheck            = -5.68905D05
    dHeatCapacityCheck     = 84.2722
    nSpeciesTest           = 3
    iSpeciesIndexTest      = [1, 2, 3] !Cs+:VA, Cs2Te, Te
    dMolFractionTest       = [3.4222D-01, 4.8666D-01, 1.7111D-01]
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
        print *, 'TestThermo51: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo51: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo51
