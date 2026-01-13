
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo31.F90
    !> \brief   Spot test - W-Au-Ar-O_02.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer           Description of change
    !    ----          ----------           ---------------------
    !    05/14/2013    M.H.A. Piro          Original code
    !    08/31/2018    B.W.N. Fitzpatrick   Change to a fictive database
    !    04/17/2024    A.E.F. Fitzsimmons   Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons   SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive W-Au-Ar-O_02 system.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo16

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
    cThermoFileName        = DATA_DIRECTORY // 'WAuArO-2.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000D0
    dElementMass(74)       = 1D0        ! W
    dElementMass(79)       = 3D0        ! Au
    dElementMass(18)       = 5D0        ! Ar
    dElementMass(8)        = 2D0        ! O

    ! Init test values
    dGibbsCheck            = 6.76880D05
    dHeatCapacityCheck     = 251.212
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [1, 2, 3, 4] !Au, W, O, Ar
    dMolFractionTest       = [0.272727D0, 9.09090D-02, 0.18181D0, 0.45454D0 ]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity
    
    ! Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    ! Check results:
    if (lPass) then
        ! The test passed:
        print *, 'TestThermo16: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo16: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if
    
    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo16
