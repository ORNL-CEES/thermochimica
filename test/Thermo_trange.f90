

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
    ! The purpose of this application test is to ensure that the subroutine GetOutputSiteFraction.f90 works
    ! correctly. The purpose of that subroutine is to get the site fraction of a particular constituent on
    ! a particular sublattice of a particular phase, which could be used by some other multi-physics
    ! code.  For example, simulating solid state diffusion in some multi-physics code may need the site
    ! fraction of interstitials of some component in a multi-subllatice phase.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo59

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    integer       :: INFO
    integer       :: iSublatticeOut, iConstituentOut
    character(25) :: cSolnOut
    real(8)       :: dSiteFractionOut, A, B

    integer       :: i, j, k, l
    character(3)  :: cElName
    real(8)       :: dDataOut, T1, T2, DeltT, dOtoMRatio

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = '../data/Pu_U_O_CEA.dat'

    T1= 500D0
    T2=2000D0
    DeltT=(T2-T1)/15.0D0
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    do i=1,16
       dTemperature = T1 + (i-1)*DeltT
    
    ! Initialize variables:
    ! dTemperature            = 1234D0
    dPressure               = 1D0    
    dElementMass            = 0D0
    dElementMass(8)         = 2.05D0
    ! dElementMass(8)         = 1.95D0
    dElementMass(92)        = 1D0
    ! dElementMass(8)         = 81260.0D0 * 2.05D0/2.0D0
    ! dElementMass(92)        = 40630.0D0

    ! Specify output variables to be fetched:
    cSolnOut          = 'O2ZRU_C'
    dSiteFractionOut  = 0D0

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Get the site fraction of vacanies on the normal anion sublattice:
    iSublatticeOut  = 2
    iConstituentOut = 2
    call GetOutputSiteFraction(cSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)
    A = dSiteFractionOut

    ! Get the site fraction of oxygen interstitials:
    iSublatticeOut  = 3
    iConstituentOut = 1
    call GetOutputSiteFraction(cSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)
    B = dSiteFractionOut

    !-------------------------------------------------------------------------------------------------------------
    !
    ! Commented out the following lines for the purposes of running application tests.
    ! These lines should indicate how one could use these variables.
    !
    !-------------------------------------------------------------------------------------------------------------

    !print *, ' The site fraction of vacanies on the normal anion sublattice is ', A, ', and the site fraction'
    !print *, ' of oxygen interstitials is ', B, '.'
    !print *


    iPrintResultsMode = 2
    ! Perform post-processing of results:
    ! if (iPrintResultsMode > 0) call PrintResults

    call CompOtoMRatio(cSolnOut, dOtoMRatio, INFO)
    
    cElName='O'
    call GetOutputChemPot(cElName, dDataOut, INFO)
    write(*,'(A3,f6.1,A4,E21.14,A3,A1,A3,E12.5,A4,E12.5,A4,E21.14)') &
         'T= ', dTemperature, ' O/M= ', dOtoMRatio, ' G(', TRIM(cElName), &
         ')= ', dDataOut, ' v= ', A, ' i= ', B

    do j = 1, nSolnPhases
       ! Relative and absolute solution phase indices, respectively:
       k      = nElements - j + 1
       l      = -iAssemblage(k)
       print '(A15,E12.4)', cSolnPhaseName(l), dMolesPhase(k)
    enddo

    do k = 1, nConPhases
       j = iAssemblage(k)
        print '(A15,E12.4)', cSpeciesName(iAssemblage(k)), dMolesPhase(k)
    enddo
    
    ! Reset Thermochimica:
    call ResetThermo

    enddo
 
    call EXIT(STATUS)
end program TestThermo59
