
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ThermoOutput.f90
    !> \brief   This subroutine determines which values are to be provided as output.
    !> \author  M.H.A. Piro
    !> \date    May 8, 2012
    !> \sa      Thermochimica.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose
    ! =======
    !
    !> \details The purpose of this subroutine is to store particular information that is requested as output.
    !! This is done to avoid unnecessary use of memory.  Note that it is possible that information is requested
    !! that is not in the system.  This subroutine first checks to see if the requested information is even
    !! in the system and then it checks to see if the particular species or phase is stable at equilibrium.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! GetSolnPhaseIndex     An integer scalar representing the absolute index of a solution phase.
    ! dSolnPhaseMolesOut    A double real vector representing the number of moles of particular solution phases
    !                        requested as output.  If this phase is not present at equilibrium, it is assigned
    !                        a value of zero.
    ! dPureConPhaseMolesOut A double real vector representing the number of moles of particular pure condensed 
    !                        phases requested as output.  If this phase is not present at equilibrium, it is 
    !                        assigned a value of zero.
    ! dSpeciesMoleFractionOut A double real vector representing the mole fraction of particular soluton species 
    !                        requested as output.  If this species is not present at equilibrium, it is 
    !                        assigned a value of zero.
    ! IsSolnPhaseInSys      A character string indicating whether a particular solution phase is in the system
    !                        ('yes') or not ('no').
    ! IsPureConPhaseInSys   A character string indicating whether a particular pure condensed phase is in the 
    !                        system ('yes') or not ('no').
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ThermoOutput
    
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none
    
    integer::   i, j, k, GetSolnPhaseIndex
    logical::   IsSolnPhaseInSys, IsPureConPhaseInSys
    
    
    ! Initialize variables:
    dSolnPhaseMolesOut      = 0D0
    dPureConPhaseMolesOut   = 0D0
    dSpeciesMoleFractionOut = 0D0

    ! Only proceed if there haven't been any previously encountered errors:
    if (INFOThermo == 0) then

        ! ----------------
        ! SOLUTION PHASES:
        ! ----------------

        LOOP_SolnPhasesOut: do i = 1, nSolnPhasesOut
            
            ! Check if this solution phase is in the system:
            if (IsSolnPhaseInSys(cSolnPhaseNameOut(i)) .EQV. .TRUE. ) then
                ! The solution phase requested is in the system. Check if it is stable at equilibrium:
                LOOP_nSolnPhases: do j = 1, nSolnPhases 
                    k = -iAssemblage(nElements - j + 1)         ! Absolute solution phase index in system.
                    if (cSolnPhaseNameOut(i) == cSolnPhaseName(k)) then
                        ! This solution phase is expected to be stable at equilibrium.  
                        ! Store the number of moles of this phase:
                        dSolnPhaseMolesOut(i) = dMolesPhase(nElements - j + 1)
                        exit LOOP_nSolnPhases
                    end if
                end do LOOP_nSolnPhases
            else
                ! The solution phase requested as output is not in the system.  Placeholder if an error should be 
                ! reported at a later time.  
            end if          
        end do LOOP_SolnPhasesOut
  
  
        ! -----------------
        ! SOLUTION SPECIES:
        ! -----------------
        
        LOOP_SpeciesOut: do i = 1, nSpeciesOut
            ! First, make sure that the solution phase that this species belongs to is in the system:
            if (IsSolnPhaseInSys(cSpeciesPhaseOut(i)) .EQV. .TRUE. ) then
                ! Get the solution phase index:
                k = GetSolnPhaseIndex(cSpeciesPhaseOut(i))
                
                if (k == 0) cycle LOOP_SpeciesOut
                
                cSpeciesNameOut(i) = ' ' // cSpeciesNameOut(i)
                
                ! Loop through all species in this solution phase:
                LOOP_SpeciesInSoln: do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)                    
                    if (cSpeciesNameOut(i) == cSpeciesName(j)) then
                        dSpeciesMoleFractionOut(i) = dMolFraction(j)
                        exit LOOP_SpeciesInSoln
                    end if
                end do LOOP_SpeciesInSoln
            else
                ! The solution phase that this species belongs to is not part of the system.
                ! This is a placeholder for later.
            end if
        end do LOOP_SpeciesOut
        
        
        ! ----------------------
        ! PURE CONDENSED PHASES:
        ! ----------------------
        
        LOOP_PureConPhasesOut: do i = 1, nPureConPhaseOut
            
            ! Check if this pure condensed phase is in the system:
            if (IsPureConPhaseInSys(cPureConPhaseNameOut(i)) .EQV. .TRUE.) then
                ! The pure condensed phase requested as output is in the system.
                ! Check if it is stable at equilibrium:
                LOOP_nConPhases: do j = 1, nConPhases
                    k = iAssemblage(j)      ! Absolute pure condensed phase index in system.
                    if (cPureConPhaseNameOut(i) == cSpeciesName(k)) then
                        ! This pure condensed phase is expected to be stable at equilibrium.
                        ! Store the number of moles of this phase:
                        dPureConPhaseMolesOut(i) = dMolesPhase(j)
                        exit LOOP_nConPhases
                    end if
                end do LOOP_nConPhases
            else 
                ! The pure condensed phase requested as output is not in the system.  Placeholder if an error 
                ! should be reported at a later time.  
            end if
        end do LOOP_PureConPhasesOut

    else
        ! Provide null values if an error has been detected.
    end if
    
    return
            
end subroutine ThermoOutput


    !-------------------------------------------------------------------------------------------------------------
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose
    ! =======
    !
    ! The purpose of this function is to determine whether a particular solution phase is in the system.  Note
    ! that this is not to be confused with a solution phase that is expected to be stable at equilibrium.
    !
    !-------------------------------------------------------------------------------------------------------------



function IsSolnPhaseInSys(cSolnPhaseNameOut)
    
    USE ModuleThermo
    
    implicit none

    integer::                 i
    character(*),intent(in):: cSolnPhaseNameOut
    logical::                 IsSolnPhaseInSys


    ! Initialize variables:
    IsSolnPhaseInSys = .FALSE.

    ! Check to see if this phase is even in the system:
    do i = 1, nSolnPhasesSys
        if (cSolnPhaseName(i) == cSolnPhaseNameOut) then
            IsSolnPhaseInSys = .TRUE.
            exit 
        end if
    end do   

end function IsSolnPhaseInSys



    !-------------------------------------------------------------------------------------------------------------
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose
    ! =======
    !
    ! The purpose of this function is to get the index of a requested solution phase.
    !
    !-------------------------------------------------------------------------------------------------------------



function GetSolnPhaseIndex(cSolnPhaseNameOut)
    
    USE ModuleThermo
    
    implicit none

    integer::                 i, GetSolnPhaseIndex
    character(*),intent(in):: cSolnPhaseNameOut


    ! Initialize variables:
    GetSolnPhaseIndex = 0

    ! Check to see if this phase is even in the system:
    do i = 1, nSolnPhasesSys
        if (cSolnPhaseName(i) == cSolnPhaseNameOut) then
            GetSolnPhaseIndex = i
            exit 
        end if
    end do   

end function GetSolnPhaseIndex


    !-------------------------------------------------------------------------------------------------------------
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose
    ! =======
    !
    ! The purpose of this function is to determine whether a particular pure condensed phase is in the system.  
    ! Note that this is not to be confused with a pure condensed phase that is expected to be stable at 
    ! equilibrium.
    !
    !-------------------------------------------------------------------------------------------------------------


function IsPureConPhaseInSys(cPureConPhaseNameOut)
    
    USE ModuleThermo
    
    implicit none
 
    integer::                 i 
    character(*),intent(in):: cPureConPhaseNameOut
    character(26)::           cTemp
    logical::                 IsPureConPhaseInSys


    ! Initialize variables:
    IsPureConPhaseInSys = .FALSE.
    
    cTemp = ' ' // cPureConPhaseNameOut

    ! Check to see if this phase is even in the system:
    do i = nSpeciesPhase(nSolnPhasesSys) + 1, nSpecies

        if (cSpeciesName(i) == cTemp) then
            IsPureConPhaseInSys = .TRUE.
            exit 
        end if
    end do   

end function IsPureConPhaseInSys
