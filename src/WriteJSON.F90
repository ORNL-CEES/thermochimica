subroutine WriteJSON(append)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    logical, intent(in) :: append
    logical :: exist
    integer :: i

    ! Only proceed for a successful calculation:
    IF_PASS: if (INFOThermo /= 0) return

    inquire(file= DATA_DIRECTORY // '../thermoout.json', exist=exist)
    if (append .AND. exist) then
        open(1, file= DATA_DIRECTORY // '../thermoout.json', status='OLD', &
            position='append', action='write')
    else
        open(1, file= DATA_DIRECTORY // '../thermoout.json', status='REPLACE', &
              action='write')
    end if

    write(1,*) '{'

    ! Print the results for solution phases:
    call WriteJSONSolnPhase

    ! Print the results for pure condensed phases:
    call WriteJSONPureConPhase

    write(1,*) '  "elements": {'
    do i = 1, nElements
        write(1,*) '    "', TRIM(cElementName(i)), '": {'
        write(1,*) '      "moles": ', dMolesElement(i), ','
        write(1,*) '      "element potential": ', dElementPotential(i) * dIdealConstant * dTemperature
        if (i < nElements) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if
    end do
    write(1,*) '  },'

    write(1,*) '  "temperature": ', dTemperature, ','
    write(1,*) '  "pressure": ', dPressure, ','
    write(1,*) '  "Integral Gibbs energy": ', dGibbsEnergySys, ','
    write(1,*) '  "Functional norm": ', dGEMFunctionNorm, ','
    write(1,*) '  "solution phases": ', nSolnPhases, ','
    write(1,*) '  "pure condensed phases": ', nConPhases

    write(1,*) '}'

    close (1)

    return

end subroutine WriteJSON

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

subroutine WriteJSONSolnPhase

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    integer :: c, i, j, k, l, s, iFirst, iLast, iChargedPhaseID
    real(8) :: Tcritical, B, StructureFactor

    write(1,*) '  "solution phases": {'

    ! Loop through solution phases:
    LOOP_SolnStable: do j = 1, nSolnPhasesSys

        ! First and last solution species indices, respectively:
        iFirst = nSpeciesPhase(j-1) + 1
        iLast  = nSpeciesPhase(j)

        ! Print solution phase name:
        write(1,*) '    "', TRIM(ADJUSTL(cSolnPhaseName(j))), '": {'
        l = 0
        do k = 1, nElements
            if (-iAssemblage(k) == j) l = k
        end do
        if (l > 0) then
            write(1,*) '      "moles": ', dMolesPhase(l), ','
        else
            write(1,*) '      "moles": 0.0,'
        end if

        if ((cSolnPhaseType(j) == 'SUBLM') .OR. (cSolnPhaseType(j) == 'RKMPM')) then
            Tcritical = 0D0
            B = 0D0
            call CompMagneticTemperatureMoment(j,Tcritical,B)
            StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
            if (Tcritical < 0D0) then
                Tcritical = -Tcritical * StructureFactor
                write(1,*) '      "Neel temperature": ', Tcritical, ','
            else
                write(1,*) '      "Curie temperature": ', Tcritical, ','
            end if
            if (B < 0D0) then
                B         = -B * StructureFactor
            end if
            write(1,*) '      "Magnetic moment per atom": ', B, ','
        end if

        ! if ((cSolnPhaseType(l)) == 'SUBG' .OR. (cSolnPhaseType(l) == 'SUBQ')) then
        !     call CalculateCompositionSUBG(iSolnIndex=l,dMolesPairs=dMolesPairs,lPrint=.TRUE.)
        ! end if

        write(1,*) '      "species": {'

        ! Print species:
        do i = iFirst, iLast
            write(1,*) '        "', TRIM(ADJUSTL(cSpeciesName(i))), '": {'
            write(1,*) '          "mole fraction":', dMolFraction(i), ","
            write(1,*) '          "chemical potential":', dChemicalPotential(i)*dIdealConstant*dTemperature, ","
            write(1,*) '          "stoichiometry": [', (dStoichSpecies(i,c), ',', c = 1,nElements-1), &
                                                        dStoichSpecies(i,nElements), ']'
            if (i < iLast) then
                write(1,*) '        },'
            else
                write(1,*) '        }'
            end if
        end do
        write(1,*) '      }'

        if (j < nSolnPhasesSys) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if

        ! Check if this phase is represented by the Compound Energy Formalism:
        IF_SUBL: if ((csolnPhaseType(j) == 'SUBL').OR.(csolnPhaseType(j) == 'SUBLM')) then
            ! Store the index # of the charged phase:
            iChargedPhaseID = iPhaseSublattice(j)

            ! Loop through sublattices:
            LOOP_Sub: do s = 1, nSublatticePhase(iChargedPhaseID)

                write(1,*) 'Sublattice ', s, '; stoichiometric coefficient: ', &
                    dStoichSublattice(iChargedPhaseID,s)

                ! Loop through constituents in sublattice s:
                LOOP_SubCon: do c = 1, nConstituentSublattice(iChargedPhaseID,s)

                    write(1,*) cConstituentNameSUB(iChargedPhaseID,s,c), &
                            dSiteFraction(iChargedPhaseID,s,c)
                end do LOOP_SubCon

                if (s /= nSublatticePhase(iChargedPhaseID)) write(1,*)
            end do LOOP_Sub

        end if IF_SUBL

    end do LOOP_SolnStable

    write(1,*) '  },'

end subroutine WriteJSONSolnPhase

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

subroutine WriteJSONPureConPhase

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    integer                                 :: c, i, k, l

    write(1,*) '  "pure condensed phases": {'

    ! Loop through stable pure condensed phases:
    LOOP_PureStable: do i = nSpeciesPhase(nSolnPhasesSys) + 1, nSpeciesPhase(nSolnPhasesSys) + nConPhasesSys
        ! j = iTempVec(i)
        write(1,*) '    "', TRIM(ADJUSTL(cSpeciesName(i))), '": {'
        l = 0
        do k = 1, nElements
            if (iAssemblage(k) == i) l = k
        end do
        if (l > 0) then
            write(1,*) '      "moles":', dMolesPhase(l), ","
        else
            write(1,*) '      "moles": 0.0', ","
        end if
        write(1,*) '      "chemical potential":', dStdGibbsEnergy(i)*dIdealConstant*dTemperature, ","
        write(1,*) '      "stoichiometry": [', (dStoichSpecies(i,c), ',', c = 1,nElements-1), &
                                                    dStoichSpecies(i,nElements), ']'
        if (i < nSpeciesPhase(nSolnPhasesSys) + nConPhasesSys) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if
    end do LOOP_PureStable

    write(1,*) '  },'

end subroutine WriteJSONPureConPhase
