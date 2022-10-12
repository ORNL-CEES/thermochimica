subroutine WriteJSON(append)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    logical, intent(in) :: append
    logical :: exist
    integer :: i, c, nElectron, its

    inquire(file= DATA_DIRECTORY // '../outputs/thermoout.json', exist=exist)
    if (append .AND. exist) then
        open(1, file= DATA_DIRECTORY // '../outputs/thermoout.json', &
              status='OLD', position='append', action='write')
    else
        open(1, file= DATA_DIRECTORY // '../outputs/thermoout.json', &
              status='REPLACE', action='write')
    end if

    write(1,*) '{'

    ! Only proceed for a successful calculation:
    if (INFOThermo /= 0) then
        write(1,*) '}'
        close (1)
        return
    end if

    ! Print the results for solution phases:
    call WriteJSONSolnPhase

    ! Print the results for pure condensed phases:
    call WriteJSONPureConPhase

    nElectron = 0
    do c = 1, nElements
        if (cElementName(c) == 'e-') nElectron = nElectron + 1
    end do

    write(1,*) '  "elements": {'
    do i = 1, nElements - nElectron
        write(1,*) '    "', TRIM(cElementName(i)), '": {'
        write(1,*) '      "moles": ', dMolesElement(i), ','
        write(1,'(A28,ES25.16E3)') '      "element potential": ', dElementPotential(i) * dIdealConstant * dTemperature
        if (i < nElements - nElectron) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if
    end do
    write(1,*) '  },'

    write(1,*) '  "temperature": ', dTemperature, ','
    write(1,*) '  "pressure": ', dPressure, ','
    write(1,'(A28,ES25.16E3,A1)') '  "integral Gibbs energy": ', dGibbsEnergySys, ','
    write(1,'(A14,ES25.16E3,A1)') '  "entropy": ', dEntropy, ','
    write(1,'(A15,ES25.16E3,A1)') '  "enthalpy": ', dEnthalpy, ','
    write(1,'(A20,ES25.16E3,A1)') '  "heat capacity": ', dHeatCapacity, ','
    write(1,*) '  "functional norm": ', dGEMFunctionNorm, ','
    its = iterGlobal
    if (lRetryAttempted) its = its + iterGlobalMax
    write(1,*) '  "GEM iterations": ', iterGlobal, ','
    write(1,*) '  "# solution phases": ', nSolnPhases, ','
    write(1,*) '  "# pure condensed phases": ', nConPhases

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

    implicit none

    integer :: c, i, j, k, l, s, iFirst, iLast, iChargedPhaseID, nElectron
    real(8) :: Tcritical, B, StructureFactor, dTempMolesPhase, dTotalElements, dCurrentElement
    character(16) :: intStr
    character(25) :: tempSpeciesName

    write(1,*) '  "solution phases": {'

    ! Loop through solution phases:
    do j = 1, nSolnPhasesSys

        ! First and last solution species indices, respectively:
        iFirst = nSpeciesPhase(j-1) + 1
        iLast  = nSpeciesPhase(j)

        ! Print solution phase name:
        if (lMiscibility(j)) then
            write(1,'(A,A,A,I0,A)') '     "', TRIM(ADJUSTL(cSolnPhaseName(j))), '#', j, '": {'
        else
            write(1,*) '    "', TRIM(ADJUSTL(cSolnPhaseName(j))), '": {'
        end if
        write(1,*) '      "phase model": "', TRIM(ADJUSTL(cSolnPhaseType(j))), '",'
        l = 0
        do k = 1, nElements
            if (-iAssemblage(k) == j) l = k
        end do
        if (l > 0) then
            write(1,'(A16,ES25.16E3,A1)') '      "moles": ', dMolesPhase(l), ','
            dTempMolesPhase = dMolesPhase(l)
        else
            write(1,*) '      "moles": 0.0,'
            dTempMolesPhase = 0D0
        end if

        ! Print driving force
        write(1,'(A24,ES25.16E3,A1)') '      "driving force": ', dDrivingForceSoln(j), ','

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

        if ((cSolnPhaseType(j)) == 'SUBG' .OR. (cSolnPhaseType(j) == 'SUBQ')) then
            call WriteJSONMQM(j)
        else
            write(1,*) '      "species": {'
            ! Print species:
            do i = iFirst, iLast
                write(tempSpeciesName,100) cSpeciesName(i)
                100 FORMAT (A25)
                write(1,*) '        "', TRIM(ADJUSTL(tempSpeciesName)), '": {'
                write(1,*) '          "mole fraction":', dMolFraction(i), ','
                write(1,'(A19,ES25.16E3,A1)') '          "moles":', dMolFraction(i)*dTempMolesPhase, ','
                write(1,*) '          "chemical potential":', dChemicalPotential(i)*dIdealConstant*dTemperature, ','
                write(1,*) '          "stoichiometry": [', (dStoichSpecies(i,c), ',', c = 1,nElements-1), &
                                                            dStoichSpecies(i,nElements), ']'
                if (i < iLast) then
                    write(1,*) '        },'
                else
                    write(1,*) '        }'
                end if
            end do

            ! Check if this phase is represented by the Compound Energy Formalism:
            if ((csolnPhaseType(j) == 'SUBL').OR.(csolnPhaseType(j) == 'SUBLM')) then
                write(1,*) '      },'
                write(1,*) '      "sublattices": {'
                ! Store the index # of the charged phase:
                iChargedPhaseID = iPhaseSublattice(j)

                ! Loop through sublattices:
                do s = 1, nSublatticePhase(iChargedPhaseID)
                    write(intStr,*) s
                    write(1,*) '        "sublattice ', TRIM(ADJUSTL(intStr)), '": {'
                    write(1,*) '          "stoichiometric coefficient": ', dStoichSublattice(iChargedPhaseID,s), ','
                    write(1,*) '          "constituents": {'
                    ! Loop through constituents in sublattice s:
                    LOOP_SubCon: do c = 1, nConstituentSublattice(iChargedPhaseID,s)
                        if (c < nConstituentSublattice(iChargedPhaseID,s)) then
                            write(1,*) '            "', TRIM(ADJUSTL(cConstituentNameSUB(iChargedPhaseID,s,c))), '": ', &
                                    dSiteFraction(iChargedPhaseID,s,c), ','
                        else
                            write(1,*) '            "', TRIM(ADJUSTL(cConstituentNameSUB(iChargedPhaseID,s,c))), '": ', &
                                    dSiteFraction(iChargedPhaseID,s,c)
                        end if
                    end do LOOP_SubCon
                    write(1,*) '          }'
                    if (s < nSublatticePhase(iChargedPhaseID)) then
                        write(1,*) '        },'
                    else
                        write(1,*) '        }'
                    end if
                end do
            end if
            write(1,*) '      },'
        end if

        write(1,*) '      "elements": {'
        dTotalElements = 0D0
        nElectron = 0
        do c = 1, nElements
            do i = iFirst, iLast
                dTotalElements = dTotalElements + dStoichSpecies(i,c)*dMolFraction(i)
            end do
            if (cElementName(c) == 'e-') nElectron = nElectron + 1
        end do
        do c = 1, nElements - nElectron
            dCurrentElement = 0D0
            do i = iFirst, iLast
                dCurrentElement = dCurrentElement + dStoichSpecies(i,c)*dMolFraction(i)
            end do
            write(1,*) '        "', TRIM(cElementName(c)), '": {'
            write(1,'(A39,ES25.16E3,A1)') '          "moles of element in phase":', dCurrentElement*dTempMolesPhase, ','
            write(1,*) '          "mole fraction of phase by element":', dCurrentElement / dTotalElements, ','
            write(1,*) '          "mole fraction of element by phase":', dCurrentElement*dTempMolesPhase / dMolesElement(c)
            if (c < nElements - nElectron) then
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

    end do

    write(1,*) '  },'

end subroutine WriteJSONSolnPhase

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

subroutine WriteJSONPureConPhase

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: c, i, j, k, l, nElectron
    real(8) :: dTempMolesPhase, dTotalElements, dCurrentElement, dDriving

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
            write(1,'(A16,ES25.16E3,A1)') '      "moles": ', dMolesPhase(l), ','
            dTempMolesPhase = dMolesPhase(l)
        else
            write(1,*) '      "moles": 0.0,'
            dTempMolesPhase = 0D0
        end if
        write(1,*) '      "chemical potential":', dStdGibbsEnergy(i)*dIdealConstant*dTemperature, ','

        ! Calculate driving force
        dDriving = 0D0
        do j = 1, nElements
            dDriving = dDriving + dElementPotential(j) * dStoichSpecies(i,j)
        end do
        dDriving = (dStdGibbsEnergy(i) - dDriving) / dSpeciesTotalAtoms(i)

        ! Write driving force
        write(1,'(A24,ES25.16E3,A1)') '      "driving force": ', dDriving, ','

        write(1,*) '      "stoichiometry": [', (dStoichSpecies(i,c), ',', c = 1,nElements-1), &
                                                    dStoichSpecies(i,nElements), '],'
        write(1,*) '      "elements": {'
        dTotalElements = 0D0
        nElectron = 0
        do c = 1, nElements
            dTotalElements = dTotalElements + dStoichSpecies(i,c)
            if (cElementName(c) == 'e-') nElectron = nElectron + 1
        end do
        do c = 1, nElements - nElectron
            dCurrentElement = dStoichSpecies(i,c)
            write(1,*) '        "', TRIM(cElementName(c)), '": {'
            write(1,'(A39,ES25.16E3,A1)') '          "moles of element in phase":', dCurrentElement*dTempMolesPhase, ','
            write(1,*) '          "mole fraction of phase by element":', dCurrentElement / dTotalElements, ','
            write(1,*) '          "mole fraction of element by phase":', dCurrentElement*dTempMolesPhase / dMolesElement(c)
            if (c < nElements - nElectron) then
                write(1,*) '        },'
            else
                write(1,*) '        }'
            end if
        end do
        write(1,*) '      }'
        if (i < nSpeciesPhase(nSolnPhasesSys) + nConPhasesSys) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if
    end do LOOP_PureStable

    write(1,*) '  },'

end subroutine WriteJSONPureConPhase

subroutine WriteJSONMQM(iSolnIndex)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, m, c, a, b, x, y
    integer :: iSolnIndex, iSPI, nPhaseElements
    integer :: iFirst, iLast, nSub1, nSub2, iMax
    real(8) :: dSum, dMax, dSumElementQuads, dSumElementPairs, dMolesPairs, dTempMolesPhase
    real(8) :: dZa, dZb, dZx, dZy
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi
    real(8), allocatable, dimension(:,:) :: dXij, dNij
    ! X_ij/kl corresponds to dMolFraction

    ! Only proceed if the correct phase type is selected:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ')) return

    ! Define temporary variables for sake of convenience:
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    iSPI = iPhaseSublattice(iSolnIndex)
    nSub1 = nConstituentSublattice(iSPI,1)
    nSub2 = nConstituentSublattice(iSPI,2)

    ! Allocate allocatable arrays:
    if (allocated(dXi)) deallocate(dXi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(dXij)) deallocate(dXij)
    if (allocated(dNij)) deallocate(dNij)
    j = iLast - iFirst + 1
    nPhaseElements = nSub1 + nSub2
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(dXij(nSub1,nSub2))
    allocate(dNij(nSub1,nSub2))

    ! Initialize variables:
    dXi                               = 0D0
    dYi                               = 0D0
    dNi                               = 0D0
    dXij                              = 0D0
    dNij                              = 0D0

    ! Compute X_i and Y_i
    ! Do cations first:
    dSum = 0D0
    do i = 1, nSub1
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSPI,k,1)
            dZb = dCoordinationNumber(iSPI,k,2)
            if (i == iPairID(iSPI,k,1))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZa)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2)
            end if
            if (i == iPairID(iSPI,k,2))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZb)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2)
            end if
        end do
        dSum = dSum + dNi(i)
    end do
    write(1,*) '      "cations": {'
    do i = 1, nSub1
        dXi(i) = dNi(i) / dSum
        write(1,*) '        "', TRIM(ADJUSTL(cConstituentNameSUB(iSPI,1,i))), '": {'
        write(1,*) '          "mole fraction": ', dXi(i), ','
        write(1,*) '          "charge": ', dSublatticeCharge(iSPI,1,i), ','
        write(1,*) '          "chemical group": ', iChemicalGroup(iSPI,1,i)

        if (i < nSub1) then
            write(1,*) '        },'
        else
            write(1,*) '        }'
        end if
    end do
    write(1,*) '      },'
    ! Do anions now:
    dSum = 0D0
    do i = 1, nSub2
        j = i + nSub1
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZx = dCoordinationNumber(iSPI,k,3)
            dZy = dCoordinationNumber(iSPI,k,4)
            if (j == iPairID(iSPI,k,3))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZx)
                dYi(j) = dYi(j) + (dMolFraction(l) / 2)
            end if
            if (j == iPairID(iSPI,k,4))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZy)
                dYi(j) = dYi(j) + (dMolFraction(l) / 2)
            end if
        end do
        dSum = dSum + dNi(j)
    end do
    write(1,*) '      "anions": {'
    do i = 1, nSub2
        j = i + nSub1
        dXi(j) = dNi(j) / dSum
        write(1,*) '        "', TRIM(ADJUSTL(cConstituentNameSUB(iSPI,2,i))), '": {'
        write(1,*) '          "mole fraction": ', dXi(j), ','
        write(1,*) '          "charge": ', -dSublatticeCharge(iSPI,2,i), ','
        write(1,*) '          "chemical group": ', iChemicalGroup(iSPI,2,i)
        if (i < nSub2) then
            write(1,*) '        },'
        else
            write(1,*) '        }'
        end if
    end do
    write(1,*) '      },'

    dSum = 0D0
    do m = 1, nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSPI,k,1)
            dZb = dCoordinationNumber(iSPI,k,2)
            if ((i == iPairID(iSPI,k,1)) .AND. ((j + nSub1) == iPairID(iSPI,k,3)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZa) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,1)) .AND. ((j + nSub1) == iPairID(iSPI,k,4)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZa) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,2)) .AND. ((j + nSub1) == iPairID(iSPI,k,3)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZb) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,2)) .AND. ((j + nSub1) == iPairID(iSPI,k,4)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZb) / dConstituentCoefficients(iSPI,m,1)
            end if
        end do
        dSum = dSum + dNij(i,j)
    end do

    ! Use most abundant element in phase to normalize
    dMax = 0D0
    do i = 1, nElements
        dSumElementQuads = 0D0
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dSumElementQuads = dSumElementQuads + dStoichSpecies(l,i)*dMolesSpecies(l)
        end do
        if (dSumElementQuads > dMax) then
            dMax = dSumElementQuads
            iMax = i
        end if
    end do
    dSumElementQuads = dMax

    if (dMax == 0) iMax = 1
    dSumElementPairs = 0D0
    do m = 1, nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        dSumElementPairs = dSumElementPairs + dStoichPairs(iSPI,m,iMax)*dNij(i,j)
    end do

    dMolesPairs = dSum*dSumElementQuads/dSumElementPairs
    write(1,*) '      "moles of endmembers": ', dMolesPairs, ','
    write(1,*) '      "endmembers": {'
    do m = 1, nPairsSRO(iSPI,1)
        write(1,*) '        "', TRIM(ADJUSTL(cPairName(iSPI,m))), '": {'
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        dXij(i,j) = dNij(i,j) / dSum
        write(1,*) '          "mole fraction": ', dXij(i,j), ','
        write(1,*) '          "constituents": [ "', TRIM(ADJUSTL(cConstituentNameSUB(iSPI,1,i))), '", "', &
                                                    TRIM(ADJUSTL(cConstituentNameSUB(iSPI,2,j))), '" ],'
        write(1,*) '          "stoichiometric coefficients": [', dConstituentCoefficients(iSPI,m,1), ',', &
                                                                 dConstituentCoefficients(iSPI,m,2), ']'
        if (m < nPairsSRO(iSPI,1)) then
            write(1,*) '          },'
        else
            write(1,*) '          }'
        end if
    end do
    write(1,*) '      },'

    write(1,*) '      "quadruplets": {'
    l = 0
    do k = 1, nElements
        if (-iAssemblage(k) == iSolnIndex) l = k
    end do
    if (l > 0) then
        dTempMolesPhase = dMolesPhase(l)
    else
        dTempMolesPhase = 0D0
    end if
    ! Print species:
    do i = iFirst, iLast
        k = i + 1 - iFirst
        a = iPairID(iSPI, k, 1)
        b = iPairID(iSPI, k, 2)
        x = iPairID(iSPI, k, 3) - nConstituentSublattice(iSPI,1)
        y = iPairID(iSPI, k, 4) - nConstituentSublattice(iSPI,1)
        write(1,*) '        "', TRIM(ADJUSTL(cSpeciesName(i))), '": {'
        write(1,*) '          "mole fraction":', dMolFraction(i), ","
        write(1,'(A19,ES25.16E3,A1)') '          "moles":', dMolFraction(i)*dTempMolesPhase, ","
        write(1,*) '          "chemical potential":', dChemicalPotential(i)*dIdealConstant*dTemperature, ','
        write(1,*) '          "constituents": [ "', TRIM(ADJUSTL(cConstituentNameSUB(iSPI,1,a))), '", "', &
                                                    TRIM(ADJUSTL(cConstituentNameSUB(iSPI,1,b))), '", "', &
                                                    TRIM(ADJUSTL(cConstituentNameSUB(iSPI,2,x))), '", "', &
                                                    TRIM(ADJUSTL(cConstituentNameSUB(iSPI,2,y))), '" ],'
        write(1,*) '          "coordination numbers": [', (dCoordinationNumber(iSPI,k,c), ',', c = 1,3), &
                                                    dCoordinationNumber(iSPI,k,4), '],'
        write(1,*) '          "stoichiometry": [', (dStoichSpecies(i,c), ',', c = 1,nElements-1), &
                                                    dStoichSpecies(i,nElements), ']'
        if (i < iLast) then
            write(1,*) '        },'
        else
            write(1,*) '        }'
        end if
    end do
    write(1,*) '      },'

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,dXij,dNij)

end subroutine WriteJSONMQM
