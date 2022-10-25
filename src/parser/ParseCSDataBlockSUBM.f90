

!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBM.f90
    !> \brief   Parse the data block section corresponding to a SUBM phase of a ChemSage data-file.
    !> \author  M. Poschmann
    !> \date    Oct. 21, 2022
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \sa      ParseCSInterpolationOverrides.f90
    !
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   11/21/2022      M. Poschmann    Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "SUBM" phase, which a two-sublattice model.
    !!
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! nSpeciesCS                Number of species in the system (combined solution species and pure
    !                            separate phases).
    ! nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    ! iSpeciesAtomsCS           Integer matrix representing the number of atoms of a particular
    !                            elements in a species.
    ! iParticlesPerMoleCS       An integer vector containing the number of particles per mole of the
    !                            constituent species formula mass.  The default value is 1.
    ! cSolnPhaseNameCS          The name of a solution phase.
    ! cSolnPhaseTypeCS          The type of a solution phase.
    ! cSolnPhaseTypeSupport     A character array representing solution phase types that are supported.
    ! iRegularParamCS           An integer matrix representing the parameter index for the first dimension
    !                            and the mixing terms on the second dimension.  For the second dimension, the
    !                            first coefficient indicates whether the parameter is a binary or ternary term (n),
    !                            the next n coefficients correspond to the constituent indices, and the last
    !                            coefficient corresponds to the exponent.
    !
!-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockSUBM( i )

    USE ModuleParseCS

    implicit none

    integer                     :: i, j, k, nPairs, nCSCS, nTotalConst, c, s, iTemp
    character(75) :: cTempConstituent

    nCSCS = nCountSublatticeCS

    ! This line contains N integers (where N is the number of sublattices)
    ! where each integer represents the number of constituents on the respective
    ! sublattice. There are always two sublattices for SUBG phases.
    read (1,*,IOSTAT = INFO) nConstituentSublatticeCS(nCSCS,1:2)
    nConstituentSublatticeCS(nCSCS,1:2) = nConstituentSublatticeCS(nCSCS,1:2)
    nSublatticePhaseCS(nCSCS) = 2
    nTotalConst = nConstituentSublatticeCS(nCSCS,1)+nConstituentSublatticeCS(nCSCS,2)
    allocate(dStoichConstituentCS(nTotalConst,nElementsCS))
    dStoichConstituentCS = 0D0

    nPairs = nConstituentSublatticeCS(nCSCS,1) * nConstituentSublatticeCS(nCSCS,2)

    ! Read in names of constituents:
    do s = 1, 2
        do c = 1, CEILING(REAL(nConstituentSublatticeCS(nCSCS,s))/3D0)
            read (1,112,IOSTAT = INFO) cTempConstituent
            112 FORMAT (A75)
            do j = 1, 3
                if ((3*(c-1)+j) <= nConstituentSublatticeCS(nCSCS,s)) then
                    cConstituentNameSUBCS(nCSCS,s,3*(c-1)+j) = TRIM(ADJUSTL(cTempConstituent((j-1)*25+1:j*25)))
                end if
            end do
        end do
    end do

    do j = 1, 2
        do k = 1, nConstituentSublatticeCS(nCSCS,j)
            cConstituentNameSUBCS(nCSCS,j,k) = TRIM(ADJUSTL(cConstituentNameSUBCS(nCSCS,j,k)))
        end do
    end do

    ! Read in the charge of each constituent on the first sublattice.
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCSCS,1,1:nConstituentSublatticeCS(nCSCS,1))

    ! Chemical groups on sublattice 1:
    read (1,*,IOSTAT = INFO) iChemicalGroupCS(nCSCS,1,1:nConstituentSublatticeCS(nCSCS,1))

    ! Read in the charge of each constituent on the second sublattice.
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCSCS,2,1:nConstituentSublatticeCS(nCSCS,2))

    ! Chemical groups on sublattice 2:
    read (1,*,IOSTAT = INFO) iChemicalGroupCS(nCSCS,2,1:nConstituentSublatticeCS(nCSCS,2))

    ! This entry appears to represent the IDs matching constituents on the first sublattice to species:
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCSCS, 1, 1:nPairs)

    ! This entry appears to represent the IDs matching constituents on the second sublattice to species:
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCSCS, 2, 1:nPairs)

    ! Loop through excess mixing parameters:
    j = 0
    LOOP_ExcessMixingSUBM: do
        j = j + 1
        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the parameter listing is marked by "0",
        ! or a negative number indicating the number of extra parameter lines.
        ! These lines indicate overwriting of default interpolation schemes.
        if (iRegularParamCS(nParamCS+1,1) <= 0) then
            call ParseCSInterpolationOverrides(i)
            exit LOOP_ExcessMixingSUBM
        end if

        ! Check if the parameter is binary or ternary:
        if ((iRegularParamCS(nParamCS+1,1) == 3) .OR. (iRegularParamCS(nParamCS+1,1) == 4)) then

            ! Count the number of parameters:
            nParamCS = nParamCS + 1

            ! Read in the list of constituent indices and parameters involved in this excess term:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:iRegularParamCS(nParamCS,1)*2+1), dRegularParamCS(nParamCS,1:6)

            ! Reorder parameters to put constituents contiguously
            iTemp = iRegularParamCS(nParamCS,iRegularParamCS(nParamCS,1)*2+1)
            do k = iRegularParamCS(nParamCS,1)*2+1, iRegularParamCS(nParamCS,1)+2, -1
                iRegularParamCS(nParamCS,k) = iRegularParamCS(nParamCS,k-1)
            end do
            iRegularParamCS(nParamCS,iRegularParamCS(nParamCS,1)+1) = iTemp


        else
            !! This parameter is not recognized; record an error.
            INFO = 10000 + 1000*j + i
            return
        end if

    end do LOOP_ExcessMixingSUBM

    ! Report an error if necessary:
    if (INFO /= 0) INFO = 1600 + i

    return

end subroutine ParseCSDataBlockSUBM
