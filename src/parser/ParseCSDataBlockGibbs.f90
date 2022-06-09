

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockGibbs.f90
    !> \brief   Parse the coefficients of the Gibbs energy equations in the datablock section of a ChemSage
    !!           data-file.
    !> \author  M.H.A. Piro
    !> \date    Apr. 3, 2018
    !> \sa      ParseCSDataBlock.f90
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
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/06/2011      M.H.A. Piro         Original code
    !   10/26/2011      M.H.A. Piro         Allowed for a formatted character string when reading in a species name
    !   03/01/2012      M.H.A. Piro         Corrected array bound mismatch error in reading stocihiometry
    !                                       coefficients.
    !   06/30/2016      M.H.A. Piro         Allowed for data-files containing species with type "1" Gibbs energy
    !                                        equations to be parsed.  Careful: I don't know if this is reliable.
    !   10/12/2016      M.H.A. Piro         A correction was applied when parsing type 1 species.  Type 4 and 16
    !                                        equations end with the possibility of additional terms, whereby no
    !                                        additional terms is 0, but type 1 equations do not.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the coefficients of the Gibbs energy equations in
    !! the "data block" section of a ChemSage data-file and store the coefficients for computation
    !! in another program.  Details about the INFO error codes are given in ParseCSDataFile.f90.
    !! Details about entry numbers are explained in the "ChemApp Programmer's manual".
    !!
    !! Chemical species may include solution species, pure separate phases and what I call "dummy
    !! species".  Dummy species do not have any physical significance and are only included in ChemSage
    !! datafiles for numerical purposes.  Dummy species are always pure species (i.e., U, O or Pu;
    !! as opposed to UO2 and Pu5O8) and are included when the database does not include any pure
    !! species. Dummy species are identified in ChemSage data-files as species that follow pure
    !! separate phases that do not end their character string with ")".  All pure separate phases will
    !! end their character string with "(s)".
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   i                   An integer scalar representing the phase index.
    !> \param[in]   j                   An integer scalar representing the species index.
    !> \param[in]   iCounterGibbsEqn    An integer scalar representing the index of the Gibbs energy equation.
    !
    ! INFO                              A scalar integer that indicates a successful exit or identifies an error.
    ! nGibbsEqSpecies                   Number of Gibbs energy equations for a particular species.
    ! iSpeciesAtomsCS                   Integer matrix representing the number of atoms of a particular elements
    !                                    in a species (i.e., stoichiometry matrix) from the data-file.
    ! iDummy                            An integer scalar dummy variable.
    ! cSpeciesNameCS                    A character vector representing the name of a species from the data-file.
    ! dGibbsCoeffSpeciesTemp            Temporary double array of coefficients for a Gibbs energy equation.
    ! dGibbsMagneticCS                  A double real matrix representing magnetic contributions to the molar
    !                                    Gibbs energy term (nSpecies,4).  Second dimension: 1) Critical
    !                                    temperature; 2) B, 3) Structure factor; 4 p.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockGibbs(i,j,iCounterGibbsEqn)

    USE ModuleParseCS

    implicit none

    integer               :: i, j, k, l, m, iDummy, iCounterGibbsEqn, iGibbsEqType
    real(8)               :: dTemp
    real(8),dimension(15) :: dTempVec
    real(8),parameter     :: dGibbsDummy = 1D6
    character(5)          :: cDummy
    character(26)         :: cSpeciesNameDummy
    character(70)         :: cSpeciesNameLine, cSubstring


    ! Initialize variables:
    k = 0
    l = 0
    dTempVec = 0D0

    ! Entry 3: Read name of constituent species:
    read (1,'(a)',IOSTAT = INFO) cSpeciesNameLine
    cSpeciesNameLine = TRIM(ADJUSTL(cSpeciesNameLine))
    k = 26
    cSpeciesNameCS(j) = TRIM(ADJUSTL(cSpeciesNameLine(:k)))

    ! Create substring of remaining
    cSubstring = TRIM(ADJUSTL(cSpeciesNameLine(k:)))
    k = SCAN(cSubstring, ' ')
    if (k > 1) then
      ! Ignore non-numerics on this line
        l = SCAN(cSubstring, '-.0123456789')
        if (l == 1) then
            select case (cSolnPhaseTypeCS(i))
            case ('SUBG','SUBQ','SUBL','SUBI')
                print *, 'Max_x used in multi-sublattice phase ', cSolnPhaseNameCS(i), ' but not yet implemented.'
                print *, 'Max_x penalty term will be ignored.'
            case default
                read(cSubstring(:k-1), *) dMaxXCS(j)
                read(cSubstring(k:), *) dPenaltyXCS(j)
          end select
        end if
    end if

    ! Check to see if there is more than one particle per constituent formula mass.
    k = SCAN(cSpeciesNameCS(j), ':', back = .TRUE.)

    ! Record the number of particles per constituent formula mass.
    if (k /= 0) then
        cSpeciesNameDummy = cSpeciesNameCS(j)
        cDummy = cSpeciesNameDummy(k-1:k)

        if (cDummy(1:2) == '):') then
            cDummy = cSpeciesNameDummy(k+1:k+2)
            read( cDummy, '(I1)' )  iParticlesPerMoleCS(j)
        end if
    end if

    ! Check for a dummy species (a description of a "dummy species" is given in the header).
    l = SCAN(cSpeciesNameCS(j),'#', back = .TRUE.)

    ! The phase index of a dummy species is set to -1:
    if ((i == 0).AND.(l > 0)) then
        ! If it contains '#', assume a dummy
        iPhaseCS(j) = -1
    end if

    if (INFO /= 0) then
        INFO = 1300 + i
        return
    end if

    ! Entry 4: Read thermodynamic data for constituent species:
    read (1,*, IOSTAT = INFO) iGibbsEqType, nGibbsEqSpecies(j), dStoichSpeciesCS(j,1:nElementsCS)

    if (.NOT.((iGibbsEqType == 4).OR.(iGibbsEqType == 16).OR.(iGibbsEqType == 1).OR.(iGibbsEqType == 13))) then
        ! The type of Gibbs energy equation is not supported.
        INFO = 1400 + i
        return
    end if

    ! Note: ChemSage data-files store the stoichiometry coefficients as real variables, which may
    ! be less than one when there are more than 1 particles/mole.  I am going to convert the stoichiometry
    ! to an integer to minimize numerical error while still recording the number of particles per mole.
    dTemp = DFLOAT(iParticlesPerMoleCS(j))
    dStoichSpeciesCS(j,1:nElementsCS) = dStoichSpeciesCS(j,1:nElementsCS) * dTemp

    if (INFO /= 0) then
        INFO = 1400 + i
        return
    end if

    LOOP_GibbsEquations: do m = 1, nGibbsEqSpecies(j)

        ! Entry 5: Read upper temperature limit and the coefficients for the Gibbs energy equation for this species:
        iCounterGibbsEqn = iCounterGibbsEqn + 1

        if (iCounterGibbsEqn > SIZE(dGibbsCoeffSpeciesTemp, DIM = 2)) then
            INFO = 5
            return
        end if

        ! Read in the standard Gibbs energy terms:
        read (1,*, IOSTAT = INFO) dGibbsCoeffSpeciesTemp(1:7,iCounterGibbsEqn)

        l = 0
        if ((iGibbsEqType == 16).OR.(iGibbsEqType == 4)) then
            ! Note that Gibbs energy equations designated type 4 and 16 have an additional
            ! line, but type 1 does not.

            ! Read in the number of additional terms to the standard Gibbs energy equation:
            read (1,*,IOSTAT = INFO) l

            ! Go back to entry 5:
            backspace(UNIT = 1)
       end if

        if (l == 1) then
            read (1,*, IOSTAT = INFO) iDummy, dGibbsCoeffSpeciesTemp(8:9, iCounterGibbsEqn)
        elseif (l == 2) then
            read (1,*, IOSTAT = INFO) iDummy, dGibbsCoeffSpeciesTemp(8:11,iCounterGibbsEqn)
        elseif (l == 3) then
            read (1,*, IOSTAT = INFO) iDummy, dGibbsCoeffSpeciesTemp(8:13,iCounterGibbsEqn)
        elseif (l == 0) then
            ! Do nothing... except still have to read a line for 4 or 16.
            if ((iGibbsEqType == 16).OR.(iGibbsEqType == 4)) then
                read (1,*,IOSTAT = INFO) l
            end if
        else
            INFO = 1500 + i
            return
        end if

    end do LOOP_GibbsEquations      ! End loop of variable m

    ! Check if the equation has magnetic contributions:
    if ((iGibbsEqType == 16).OR.(iGibbsEqType == 13)) then
        ! This species has magnetic contributions to the Gibbs energy equation:

        ! Check if this is a pure condensed phase or solution phase:
        if (i == 0) then
            ! Read all four coefficients for a pure condensed phase:
            read (1,*,IOSTAT = INFO) dGibbsMagneticCS(j,1:4)
            !nMagneticTermsCS = nMagneticTermsCS + 1
        else
            ! Read the first two coefficients for a solution phase:
            read (1,*,IOSTAT = INFO) dGibbsMagneticCS(j,1:2)
        end if

        ! Report an error if necessary:
        if (INFO /= 0) then
            INFO = 1500 + j
            return
        end if
    end if

    ! Assign an arbitrarily large positive value for the Gibbs energy equations for dummy species:
    if (iPhaseCS(j) == -1) then
        dGibbsCoeffSpeciesTemp(2,iCounterGibbsEqn) = dGibbsDummy
    end if

    return

end subroutine ParseCSDataBlockGibbs
