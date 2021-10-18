

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBL.f90
    !> \brief   Parse the data block section corresponding to a SUBL phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Dec. 21, 2012
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
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
    !   12/21/2012      M.H.A. Piro     Original code
    !   02/05/2013      M.H.A. Piro     Fixed bug in parsing constituent indices with multiple lines and
    !                                    constituent names on multiple lines.
    !   02/14/2013      M.H.A. Piro     Fixed bug in parsing multiple lines for mixing parameters (happy
    !                                    valentine's day).
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "SUBL" phase (Compound Energy Formalism).
    !!
    !!
    !! The coefficients for the standard molar Gibbs energy equations for pure species originate from a ChemSage
    !! data-file that was parsed from the ParseCSDataFile program.  The format for the coefficients follow:
    !!
    !! \f$ g_i^{\circ} = A + BT + CTln(T) + DT^2 + ET^3 + F/T + (GT^U + HT^V + IT^W + JT^X + KT^Y + LT^Z) \f$
    !!
    !
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] i                 An integer scalar representing the solution phase index.
    !> \param[in] nCountSublatticeCS  An integer scalar representing a cummulative count of the number of phases
    !                            containing a sublattice.
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! cSolnPhaseNameCS          A character vector representing the name of each solution phase.
    ! cSolnPhaseTypeCS          A character vector representing the type of each solution phase.
    ! nSublatticePhaseCS        An integer vector representing the number of sublattices per phase.
    ! dStoichSublatticeCS       A double real matrix representing the stoichiometry coefficient for each
    !                            sublattice of a phase.
    ! nConstituentSublatticeCS  An integer matrix representing the number of constituents per sublattice of a
    !                            particular phase.
    ! cConstituentNameSUBCS   A character matrix representing the names of each constituent on each sublattice
    !                            for a particular phase.
    ! iConstituentSublatticeCS  An integer matrix representing the constituent indices on each sublattice.
    ! iRegularParamCS           An integer matrix representing the index of the mixing parameter and important
    !                            information for mixing.  For iRegularParamCS(i,j) for a SUBL phase, i refers to
    !                            the mixing parameter index and j refers to a vector specific to that mixing term.
    !                            The first coefficient along j refers to the number of constituents (n) involved for
    !                            that parameter, coefficients 2:n+1 represent the indices of constituents,
    !                            coefficient n+2 represents the order of the mixing parameter.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockSUBL(i)

    USE ModuleParseCS

    implicit none

    integer :: i, j, k, n, s, c
    logical :: lTripleTerm
    character(75) :: cTempConstituent

    ! Initialize variables:
    n = 0

    ! Read in the number of sublattices per phase:
    read (1,*,IOSTAT = INFO) nSublatticePhaseCS(nCountSublatticeCS)

    ! Report an error if the number of sublattices is out of range:
    if ((nSublatticePhaseCS(nCountSublatticeCS) > nMaxSublatticeCS).OR. &
        (nSublatticePhaseCS(nCountSublatticeCS) <= 0)) then
        INFO = 33
        return
    end if

    ! Read in the stoichiometry coefficients for each sublattice:
    read (1,*,IOSTAT = INFO) dStoichSublatticeCS(nCountSublatticeCS,1:nSublatticePhaseCS(nCountSublatticeCS))

    ! Read in the number of constituents for each sublattice:
    read (1,*,IOSTAT = INFO) nConstituentSublatticeCS(nCountSublatticeCS,1:nSublatticePhaseCS(nCountSublatticeCS))

    ! Read in the name of each constituent for each sublattice:
    LOOP_SUBL_CONST_NAME: do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
        do c = 1, CEILING(REAL(nConstituentSublatticeCS(nCountSublatticeCS,s))/3D0)
            read (1,112,IOSTAT = INFO) cTempConstituent
            112 FORMAT (A75)
            do j = 1, 3
                if ((3*(c-1)+j) <= nConstituentSublatticeCS(nCountSublatticeCS,s)) then
                    cConstituentNameSUBCS(nCountSublatticeCS,s,3*(c-1)+j) = TRIM(ADJUSTL(cTempConstituent((j-1)*25+1:j*25)))
                end if
            end do
        end do
    end do LOOP_SUBL_CONST_NAME

    ! Record an error if necessary:
    if (INFO /= 0) then
        INFO = 2100 + i
        return
    end if

    ! Read in the constituent indices for each component on each sublattice:
    LOOP_SUBL_CONST_ID: do s = 1, nSublatticePhaseCS(nCountSublatticeCS)

        ! Number of components for this phase:
        n = nSpeciesPhaseCS(i) - nSpeciesPhaseCS(i-1)

        read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublatticeCS, s, 1:n)

    end do LOOP_SUBL_CONST_ID

    ! Record an error and return if necessary:
    if (INFO /= 0) then
        INFO = 2200 + i
        return
    end if

    ! SUBLM phases include magnetic mixing terms right before non-ideal mixing terms.
    ! The end of the list of mixing terms is indicated by a "0".
    if (cSolnPhaseTypeCS(i) == 'SUBLM') then! call ParseCSMagneticMixing(i)
        LOOP_MagneticMixingSUBL: do

            ! Read in number of constituents involved in parameter:
            read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS+1,1)

            ! The end of the section of mixing terms is labelled "0".
            if (iMagneticParamCS(nMagParamCS+1,1) == 0) exit LOOP_MagneticMixingSUBL

            ! Update counter of the number of parameters:
            nMagParamCS = nMagParamCS + 1

            ! Read in the list of constituent indices involved in this parameter:
            j = iMagneticParamCS(nMagParamCS,1)
            read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS,2:j+2)

            ! Read in the mixing parameter:
            read (1,*,IOSTAT = INFO) dMagneticParamCS(nMagParamCS,1:2)

            ! Store number of mixing parameters per component array:
            n = iMagneticParamCS(nMagParamCS, j+2)

            ! Correct the indexing scheme (order of mixing) for the first parameter:
            iMagneticParamCS(nMagParamCS, j+2) = 0

            ! Loop through number of mixing terms per component array:
            do k = 2, n

                ! Update counter of the number of parameters:
                nMagParamCS = nMagParamCS + 1

                ! Add additional terms:
                iMagneticParamCS(nMagParamCS,1:nParamMax*2+1) = iMagneticParamCS(nMagParamCS-1,1:nParamMax*2+1)

                ! Correct the indexing scheme (order of mixing parameter):
                iMagneticParamCS(nMagParamCS, j+2) = k - 1

                ! Read mixing terms:
                read (1,*,IOSTAT = INFO) dMagneticParamCS(nMagParamCS,1:2)

            end do
        end do LOOP_MagneticMixingSUBL
    end if

    ! Loop through excess mixing parameters:
    LOOP_ExcessMixingSUBL: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the section of mixing terms is labelled "0".
        if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingSUBL

        ! Update counter of the number of parameters:
        nParamCS = nParamCS + 1

        ! Read in the list of constituent indices involved in this parameter:
        j = iRegularParamCS(nParamCS,1)
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:j+2)

        ! Read in the mixing parameter:
        read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)

        ! Store number of mixing parameters per component array:
        n = iRegularParamCS(nParamCS, j+2)

        ! Correct the indexing scheme (order of mixing) for the first parameter:
        iRegularParamCS(nParamCS, j+2) = 0

        ! lTripleTerm is false unless the specified conditions are met
        lTripleTerm = .false.
        ! Some mixing cases require three interaction parameters to be created
        ! if the term is L_0
        ! If one mixing term create three mixing terms
        if ((iRegularParamCS(nParamCS,1)+1-nSublatticePhaseCS(nCountSublatticeCS) == 3) .AND. (n == 1)) then
            ! Adds 2 mixing terms
            n = 3
            ! If lTripleTerm is true then interaction parameters will be added
            lTripleTerm = .TRUE.
        end if

        ! Loop through number of mixing terms per component array:
        LOOP_ParamArray: do k = 2, n

            ! Update counter of the number of parameters:
            nParamCS = nParamCS + 1

            ! Add additional terms:
            iRegularParamCS(nParamCS,1:nParamMax*2+1) = iRegularParamCS(nParamCS-1,1:nParamMax*2+1)

            ! Correct the indexing scheme (order of mixing parameter):
            iRegularParamCS(nParamCS, j+2) = k - 1

            ! Re-read in mixing parameters
            if (lTripleTerm) backspace(UNIT = 1)

            ! Read mixing terms:
            read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)

        end do LOOP_ParamArray
    end do LOOP_ExcessMixingSUBL

    ! Report an error and return if necessary:
    if (INFO /= 0) then
        INFO = 2400 + i
        return
    end if

    return

end subroutine ParseCSDataBlockSUBL
