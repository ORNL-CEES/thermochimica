

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBL.f90
    !> \brief   Parse the data block section corresponding to a SUBL phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Dec. 21, 2012
    !> \todo    Add capability to parse magnetic mixing terms in a SUBLM phase.
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
    !> \param[in] nCountSublattice  An integer scalar representing a cummulative count of the number of phases
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


subroutine ParseCSDataBlockSUBL(i, nCountSublattice)

    USE ModuleParseCS

    implicit none

    integer                   :: i, j, k, l, m, n, p, s, nCountSublattice
    character(8),dimension(3) :: cDummyVec


    ! Initialize variables:
    n = 0

    ! Read in the number of sublattices per phase:
    read (1,*,IOSTAT = INFO) nSublatticePhaseCS(nCountSublattice)

    ! Report an error if the number of sublattices is out of range:
    if ((nSublatticePhaseCS(nCountSublattice) > nMaxSublatticeCS).OR. &
        (nSublatticePhaseCS(nCountSublattice) <= 0)) then
        INFO = 33
        return
    end if

    ! Read in the stoichiometry coefficients for each sublattice:
    read (1,*,IOSTAT = INFO) dStoichSublatticeCS(nCountSublattice,1:nSublatticePhaseCS(nCountSublattice))

    ! Read in the number of constituents for each sublattice:
    read (1,*,IOSTAT = INFO) nConstituentSublatticeCS(nCountSublattice,1:nSublatticePhaseCS(nCountSublattice))

    ! Read in the name of each constituent for each sublattice:
    LOOP_SUBL_CONST_NAME: do s = 1, nSublatticePhaseCS(nCountSublattice)

        ! Number of constituents on last line (per sublattice):
        k = MOD(nConstituentSublatticeCS(nCountSublattice,s),3)

        ! Number of full lines of constituents (per sublattice):
        l = (nConstituentSublatticeCS(nCountSublattice,s) - k) / 3

        ! Loop through full lines of constituent names:
        do m = 1, l
            read (1,*,IOSTAT = INFO) cDummyVec(1:3)
            p = (m-1)*3 + 1
            cConstituentNameSUBCS(nCountSublattice, s, p:p+2) = cDummyVec(1:3)
        end do

        ! Read in the last line of constituent names if there is less than three constituents:
        if (k /= 0) then
            read (1,*,IOSTAT = INFO) cDummyVec(1:k)
            p = l*3 + 1
            cConstituentNameSUBCS(nCountSublattice, s, p:p+k-1) = cDummyVec(1:k)
        end if

    end do LOOP_SUBL_CONST_NAME

    ! Record an error if necessary:
    if (INFO /= 0) then
        INFO = 2100 + i
        return
    end if

    ! Read in the constituent indices for each component on each sublattice:
    LOOP_SUBL_CONST_ID: do s = 1, nSublatticePhaseCS(nCountSublattice)

        ! Number of components for this phase:
        n = nSpeciesPhaseCS(i) - nSpeciesPhaseCS(i-1)

        ! Number of consituents on last line (per sublattice):
        k = MOD(n, 19)

        ! Number of full lines of constituent indices (per sublattice):
        l = (n - k) / 19

        ! Loop through full lines of consitutent indices:
        do m = 1, l
            p = (m-1)*19 + 1
            read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublattice, s, p:p+18)
        end do

        ! Read in the last line of the constituent indices if there is less than 19 constituents:
        if (k /= 0) then
            p = (l)*19 + 1
            read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublattice, s, p:p+k-1)
        end if

    end do LOOP_SUBL_CONST_ID

    ! Record an error and return if necessary:
    if (INFO /= 0) then
        INFO = 2200 + i
        return
    end if

    ! SUBLM phases include magnetic mixing terms right before non-ideal mixing terms.
    ! NOTE: I HAVE NOT DEVELOPED THE CAPABILITY TO READ IN THESE MAGNETIC MIXING TERMS
    ! The end of the list of mixing terms is indicated by a "0".
    if (cSolnPhaseTypeCS(i) == 'SUBLM') then

        read (1,*,IOSTAT = INFO) j

        ! Record an error:
        if ((INFO /= 0).OR.(j /= 0)) then
            INFO = 2300 + i
            return
        end if
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

        ! Store number of mixing paramters per component array:
        n = iRegularParamCS(nParamCS, j+2)

        ! Correct the indexing scheme (order of mixing) for the first parameter:
        iRegularParamCS(nParamCS, j+2) = 0

        ! Loop through number of mixing terms per component array:
        LOOP_ParamArray: do k = 2, n

            ! Update counter of the number of parameters:
            nParamCS = nParamCS + 1

            ! Add additional terms:
            iRegularParamCS(nParamCS,1:nParamMax*2+1) = iRegularParamCS(nParamCS-1,1:nParamMax*2+1)

            ! Correct the indexing scheme (order of mixing parameter):
            iRegularParamCS(nParamCS, j+2) = k - 1

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
