

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBI.f90
    !> \brief   Parse the data block section corresponding to a SUBL phase of a ChemSage data-file.
    !> \author  M. Poschmann
    !> \date    Jan. 19, 2021
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
    !   01/19/2021      M. POschmann     Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "SUBI" phase (Ionic Liquid Model).
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


subroutine ParseCSDataBlockSUBI(i)

    USE ModuleParseCS

    implicit none

    integer                   :: i, j, k, n, s, iDummy


    ! Initialize variables:
    n = 0

    ! Two sublattices per phase:
    nSublatticePhaseCS(nCountSublatticeCS) = 2

    ! Report an error if the number of sublattices is out of range:
    if ((nSublatticePhaseCS(nCountSublatticeCS) > nMaxSublatticeCS).OR. &
        (nSublatticePhaseCS(nCountSublatticeCS) <= 0)) then
        INFO = 33
        return
    end if

    ! Read in the number of constituents for each sublattice:
    read (1,*,IOSTAT = INFO) nConstituentSublatticeCS(nCountSublatticeCS,1:nSublatticePhaseCS(nCountSublatticeCS))
    ! Copy to nSublatticeElementsCS to use the charge array
    ! NOTE: This is bad, these two arrays should be the same thing, but got implemented
    ! separately for different models
    nSublatticeElementsCS(nCountSublatticeCS,1:nSublatticePhaseCS(nCountSublatticeCS)) = &
      nConstituentSublatticeCS(nCountSublatticeCS,1:nSublatticePhaseCS(nCountSublatticeCS))

    ! Read in the name of each constituent for each sublattice:
    do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
        read (1,*,IOSTAT = INFO) cConstituentNameSUBCS(nCountSublatticeCS,s,1:nConstituentSublatticeCS(nCountSublatticeCS,s))
    end do

    ! Record an error if necessary:
    if (INFO /= 0) then
        INFO = 2100 + i
        return
    end if

    ! I think this is # of sublattices (always 2), but need that earlier and no reason to overwrite
    read (1,*,IOSTAT = INFO) iDummy

    ! Read charges of constituents
    do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
        read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCountSublatticeCS,s,1:nConstituentSublatticeCS(nCountSublatticeCS,s))
    end do

    ! Read in the constituent indices for each component on each sublattice:
    do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
        ! Number of components for this phase:
        n = nSpeciesPhaseCS(i) - nSpeciesPhaseCS(i-1)

        read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublatticeCS, s, 1:n)
    end do

    ! Record an error and return if necessary:
    if (INFO /= 0) then
        INFO = 2200 + i
        return
    end if

    ! Loop through excess mixing parameters:
    LOOP_ExcessMixingSUBI: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the section of mixing terms is labelled "0".
        if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingSUBI

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

        ! Loop through number of mixing terms per component array:
        do k = 2, n

            ! Update counter of the number of parameters:
            nParamCS = nParamCS + 1

            ! Add additional terms:
            iRegularParamCS(nParamCS,1:nParamMax*2+1) = iRegularParamCS(nParamCS-1,1:nParamMax*2+1)

            ! Correct the indexing scheme (order of mixing parameter):
            iRegularParamCS(nParamCS, j+2) = k - 1

            ! Read mixing terms:
            read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)

        end do
    end do LOOP_ExcessMixingSUBI

    ! Report an error and return if necessary:
    if (INFO /= 0) then
        INFO = 2400 + i
        return
    end if

    return

end subroutine ParseCSDataBlockSUBI
