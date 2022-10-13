

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockQKTO.f90
    !> \brief   Parse the data block section corresponding to a QKTO phase of a ChemSage data-file.
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
    !   10/06/2011      M.H.A. Piro     Original code
    !   12/21/2012      M.H.A. Piro     Relocate to independent f90 file.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "QKTO" phase (Quasi-chemical Kohler-TOop).  The molar excess Gibbs energy of mixing is
    !! generally given as \f$ \Delta g_i^{ex} = \sum_{p=1} x_i^m x_j^n (^pL) \f$, where i and j are the
    !! constituent indices, m and n are the exponents (provided by the model) and \f$ ^pL \f$ is the
    !! temperature dependent mixing parameter.
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
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockQKTO( i )

    USE ModuleParseCS

    implicit none

    integer :: i, j, k, l, iK, iT
    character(8),dimension(3)  :: cTempVec
    character(8)  :: cTemp
    integer, dimension(6)  :: iTemp
    integer, dimension(2)  :: iToop
    integer                :: nToop


    ! Loop through excess parameters:
    LOOP_ExcessMixingQKTO: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingQKTO
        ! The end of the parameter listing is marked by "0"
        ! or a negative number indicating the number of extra parameter lines.
        ! These lines indicate overwriting of default interpolation schemes.
        ! Not implemented yet.
        if (iRegularParamCS(nParamCS+1,1) <= 0) then
            do k = 1, -iRegularParamCS(nParamCS+1,1)
                nInterpolationOverrideCS(i) = nInterpolationOverrideCS(i) + 1
                ! Have to read interpolation override line, unfortunately the formatting is a mess.
                ! Specifically, strings and integers abut.
                read (1,*,IOSTAT = INFO) iInterpolationOverrideCS(i,1), iInterpolationOverrideCS(i,2), cTempVec(1), & 
                iTemp(1), cTempVec(2), iTemp(2), cTempVec(3), iTemp(3), iTemp(4), iInterpolationOverrideCS(i,4)

                nToop = 0
                ! Loop over strings and break them into string and integer parts
                do j = 1, 3
                    cTemp = cTempVec(j)
                    ! First check if Kohler
                    iK = INDEX(cTemp,'K')
                    ! Then if Toop
                    iT = INDEX(cTemp,'T')
                    if (iK > 0) then
                        if (j == 1) then
                            READ(cTemp(1:iK-1),*) iInterpolationOverrideCS(i,3)
                        else
                            READ(cTemp(1:iK-1),*) iTemp(j+3)
                        end if
                    else if (iT > 0) then
                        nToop = nToop + 1
                        if (j == 1) then
                            READ(cTemp(1:iT-1),*) iInterpolationOverrideCS(i,3)
                        else
                            READ(cTemp(1:iT-1),*) iTemp(j+3)
                        end if
                        ! Save constituents identified with Toop label
                        READ(cTemp(iT+1:8),*) iToop(nToop)
                    end if
                end do

                ! An asymmetric triad will have 2 Toop entries.
                ! Frustratingly, these identify the NON-CONSTANT constituents.
                ! So, back out the constant one here.
                if (nToop == 2) then
                    LOOP_checkAsym: do j = 1, 3
                        do l = 1, 2
                            if (iInterpolationOverrideCS(i,j) == iToop(l)) then
                                ! Constituent NON-CONSTANT, keep going
                                cycle LOOP_checkAsym
                            end if
                        end do
                        ! If we get here, it is the constant constituent. Save for calculation.
                        iInterpolationOverrideCS(i,5) = iInterpolationOverrideCS(i,j)
                        exit LOOP_checkAsym
                    end do LOOP_checkAsym
                else if (nToop == 0) then
                    ! Nothing, this is symmetric.
                    ! Leave iInterpolationOverrideCS(i,5) = 0 to indicate symmetric.
                else
                    ! This shouldn't be possible
                    nInterpolationOverrideCS(i) = nInterpolationOverrideCS(i) - 1
                    iInterpolationOverrideCS(i,:) = 0
                    print *, 'Interpolation override not understood, ignored.'
                    print *,  iInterpolationOverrideCS(i,1), iInterpolationOverrideCS(i,2), cTempVec(1), & 
                    iTemp(1), cTempVec(2), iTemp(2), cTempVec(3), iTemp(3), iTemp(4), iInterpolationOverrideCS(i,4)
                end if
            end do
            exit LOOP_ExcessMixingQKTO
        end if

        nParamCS = nParamCS + 1

        if (iRegularParamCS(nParamCS,1) == 2) then
            ! Binary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:5), dRegularParamCS(nParamCS,1:6)
        elseif (iRegularParamCS(nParamCS,1) == 3) then
            ! Ternary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:7), dRegularParamCS(nParamCS,1:6)
        else
            ! This is not recognized.
            INFO = 1600 + i
            return
        end if

    end do LOOP_ExcessMixingQKTO

    ! Report an error if necessary:
    if (INFO /= 0) then
        INFO = 1600 + i
        return
    end if

    return

end subroutine ParseCSDataBlockQKTO
