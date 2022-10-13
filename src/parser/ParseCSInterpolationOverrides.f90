
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSInterpolationOverrides.f90
    !> \brief   Parse overrides of ternary interpolation schemes from a ChemSage data-file.
    !> \author  M. Poschmann
    !> \date    Oct. 13, 2022
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \sa      ParseCSDataBlockQKTO.f90
    !> \sa      ParseCSDataBlockSUBG.f90
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
    !> \details While ternary interpolation schemes are usually determined by rules based on the chemical
    !! groups of the constituents involved, database developers can manually override these using some
    !! cryptic lines appended after the excess mixing parameters. Those lines are parsed here.
    !!
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nInterpolationOverrideCS  Integer array (per phase) indicating number of override lines.
    ! iInterpolationOverrideCS  Integer array (per phase, override line) containing constituent indices.
    !! Order is: active sublattice indices 1-3 (ternary), inactive sublattice index (0 for non-sublattice
    !! phases), asymmetric constituent index (0 if none).
    !
!-------------------------------------------------------------------------------------------------------------

subroutine ParseCSInterpolationOverrides( i )

    USE ModuleParseCS

    implicit none

    integer :: i, j, k, l, iK, iT, nIO
    character(8),dimension(3)  :: cTempVec
    character(8)  :: cTemp
    integer, dimension(6)  :: iTemp
    integer, dimension(2)  :: iToop
    integer                :: nToop

    ! Loop over interpolation scheme overrides
    do k = 1, -iRegularParamCS(nParamCS+1,1)
        nInterpolationOverrideCS(i) = nInterpolationOverrideCS(i) + 1
        nIO = nInterpolationOverrideCS(i)
        ! Have to read interpolation override line, unfortunately the formatting is a mess.
        ! Specifically, strings and integers abut.
        read (1,*,IOSTAT = INFO) iInterpolationOverrideCS(i,nIO,1), iInterpolationOverrideCS(i,nIO,2), cTempVec(1), & 
        iTemp(1), cTempVec(2), iTemp(2), cTempVec(3), iTemp(3), iTemp(4), iInterpolationOverrideCS(i,nIO,4)

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
                    READ(cTemp(1:iK-1),*) iInterpolationOverrideCS(i,nIO,3)
                else
                    READ(cTemp(1:iK-1),*) iTemp(j+3)
                end if
            else if (iT > 0) then
                nToop = nToop + 1
                if (j == 1) then
                    READ(cTemp(1:iT-1),*) iInterpolationOverrideCS(i,nIO,3)
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
                    if (iInterpolationOverrideCS(i,nIO,j) == iToop(l)) then
                        ! Constituent NON-CONSTANT, keep going
                        cycle LOOP_checkAsym
                    end if
                end do
                ! If we get here, it is the constant constituent. Save for calculation.
                iInterpolationOverrideCS(i,nIO,5) = iInterpolationOverrideCS(i,nIO,j)
                exit LOOP_checkAsym
            end do LOOP_checkAsym
        else if (nToop == 0) then
            ! Nothing, this is symmetric.
            ! Leave iInterpolationOverrideCS(i,5) = 0 to indicate symmetric.
        else
            ! This shouldn't be possible
            nInterpolationOverrideCS(i) = nInterpolationOverrideCS(i) - 1
            iInterpolationOverrideCS(i,nIO,:) = 0
            print *, 'Interpolation override not understood, ignored.'
            print *,  iInterpolationOverrideCS(i,nIO,1), iInterpolationOverrideCS(i,nIO,2), cTempVec(1), & 
            iTemp(1), cTempVec(2), iTemp(2), cTempVec(3), iTemp(3), iTemp(4), iInterpolationOverrideCS(i,nIO,4)
        end if
    end do

end subroutine ParseCSInterpolationOverrides