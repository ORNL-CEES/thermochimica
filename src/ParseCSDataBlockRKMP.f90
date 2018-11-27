

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockRKMP.f90
    !> \brief   Parse the data block section corresponding to a RKMP phase of a ChemSage data-file.
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
    !   01/07/2013      M.H.A. Piro     Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "RKMP" phase (Redlich-Kister-Muggiano-Polynomial).  This subroutine also handles phases
    !! designated by the "RKMPM" flag, where the additional "M" indicates contributions from magnetic ordering.
    !! The molar excess Gibbs energy of mixing for this model is generally given as,
    !! \f$ \Delta g_i^{ex} = \sum_{p=1} x_i x_j \sum (x_i - x_j) (^p L_{i,j}) + x_i x_j x_k (^p L_{i,j,k}) \f$,
    !! where \f$ ^p L_{i,j} \f$ is the binary mixing parameter and \f$ ^p L_{i,j,k} \f$ is the ternary
    !! mixing parameter.
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


subroutine ParseCSDataBlockRKMP( i )

    USE ModuleParseCS

    implicit none

    integer :: i, j, k


    ! RKMPM phases include magnetic mixing terms right before non-ideal mixing terms.
    ! The end of the list of mixing terms is indicated by a "0".
    if (cSolnPhaseTypeCS(i) == 'RKMPM') then

        read (1,*,IOSTAT = INFO) j

        ! Record an error:
        if ((INFO /= 0).OR.(j /= 0)) then
            INFO = 1600 + i
            return
        end if
    end if

    ! Loop through excess mixing parameters:
    LOOP_ExcessMixingRKMP: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the parameter listing is marked by "0":
        if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingRKMP

        ! Check if the parameter is binary or ternary:
        if (iRegularParamCS(nParamCS+1,1) == 2) then

            ! Binary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,2:4)

            ! Store the number of equations to temporary memory:
            k = iRegularParamCS(nParamCS+1,4)

            ! Loop through the number of terms per constituent matrix:
            do j = 1, k
                ! Count the number of parameters:
                nParamCS = nParamCS + 1

                ! Record mixing ID's:
                iRegularParamCS(nParamCS,1:4) = iRegularParamCS(nParamCS-j+1,1:4)

                ! Set the exponent for RKMP:
                iRegularParamCS(nParamCS,4)   = j - 1
                read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)

            end do

        elseif (iRegularParamCS(nParamCS+1,1) == 3) then

            ! Ternary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,2:5)

            k = iRegularParamCS(nParamCS+1,5)

            ! Loop through the number of terms per constituent matrix:
            do j = 1, k

                ! Count the number of parameters:
                nParamCS = nParamCS + 1

                ! Record mixing ID's:
                iRegularParamCS(nParamCS,1:5) = iRegularParamCS(nParamCS-j+1,1:5)

                ! Record the last constituent index:
                iRegularParamCS(nParamCS,5)   = iRegularParamCS(nParamCS,1+j)
                read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)
            end do

        elseif (iRegularParamCS(nParamCS+1,1) == 4) then

            ! Quaternary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,2:6)

            k = iRegularParamCS(nParamCS+1,6)

            ! Loop through the number of terms per constituent matrix:
            do j = 1, k

                ! Count the number of parameters:
                nParamCS = nParamCS + 1

                ! Record mixing ID's:
                iRegularParamCS(nParamCS,1:6) = iRegularParamCS(nParamCS-j+1,1:6)

                ! Record the last constituent index:
                iRegularParamCS(nParamCS,6)   = iRegularParamCS(nParamCS,1+j)
                read (1,*,IOSTAT = INFO) dRegularParamCS(nParamCS,1:6)
            end do

        else
            ! This parameter is not recognized; record an error.
            INFO = 1600 + i
            return
        end if

    end do LOOP_ExcessMixingRKMP

    ! Report an error if necessary:
    if (INFO /= 0) INFO = 1600 + i

    return

end subroutine ParseCSDataBlockRKMP
