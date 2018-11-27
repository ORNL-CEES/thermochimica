

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

    integer :: i


    ! Loop through excess parameters:
    LOOP_ExcessMixingQKTO: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingQKTO

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
