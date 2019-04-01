
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataFile.f90
    !> \brief   Parse a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      ParseCSHeader.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ModuleParseCS.f90
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
    !   10/05/2011      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse a ChemSage data-file into a format that can be
    !! used by Thermochimica.  The intended purpose of this program is to be called once in AMP
    !! and to store the parsed information to memory to minimize the number of times (to one)
    !! that a data-file has to be read from the hard disk.
    !!
    !! First, the "header section" of the data-file is parsed to gather important information to
    !! allocate memory.  Second, the "data block" section is parsed to store thermodynamic data
    !! of all species and phases in the system.  A description of all the pertinent variables
    !! can be found in the ParseCSModule module.
    !!
    !! This subroutine is looking for a specific file specified by the variable cFileName and
    !! it will return a non-zero value for INFO if there is an error.  The error code returned
    !! by this program is represented by the integer scalar INFO and is described below:
    !!
    !! <table border="1" width="600" align="center|left">
    !! <tr>
    !!    <td> <b> INFO Value </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> 0 </td> <td>  Successful exit </td>
    !! </tr>
    !! <tr>
    !!    <td> 1 </td> <td> The specified data-file was not found in the local directory </td>
    !! </tr>
    !! <tr>
    !!    <td> 2 </td> <td> The number of elements in the system exceeds the maximum allowable
    !!                   number </td>
    !! </tr>
    !! <tr>
    !!    <td> 3 </td> <td> The number of solution phases in the system exceeds the maximum
    !!                   allowable number </td>
    !! </tr>
    !! <tr>
    !!    <td> 4 </td> <td> The solution phase type is not supported </td>
    !! </tr>
    !! <tr>
    !!    <td> 5 </td> <td> The number of Gibbs energy equations has exceeded the amount allocated
    !!                   to memory </td>
    !! </tr>
    !! <tr>
    !!    <td> 1x </td> <td> Error reading line x of the header section (e.g., INFO = 13 means
    !!                   that line 3 of the header section failed). The line numbering convention
    !!                   from the ChemApp programmer's library is used here </td>
    !! </tr>
    !! <tr>
    !!    <td> 1xab </td> <td> Error reading entry x of the data block section, where ab refers to the
    !!                   number of the solution phase that is causing the problem (e.g.,
    !!                   INFO = 1234 means that entry 2 of the data block section had an error
    !!                   with solution phase 34). The entry numbering convention from the ChemApp
    !!                   programmer's library is used here. </td>
    !! </tr>
    !! </table>
    !
    !
    ! Limitations:
    ! ============
    !
    !> \details This program is not capable of parsing all types of ChemSage data-files.  The following is
    !! a summary of known limitations:
    !!   - Only IDMX and QKTO solution phases can be considered.
    !!   - Pure separate phases can be considered.
    !!   - Aqueous phases cannot be considered.
    !!   - Only Gibbs energy equations of the '4' type can be considered.
    !!   - The maximum number of solution phases is 42.  The only reason why this is constrained is
    !!     because this is the maximum number of solution phases that FactSage can handle.
    !!   - The maximum number of elements is 118 (limited by the periodic table)
    !!   - The maximum number of components in a sub-system for a particular interaction parameter
    !!     is set to 4 (i.e., quaternary).  Note that this limitation is for a sub-system, not for
    !!     a solution phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]  cFileName    A character string representing the path and name of the ChemSage data-file.
    !  INFO         An integer scalar representing a successful exit or an error.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataFile(cFileName)

    USE ModuleParseCS
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    character(*)::  cFileName


    ! Initialize variables:
    INFOThermo = 0
    INFO       = 0

    ! Attempt to open the ChemSage datafile:
    open (UNIT = 1, FILE = cFileName, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)

    ! Record an error if there are issues opening the data-file:
    if (INFO /= 0) then
        INFOThermo = 6
        return
    end if

    ! Parse the "header section" of the data-file:
    if (INFOThermo == 0) call ParseCSHeader

    if (INFO /= 0) INFOThermo = INFO

    ! Parse the "data block section" of the data-file:
    if (INFOThermo == 0) call ParseCSDataBlock

    if (INFO /= 0) INFOThermo = INFO

    ! Attempt to close the data-file:
    close (UNIT = 1, STATUS = 'keep', IOSTAT = INFO)

    ! Make sure that there aren't any issues closing the data-file:
    if (INFO /= 0) INFOThermo = 6

    return

end subroutine ParseCSDataFile

