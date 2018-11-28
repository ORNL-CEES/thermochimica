
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseInput.f90
    !> \brief   Parse the contents of an input script file.
    !> \author  M. Poschmann
    !> \date    Nov. 28, 2018
    !> \sa      CheckThermoInput.f90
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   28/11/2018      M. Poschmann         Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to parse an input file for Thermochimica
    !
    !
    ! Pertinent variables:
    ! ====================
    !> \param[in]   cInputFileName     The name of the input file to be read.
    !
    ! dTemperature              Temperature (converted to K)
    ! dPresssure                Absolute hydrostatic pressure (converted to atm)
    ! dElementMass              Total mass of an element, where the coefficient corresponds to the
    !                            atomic number (e.g., dMolesElement(92) refers to uranium).
    ! cThermoInputUnits         A character vector containing the units for temperature, pressure and
    !                            mass.
    ! INFOThermo                A scalar integer that indicates a successful exit or identifies an error.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine ParseInput(cInputFileName)

  USE ModuleThermoIO

  implicit none

  integer               :: iCounter
  logical               :: lEnd
  character(40)         :: cLine

  INFO = 0
  open (UNIT = 1, FILE = cInputFileName, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)
  if (INFO /= 0) then
      INFOThermo = 7
      return
  end if

  lEnd = .FALSE.
  do while (.NOT. lEnd)
    read (1,*,IOSTAT = INFO) cLine

  enddo









end subroutine ParseInput
