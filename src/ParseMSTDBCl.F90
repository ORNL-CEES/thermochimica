subroutine ParseMSTDBCl

    USE ModuleThermoIO
    USE ModuleSS

    implicit none

    cThermoFileName = DATA_DIRECTORY // 'CLAK+Cl2.dat'

    call ParseCSDataFile(cThermoFileName)

    return

end subroutine ParseMSTDBCl
