subroutine ParseMSTDB

    USE ModuleThermoIO
    USE ModuleSS

    implicit none

    cThermoFileName = DATA_DIRECTORY // 'MSAX+F2.dat'

    call ParseCSDataFile(cThermoFileName)

    return

end subroutine ParseMSTDB
