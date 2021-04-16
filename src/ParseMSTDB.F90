subroutine ParseMSTDB

    USE ModuleThermoIO
    USE ModuleSS

    implicit none

    cThermoFileName = DATA_DIRECTORY // 'MSAX+CationVacancies.dat'

    call ParseCSDataFile(cThermoFileName)

    return

end subroutine ParseMSTDB
