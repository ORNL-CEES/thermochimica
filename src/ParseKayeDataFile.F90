subroutine ParseKayeDataFile

    USE ModuleThermoIO
    USE ModuleSS

    implicit none

    cThermoFileName = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    call ParseCSDataFile(cThermoFileName)

    return

end subroutine ParseKayeDataFile
