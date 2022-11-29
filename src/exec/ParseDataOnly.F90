program ParseDataOnly

    USE ModuleThermoIO
    USE ModuleParseCS

    implicit none

    integer :: i, j
    logical :: lMisc
    character(1024) :: cThermoFileNameTemp

    if (COMMAND_ARGUMENT_COUNT() < 1) then
        print *, 'One command-line argument (filename) required.'
        call EXIT(1)
    endif

    call GET_COMMAND_ARGUMENT(1,cThermoFileNameTemp)

    cThermoFileName = cThermoFileNameTemp

    call ParseCSDataFile(cThermoFileName)

    ! Write JSON with phase lists
    open(1, file = DATA_DIRECTORY // &
            '../outputs/phaseLists.json', status='REPLACE', action='write')

    write(1,*) '{'
    write(1,*) '  "solution phases": {'
    do i = 1, nSolnPhasesSysCS
        ! Check if a miscibility gap phase
        lMisc = .FALSE.
        if (i > 1) then
            if (cSolnPhaseNameCS(i) == cSolnPhaseNameCS(i-1)) lMisc = .TRUE.
        end if

        if (lMisc) then
            write(1,'(A,A,A,I0,A)') '    "', TRIM(ADJUSTL(cSolnPhaseNameCS(i))), '#', i, '": {'
        else
            write(1,*) '    "', TRIM(ADJUSTL(cSolnPhaseNameCS(i))), '": {'
        end if
        write(1,*) '      "species": ['
        do j = nSpeciesPhaseCS(i - 1) + 1, nSpeciesPhaseCS(i)
            if (j < nSpeciesPhaseCS(i)) then
                write(1,*) '        "', TRIM(ADJUSTL(cSpeciesNameCS(j))), '",'
            else
                write(1,*) '        "', TRIM(ADJUSTL(cSpeciesNameCS(j))), '"'
            end if
        end do

        ! Write endmember names for MQMQA phases
        if ((cSolnPhaseTypeCS(i) == 'SUBG') .OR. (cSolnPhaseTypeCS(i) == 'SUBQ')) then
            write(1,*) '      ],'
            write(1,*) '      "endmembers": ['
            
            do j = 1, nPairsSROCS(iPhaseSublatticeCS(i),1)
                if (j < nPairsSROCS(iPhaseSublatticeCS(i),1)) then
                    write(1,*) '        "', TRIM(ADJUSTL(cPairNameCS(iPhaseSublatticeCS(i),j))), '",'
                else
                    write(1,*) '        "', TRIM(ADJUSTL(cPairNameCS(iPhaseSublatticeCS(i),j))), '"'
                end if
            end do
        end if

        write(1,*) '      ]'

        if (i < nSolnPhasesSysCS) then
            write(1,*) '    },'
        else
            write(1,*) '    }'
        end if
    end do
    write(1,*) '  },'

    write(1,*) '  "pure condensed phases": ['
    do i = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
        if (i < nSpeciesCS) then
            write(1,*) '    "', TRIM(ADJUSTL(cSpeciesNameCS(i))), '",'
        else
            write(1,*) '    "', TRIM(ADJUSTL(cSpeciesNameCS(i))), '"'
        end if
    end do
    write(1,*) '  ]'
    write(1,*) '}'

    close (1)

end program
