SET(Fortran_FLAGS
    -fno-automatic
    -D"DATA_DIRECTORY='${${PACKAGE_NAME}_SOURCE_DIR}/data/'"
)
IF(${PROJECT_NAME}_GPROF})
    SET(Fortran_FLAGS
        ${Fortran_FLAGS}
        -pg
    )
ENDIF()

SET(Fortran_FLAGS_DEBUG
    -Wall
    -g
    -O0
    -fbounds-check
    -ffpe-trap=zero
)
SET(Fortran_FLAGS_RELEASE
    -O3
)

FOREACH(flag ${Fortran_FLAGS})
    STRING(REGEX MATCH "${flag}" ispresent "${CMAKE_Fortran_FLAGS}")
    IF(NOT ispresent)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flag}")
    ENDIF()
ENDFOREACH()
FOREACH(flag ${Fortran_FLAGS_DEBUG})
    STRING(REGEX MATCH "${flag}" ispresent "${CMAKE_Fortran_FLAGS_DEBUG}")
    IF(NOT ispresent)
        SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${flag}")
    ENDIF()
ENDFOREACH()
FOREACH(flag ${Fortran_FLAGS_RELEASE})
    STRING(REGEX MATCH "${flag}" ispresent "${CMAKE_Fortran_FLAGS_RELEASE}")
    IF(NOT ispresent)
        SET(CMAKE_Fortran_FLAGS_RELEASE "${CMEAK_Fortran_FLAGS_RELEASE} ${flag}")
    ENDIF()
ENDFOREACH()

UNSET(Fortran_FLAGS)
UNSET(Fortran_FLAGS_DEBUG)
UNSET(Fortran_FLAGS_RELEASE)
UNSET(flag)