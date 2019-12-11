SET(Fortran_FLAGS
    -Wall
    -g
    -O0
    -fno-automatic
    -fbounds-check
    -ffpe-trap=zero
    -D"DATA_DIRECTORY='${${PACKAGE_NAME}_SOURCE_DIR}/data/'"
)

FOREACH(flag ${Fortran_FLAGS})
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flag}")
ENDFOREACH()

UNSET(Fortran_FLAGS)
