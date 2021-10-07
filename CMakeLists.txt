PROJECT(Thermochimica)
enable_language (Fortran)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90 F90)
CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

SET(Fortran_FLAGS
    -fno-automatic
    -D"DATA_DIRECTORY='${${PACKAGE_NAME}_SOURCE_DIR}/data/'"
)

FOREACH(flag ${Fortran_FLAGS})
    STRING(REGEX MATCH "${flag}" ispresent "${CMAKE_Fortran_FLAGS}")
    IF(NOT ispresent)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flag}")
    ENDIF()
ENDFOREACH()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
