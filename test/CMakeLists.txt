set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90 F90)
SET(UNIT_TEST_NAMES
      daily
)

FOREACH(test ${UNIT_TEST_NAMES})
      ADD_SUBDIRECTORY(${test})
ENDFOREACH()

UNSET(UNIT_TEST_NAMES)
UNSET(test)

add_executable(Thermo Thermo.F90)
set_target_properties(Thermo PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(Thermo thermochimica_lib ${lapackblas_libraries})

add_executable(ThermochimicaInputScriptMode ThermochimicaInputScriptMode.F90)
set_target_properties(ThermochimicaInputScriptMode PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(ThermochimicaInputScriptMode thermochimica_lib ${lapackblas_libraries})
