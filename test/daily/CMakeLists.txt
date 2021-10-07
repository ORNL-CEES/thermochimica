FILE(GLOB THERMOCHIMICA_TESTS "*.F90")

FOREACH(test ${THERMOCHIMICA_TESTS})
    string( REPLACE ".F90" "" test ${test})
    get_filename_component(test ${test} NAME)
    add_executable(${test} ${test}.F90)
    set_target_properties(${test} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${test} thermochimica_lib ${lapackblas_libraries})
ENDFOREACH()
