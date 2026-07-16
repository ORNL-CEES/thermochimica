FILE(GLOB THERMOCHIMICA_TESTS "*.F90")

FOREACH(test ${THERMOCHIMICA_TESTS})
    string( REPLACE ".F90" "" test ${test})
    get_filename_component(test ${test} NAME)
    add_executable(${test} ${test}.F90)
    set_target_properties(${test} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${test} thermochimica_lib ${lapackblas_libraries})
ENDFOREACH()

add_executable(TestPhaseConstraintAPI
    TestPhaseConstraintAPI.C
    ${CMAKE_SOURCE_DIR}/src/Thermochimica-c.C
    ${CMAKE_SOURCE_DIR}/src/Thermochimica-cxx.C)
target_include_directories(TestPhaseConstraintAPI PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_features(TestPhaseConstraintAPI PRIVATE cxx_std_17)
target_compile_definitions(TestPhaseConstraintAPI PRIVATE
    THERMOCHIMICA_DATA_DIRECTORY="${CMAKE_SOURCE_DIR}/data/")
target_link_libraries(TestPhaseConstraintAPI thermochimica_lib ${lapackblas_libraries})
