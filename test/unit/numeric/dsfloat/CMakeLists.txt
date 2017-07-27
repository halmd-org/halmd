if(HALMD_WITH_GPU)
    cuda_add_executable(test_unit_dsfloat_dsfloat
            dsfloat.cpp
            dsfloat.cu
            )
    target_link_libraries(test_unit_dsfloat_dsfloat
            ${HALMD_TEST_LIBRARIES}
            )
    halmd_add_gpu_test(unit/numeric/dsfloat
            test_unit_dsfloat_dsfloat --log_level=test_suite
            )
endif()
