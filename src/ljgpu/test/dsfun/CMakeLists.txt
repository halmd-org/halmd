foreach(test
    add
    div
    mul
    sqrt
    mulss
    sub
    )
  add_executable(test_dsfun_${test}
    ${test}.cpp
    ${test}_kernel.cu
    )
  target_link_libraries(test_dsfun_${test}
    ${CUDA_LIBRARIES}
    )
  add_test(dsfun_${test} test_dsfun_${test})
endforeach(test)
