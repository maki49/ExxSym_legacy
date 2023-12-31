list(APPEND objects
    diago_cg.cpp
    diago_david.cpp
    hsolver_pw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      hsolver_lcao.cpp
      diago_blas.cpp
  )
  if (USE_ELPA)
    list(APPEND objects
        diago_elpa.cpp
    )
  endif ()

  if(USE_CUSOLVER_LCAO)
      list(APPEND objects
          diag_cusolver.cu
      )

    add_library(
        diag_cusolver
        OBJECT
        ${objects}
    )
    if(ENABLE_COVERAGE)
      add_coverage(diag_cusolver)
    endif()
  endif()
endif()

add_library(
    hsolver
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hsolver)
endif()

if(ENABLE_LCAO AND USE_ELPA)
  add_subdirectory(genelpa)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
  add_subdirectory(kernels/test)
endif()
