if(ENABLE_LCAO)
  add_library(
	numerical_atomic_orbitals
	OBJECT
	numerical_radial.cpp
    radial_set.cpp
    atomic_radials.cpp
    beta_radials.cpp
    radial_collection.cpp
    real_gaunt_table.cpp
  )

  if(ENABLE_COVERAGE)
	add_coverage(numerical_atomic_orbitals)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
