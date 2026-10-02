# CMake helper for creating Polytope Python virtual environments.
#
# Usage:
#
#   polytope_add_python_environment(
#     build
#     "${CMAKE_BINARY_DIR}/venv"
#     "${CMAKE_BINARY_DIR}/lib"
#     DEPENDS polytope
#   )
#
# This creates a target named polytope_python_env_build. The environment gets a
# polytope.pth file containing a relative path from the venv site-packages
# directory to the requested pth target.

function(polytope_add_python_environment name env_root pth_target)
  set(options ALL)
  set(oneValueArgs TARGET_NAME)
  set(multiValueArgs DEPENDS)

  cmake_parse_arguments(arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT Python3_EXECUTABLE)
    message(FATAL_ERROR "polytope_add_python_environment requires Python3_EXECUTABLE")
  endif()

  if(NOT POLYTOPE_SITE_PACKAGES_PATH)
    message(FATAL_ERROR "polytope_add_python_environment requires POLYTOPE_SITE_PACKAGES_PATH")
  endif()

  if(arg_TARGET_NAME)
    set(target "${arg_TARGET_NAME}")
  else()
    set(target "polytope_python_env_${name}")
  endif()

  if(NOT IS_ABSOLUTE "${env_root}")
    get_filename_component(env_root "${env_root}" ABSOLUTE
      BASE_DIR "${CMAKE_BINARY_DIR}")
  endif()

  if(NOT IS_ABSOLUTE "${pth_target}")
    get_filename_component(pth_target "${pth_target}" ABSOLUTE
      BASE_DIR "${CMAKE_BINARY_DIR}")
  endif()

  set(env_site_packages "${env_root}/${POLYTOPE_SITE_PACKAGES_PATH}")
  set(stamp "${env_root}/polytope_venv.stamp")
  set(script "${CMAKE_CURRENT_BINARY_DIR}/${target}.cmake")

  file(RELATIVE_PATH pth_entry "${env_site_packages}" "${pth_target}")

  file(WRITE "${script}" "
execute_process(
  COMMAND \"${Python3_EXECUTABLE}\" -m venv --symlinks \"${env_root}\"
  RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR \"Failed to create Python virtual environment: ${env_root}\")
endif()

file(MAKE_DIRECTORY \"${env_site_packages}\")
file(WRITE \"${env_site_packages}/polytope.pth\" \"${pth_entry}\\n\")
file(WRITE \"${stamp}\" \"ok\\n\")
")

  add_custom_command(
    OUTPUT "${stamp}"
    COMMAND "${CMAKE_COMMAND}" -P "${script}"
    DEPENDS ${arg_DEPENDS}
    VERBATIM
  )

  if(arg_ALL)
    add_custom_target("${target}" ALL
      DEPENDS "${stamp}"
    )
  else()
    add_custom_target("${target}"
      DEPENDS "${stamp}"
    )
  endif()

  set(POLYTOPE_PYTHON_ENV_${name}_TARGET "${target}" PARENT_SCOPE)
  set(POLYTOPE_PYTHON_ENV_${name}_EXECUTABLE "${env_root}/bin/python" PARENT_SCOPE)
endfunction()
