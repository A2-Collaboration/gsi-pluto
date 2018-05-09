# Install script for directory: /home/wagners/git/pluto_v6.00/plugins

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/wagners/git/pluto_v6.00/plugins/eta_decays/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/w_decays/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/strangeness/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/pion_beam/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/beamline/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/pdg_unigen/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/tools/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/scatter_mod/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/fairroot/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/elementary/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/dalitz_mod/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/brems/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/hades/cmake_install.cmake")
  include("/home/wagners/git/pluto_v6.00/plugins/nucleus_fermi/cmake_install.cmake")

endif()

