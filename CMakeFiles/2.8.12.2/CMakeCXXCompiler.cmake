set(CMAKE_CXX_COMPILER "/calculate/spack/opt/spack/linux-ubuntu14.04-x86_64/gcc-4.8/openmpi-3.0.0-y2iqy5u7gtlzkuzwb44fslhx6owkxbia/bin/mpic++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_VERSION "4.8.4")
set(CMAKE_CXX_PLATFORM_ID "Linux")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCXX 1)
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()




set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpi_cxx;mpi;stdc++;m;pthread;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/calculate/spack/opt/spack/linux-ubuntu14.04-x86_64/gcc-4.8/hwloc-1.11.7-jcvgc3lqf3xincxtgzb227rpud766ngh/lib;/calculate/spack/opt/spack/linux-ubuntu14.04-x86_64/gcc-4.8/openmpi-3.0.0-y2iqy5u7gtlzkuzwb44fslhx6owkxbia/lib;/calculate/spack/opt/spack/linux-ubuntu14.04-x86_64/gcc-4.8/dealii-8.5.1-lbsqsz6havsivngvjqbk5s76qs453eev/lib;/usr/lib/gcc/x86_64-linux-gnu/4.8;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")



