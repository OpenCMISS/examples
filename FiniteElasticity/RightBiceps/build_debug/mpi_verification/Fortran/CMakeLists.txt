
                cmake_minimum_required(VERSION 3.0)
                project(verify_mpi_compat VERSION 1.0 LANGUAGES Fortran)
                file(WRITE compiler_info.cmake "
                    set(MPI_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
                    set(MPI_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID})
                    ")
            