
                cmake_minimum_required(VERSION 3.0)
                project(verify_mpi_compat VERSION 1.0 LANGUAGES C)
                file(WRITE compiler_info.cmake "
                    set(MPI_C_COMPILER_VERSION ${CMAKE_C_COMPILER_VERSION})
                    set(MPI_C_COMPILER_ID ${CMAKE_C_COMPILER_ID})
                    ")
            