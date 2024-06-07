cmake_minimum_required(VERSION 3.25)
project(netcdf4Examples)

set(CMAKE_CXX_STANDARD 17)

find_package(netCDF REQUIRED)
find_package(MPI REQUIRED)

add_executable(netcdf4Write netcdf4Write.cpp)
target_link_libraries(netcdf4Write PRIVATE netCDF::netcdf) # should be netcdf here but need to look or netCDF

add_executable(mpiExample mpiExample.cpp)
target_link_libraries(mpiExample PRIVATE MPI::MPI_CXX)