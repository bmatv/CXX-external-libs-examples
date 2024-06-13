//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>

int main(){
    /* This will be the netCDF ID for the file and data variable. */

    // auto filepath = "/home/bogdan/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc";
    auto filepath = "simple_xy.nc";


    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("data");

    int data_out[6][12];
    // short tomoIn[10];
    tomodata.getVar(data_out);

    for (int x = 0; x < 6; x++)
    {
        for (int y = 0; y < 12; y++)
        {
            std::cout<< data_out[x][y]<<' ';
            }

        std::cout <<'\n';
        }


}