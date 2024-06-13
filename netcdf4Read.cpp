//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>


int main(){
    /* This will be the netCDF ID for the file and data variable. */

    auto filepath = "/home/bogdan/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc";
    // auto filepath = "simple_xy.nc";


    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("tomo",netCDF::NcGroup::Current);

    std::cout << "there are "<< tomodata.getDimCount() <<" dimensions in the dataset\n";
    
    // std::vector<int> data_out(6*12); // segfault without the size
    std::vector<u_int16_t> data_out(1*2722*2722); // segfault without the size
    // int data_out[1*2722*2722] {};
    // int data_out[6*12] {};
    tomodata.getVar(data_out.data());

    for (int x = 2772*1300; x < 2772*1301; x++)
    {
        std::cout<< data_out[x]<<' ';
    }

}