//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>

int main(){
    /* This will be the netCDF ID for the file and data variable. */
    int ncid, varid,dimid;

//    int data_in[NX][NY];

    // Open the file and check to make sure it's valid.
    netCDF::NcFile dataFile("simple_xy.nc", netCDF::NcFile::read);

    // There are a number of inquiry functions in netCDF which can be
    // used to learn about an unknown netCDF file. In this case we know
    // that there are 2 netCDF dimensions, 4 netCDF variables, no
    // global attributes, and no unlimited dimension.

    std::cout<<"there are "<<dataFile.getVarCount()<<" variables"<<std::endl;
    std::cout<<"there are "<<dataFile.getAttCount()<<" attributes"<<std::endl;
    std::cout<<"there are "<<dataFile.getDimCount()<<" dimensions"<<std::endl;
    std::cout<<"there are "<<dataFile.getGroupCount()<<" groups"<<std::endl;
    std::cout<<"there are "<<dataFile.getTypeCount()<<" types"<<std::endl;

    auto xdim = dataFile.getDim("x");

    auto dims = dataFile.getDims();

    for(auto& entry:dims){
        std::cout <<"second.getName() = "<< entry.second.getName() <<" first = "<< entry.first <<" second.getSize() "<<entry.second.getSize() <<'\n';
    }

    std::cout<<"getDim returs "<<xdim.getName()<<std::endl;

}