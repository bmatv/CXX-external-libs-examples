//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>
#include <cmath>


int main(){
    /* This will be the netCDF ID for the file and data variable. */

    auto filepath = "/home/bogdan/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc";
    // auto filepath = "simple_xy.nc";


    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("tomo",netCDF::NcGroup::Current);

    std::cout << "there are "<< tomodata.getDimCount() <<" dimensions in the dataset\n";
    

    unsigned int row_size = 15;
    unsigned int col_size = 15;
    
    const std::vector<size_t>& start {0,500,500}; // 3*4, 12 bytes 
    const std::vector<size_t>& count {1,row_size,col_size}; // 12 bytes

    u_int16_t arr [row_size][col_size];
    // std::vector<u_int16_t> data_out(row_size*col_size);

    // uint16_t data_out[8364144/2] {}; //check if 8192, 8388608 the actual pass is 8364032, but the limit is inconsistent
    // e.g. uint8_t data[8377999] {}; might pass and might throw a Segmentation fault

    // uint8_t data[8377999] {};

    // std::vector<u_int16_t> data_out(2722*2722);
    // int data_out[3] {};
    // tomodata.getVar(start,count,data_out.data()); // the data is in, it's a 1D array
    tomodata.getVar(start,count, arr); // the data is in, it's a 2D array


    for (int i =0;i<row_size;++i){
    for (int j=0;j<col_size;++j){
        std::cout << arr[i][j] <<' ';
    }
    std::cout << '\n';
    }

    long long accumulator = 0;
    for (int i=0;i<row_size;++i)
    for(int j=0;j<col_size;++j)
    {
    accumulator+=arr[i][j];
    }
    std::cout << "The total sum is: "<< accumulator <<'\n';
    std::cout << "The average is: "<< accumulator/(row_size*col_size) <<'\n';

    int avg_filter_size = 3;
    int filter_offset = avg_filter_size/2;
    uint16_t smooth_arr[row_size-filter_offset*2][col_size-filter_offset*2] {};


    // smoothing via averaging box
    for (int i=0 + filter_offset;i<row_size-filter_offset;++i)
    for(int j=0 + filter_offset;j<col_size-filter_offset;++j)
    {
        uint32_t tmp = 0;
        for(int x=-filter_offset;x<=filter_offset;++x)
        for(int y=-filter_offset;y<=filter_offset;++y){
            tmp += arr[i+x][j+y];
        }

        smooth_arr[i-filter_offset][j-filter_offset] = tmp / (avg_filter_size*avg_filter_size);
    }

    for (int i =0;i<row_size-filter_offset*2;++i){
    for (int j=0;j<col_size-filter_offset*2;++j){
        std::cout << smooth_arr[i][j] <<' ';
    }
    std::cout << '\n';
    }


    // sobel filter
    std::cout << "---Sobel---\n";
    double edges_arr[row_size-filter_offset*2-2][col_size-filter_offset*2-2] {};
    // overflows if int16


    for (int i=1;i<row_size-filter_offset*2-1;++i){
    for (int j=1;j<col_size-filter_offset*2-1;++j){
        auto gx = (smooth_arr[i+1][j-1] + 2*smooth_arr[i][j-1] + smooth_arr[i-1][j-1]
                +
                -smooth_arr[i+1][j+1] + -2*smooth_arr[i][j+1] + -smooth_arr[i-1][j+1]);
                
        auto gy = (smooth_arr[i-1][j+1] + 2*smooth_arr[i-1][j] + smooth_arr[i-1][j-1]
                +
                -smooth_arr[i+1][j+1] + -2*smooth_arr[i+1][j] + -smooth_arr[i+1][j-1]);

        edges_arr[i-1][j-1] = pow( (gx*gx + gy*gy),0.5);

        std::cout << edges_arr[i-1][j-1] <<' ';

        
    }
    std::cout << '\n';
    }
}