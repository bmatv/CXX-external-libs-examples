//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>
#include <cmath>

#include <vector>

int main()
{
    /* This will be the netCDF ID for the file and data variable. */

    auto filepath = "/home/bogdan/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc";

    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("tomo", netCDF::NcGroup::Current);

    std::cout << "there are " << tomodata.getDimCount() << " dimensions in the tomo dataset\n";

    const unsigned int row_size = 2722;
    const unsigned int col_size = 2722;

    const std::vector<size_t> &start{0, 0, 0};           // 3*4, 12 bytes
    const std::vector<size_t> &count{1, row_size, col_size}; // 12 bytes

    std::vector<uint16_t> arr(row_size * col_size);

    // uint16_t data_out[8364144/2] {}; //check if 8192, 8388608 the actual pass is 8364032, but the limit is inconsistent
    // e.g. uint8_t data[8377999] {}; might pass and might throw a Segmentation fault

    tomodata.getVar(start, count, arr.data()); // the data is in, it's a 2D array

    for (int i = 500; i < 510; ++i)
    {
        for (int j = 500; j < 510; ++j)
        {
            std::cout << arr[i * row_size + j] << ' ';
        }
        std::cout << '\n';
    }

    long long accumulator = 0;
    for (int i = 0; i < row_size; ++i)
        for (int j = 0; j < col_size; ++j)
        {
            accumulator += arr[i * row_size + j];
        }
    std::cout << "The total sum is: " << accumulator << '\n';
    std::cout << "The average is: " << accumulator / (row_size * col_size) << '\n';

    const int avg_filter_size = 3;
    const int filter_offset = avg_filter_size / 2;
    // std::vector<uint16_t> smooth_arr((row_size - filter_offset * 2) * (col_size - filter_offset * 2),0);
    std::vector<uint16_t> smooth_arr(row_size * col_size,0);
    // uint16_t smooth_arr[row_size-filter_offset*2][col_size-filter_offset*2] {};

    // smoothing via averaging box
    for (int i = 0 + filter_offset; i < row_size - filter_offset; ++i)
        for (int j = 0 + filter_offset; j < col_size - filter_offset; ++j)
        {
            uint32_t tmp = 0;
            for (int x = -filter_offset; x <= filter_offset; ++x)
                for (int y = -filter_offset; y <= filter_offset; ++y)
                {
                    tmp += arr[(i + x) * row_size + j + y];
                }
            // std::cout << tmp / (avg_filter_size * avg_filter_size) << ' ';
            smooth_arr[i*row_size+j] = tmp / (avg_filter_size * avg_filter_size);
        }

    for (int i=500;i<510;++i){
    for (int j=500;j<510;++j){
        std::cout << smooth_arr[i*row_size + j] <<' ';
    }
    std::cout << '\n';}

    // sobel filter
    std::cout << "---Sobel---\n";
    std::vector<double> edges_arr(row_size * col_size,0);
    // overflows if int16

    for (int i=2;i<row_size-filter_offset*2;++i){
    for (int j=2;j<col_size-filter_offset*2;++j){
        auto gx = (smooth_arr[(i+1)*row_size + j-1] + 2*smooth_arr[i*row_size + j-1] + smooth_arr[(i-1)*row_size + j-1]
                +
                -smooth_arr[(i+1)*row_size + j+1] + -2*smooth_arr[i*row_size+j+1] + -smooth_arr[(i-1)*row_size+j+1]);

        auto gy = (smooth_arr[(i-1)*row_size+j+1] + 2*smooth_arr[(i-1)*row_size+j] + smooth_arr[(i-1)*row_size+j-1]
                +
                -smooth_arr[(i+1)*row_size+j+1] + -2*smooth_arr[(i+1)*row_size+j] + -smooth_arr[(i+1)*row_size+j-1]);

        edges_arr[(i)*row_size+j] = pow( (gx*gx + gy*gy),0.5);

    }
    }


    for (int i=500;i<510;++i){
    for (int j=500;j<510;++j){
        std::cout << edges_arr[i*row_size + j] <<' ';
    }
    std::cout << '\n';}


    //TODO filter the edges, potentially with a hard cutoff

    //TODO run hough transform for the edges left (non-zero) 
}