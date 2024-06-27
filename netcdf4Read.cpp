//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <netcdf>
#include <cmath>
#include <cstdint>

#include <vector>

int main()
{
    /* This will be the netCDF ID for the file and data variable. */

    auto filepath = "/home/bogdanm/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc";

    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("tomo", netCDF::NcGroup::Current);

    std::cout << "there are " << tomodata.getDimCount() << " dimensions in the tomo dataset\n";

    const int row_size = 2722;
    const int col_size = 2722;

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

    const int avg_filter_size = 5;
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
    double maxSobel = 0;
    uint64_t countNonZero = 0;
    double nonZeroSum = 0;
    double nonZeroSquaredSum = 0;

    for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i){
    for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        auto gx = (smooth_arr[(i+1)*row_size + j-1] + 2*smooth_arr[i*row_size + j-1] + smooth_arr[(i-1)*row_size + j-1]
                +
                -smooth_arr[(i+1)*row_size + j+1] + -2*smooth_arr[i*row_size+j+1] + -smooth_arr[(i-1)*row_size+j+1]);

        auto gy = (smooth_arr[(i-1)*row_size+j+1] + 2*smooth_arr[(i-1)*row_size+j] + smooth_arr[(i-1)*row_size+j-1]
                +
                -smooth_arr[(i+1)*row_size+j+1] + -2*smooth_arr[(i+1)*row_size+j] + -smooth_arr[(i+1)*row_size+j-1]);

        auto sobelVal = pow( (gx*gx + gy*gy),0.5);

        // edges_arr[(i)*row_size+j] = sobelVal<15000?sobelVal:0; // then we could ignore the mask but need to check whether the value is Nan in case of float tomo

        if (sobelVal != 0){
            if (sobelVal>maxSobel){
                maxSobel = sobelVal;
            }
            countNonZero++;
            nonZeroSum+=sobelVal;
            nonZeroSquaredSum+=sobelVal*sobelVal;

            edges_arr[i*row_size+j] = sobelVal;

        }


        
    }
    }

    std::cout << "Maximum sobel value = " << maxSobel << '\n';
    std::cout << "NonZero values count: " << countNonZero << '\n';
    std::cout << "Zero values count: " << (row_size * col_size) - countNonZero << '\n';
    std::cout << "Average value: " << nonZeroSum/countNonZero << '\n';

    auto meanSobelVal = nonZeroSum/countNonZero;
    auto meanSumSquared = nonZeroSquaredSum/countNonZero;

    auto stdSobel = pow(meanSumSquared - meanSobelVal*meanSobelVal,0.5); // approx std. good enough

    std::cout << "STD value: " << stdSobel << '\n' << "Filter bounds:" << meanSobelVal-stdSobel << " to " << meanSobelVal+stdSobel << '\n';


    //TODO filter the edges, potentially with a hard cutoff
    long long filteredSobelCount = 0;
    long long iAvg {};
    long long jAvg {};
    long long iStd {};
    long long jStd {};

    // compute average edge i and j, and respective standard deviations. Those will be our center coordinate and search box dimensions.


    for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i)
    for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        if (edges_arr[i*row_size+j]>(meanSobelVal+stdSobel) || edges_arr[i*row_size+j] < meanSobelVal){
                edges_arr[i*row_size+j] = 0;
        }
        else { // edge found
            filteredSobelCount++;
            iAvg+=i;
            jAvg+=j;
            iStd+=i*i;
            jStd+=j*j;
            }
    }

    iAvg = iAvg/filteredSobelCount;
    jAvg = jAvg/filteredSobelCount;
    iStd = pow(iStd/filteredSobelCount - iAvg,0.5);
    jStd = pow(jStd/filteredSobelCount - jAvg,0.5);
    std::cout << "Zero values count after filtering: " <<filteredSobelCount <<'\n';
    std::cout << "Average i = " << iAvg <<"\nSTD i = "<< iStd <<'\n';
    std::cout << "Average j = " << jAvg <<"\nSTD j = "<< jStd <<'\n';

    // The question is: which metric to use for hough transform - points being +=1 if non-zero while could also add up intensities 
    // (high intensity targets could produce anomalous signals that will overwhelm the transform)
    // long long test = iAvg - iStd*0.1;



    // const double searchBoxBoundsSTDCoef = 1.0/100.0;
    // long long iSearchSize = iStd*searchBoxBoundsSTDCoef;
    // long long jSearchSize = jStd*searchBoxBoundsSTDCoef;
    long long iSearchSize = 64; // has to be even otherwise not all values will be accesible
    long long jSearchSize = 64;
    int halfDiagLen = pow(pow(iSearchSize,2) + pow(jSearchSize,2),0.5);

    std::vector<uint64_t>searchBox(iSearchSize*jSearchSize,0);

    std::cout << "centre search box coordinates (i): " << iAvg - iSearchSize/2 << " to " << iAvg + iSearchSize/2 << "\niSearchSize (Window) = " << iSearchSize <<'\n';
    std::cout << "centre search box coordinates (j): " << jAvg - jSearchSize/2 << " to " << jAvg + jSearchSize/2 << "\njSearchSize (Window) = " << jSearchSize <<'\n';

    // int radius_min = 1200, radius_max = 1205;
    int radius = 1192;



    for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i) // for every pixel of the tomogram
    for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        if (edges_arr[i*row_size+j]!=0){ // real edge, compute distance from every centre point

            //check if the edge is not within 0.5 (0.7?) std from the avg centre


            int distanceToAvgCentre = pow(pow(i - iAvg,2) + pow(j - jAvg,2),0.5);
            // for radius in radii should be here

            if (distanceToAvgCentre < 0.5*(iStd)) // needs to be distance
                continue;

            if ((distanceToAvgCentre-halfDiagLen) <= radius && (distanceToAvgCentre+halfDiagLen) >= radius){ // the edge could theoretically be part of the circle with the centre in the search box so we can continue
                for(auto iCentre = iAvg - iSearchSize/2; iCentre < iAvg +iSearchSize/2;iCentre++) 
                for(auto jCentre = jAvg - jSearchSize/2; jCentre < jAvg + jSearchSize/2;jCentre++){

                    auto sideA = pow(i - iCentre,2);
                    auto sideB = pow(j - jCentre,2);

                    int distance = (int) pow( sideA + sideB,0.5);
                    long long indexI = 0;
                    long long indexJ = 0;
                    if(distance == radius){
                        indexI = iCentre-iAvg+ iSearchSize/2;
                        indexJ = jCentre-jAvg + jSearchSize/2;
                        searchBox.at(indexI * jSearchSize + indexJ ) += 1; // searchBox[0]  // may be too much write accesses, could be better to revise the loop so a tmp var could be formed
                    }


                }

            }

            //     }

        }

    }

    // for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i)
    // for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
    
    
    for(auto i = 0; i < iSearchSize;i++){
    for(auto j = 0; j < jSearchSize;j++){
        std::cout << searchBox.at(i*jSearchSize + j) <<' ';
    }
    std::cout << '\n';
    }


    // write out the resultant netcdf after sobel filtering?




    //TODO run hough transform for the edges left (non-zero) 
}