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

    // smoothing via averaging box, might want to make smooth float, so sobel could be faster?
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
            smooth_arr[i*row_size+j] = tmp / (avg_filter_size * avg_filter_size); // int math here might be slow
        }

    for (int i=500;i<510;++i){
    for (int j=500;j<510;++j){
        std::cout << smooth_arr[i*row_size + j] <<' ';
    }
    std::cout << '\n';}

    // sobel filter
    std::cout << "---Sobel---\n";
    std::vector<float> edges_arr(row_size * col_size,0);
    // overflows if int16
    double maxSobel = 0;
    uint64_t countNonZero = 0;
    double nonZeroSum = 0;
    double nonZeroSquaredSum = 0;

    for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i){
    for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        float gx = static_cast<float>(
            smooth_arr[(i + 1) * row_size + j - 1] + 2 * smooth_arr[i * row_size + j - 1] +
            smooth_arr[(i - 1) * row_size + j - 1] + -smooth_arr[(i + 1) * row_size + j + 1] +
            -2 * smooth_arr[i * row_size + j + 1] + -smooth_arr[(i - 1) * row_size + j + 1]);

        float gy = static_cast<float>(
            smooth_arr[(i - 1) * row_size + j + 1] + 2 * smooth_arr[(i - 1) * row_size + j] +
            smooth_arr[(i - 1) * row_size + j - 1] + -smooth_arr[(i + 1) * row_size + j + 1] +
            -2 * smooth_arr[(i + 1) * row_size + j] + -smooth_arr[(i + 1) * row_size + j - 1]);

        auto sobelVal = std::sqrt(gx*gx + gy*gy);

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

    auto meanSobelVal = static_cast<float>(nonZeroSum/countNonZero);
    auto meanSumSquared = static_cast<float>(nonZeroSquaredSum/countNonZero);

    auto stdSobel = std::sqrt(static_cast<float>(meanSumSquared - meanSobelVal*meanSobelVal)); // approx std. good enough

    std::cout << "STD value: " << stdSobel << '\n' << "Filter bounds:" << meanSobelVal-stdSobel << " to " << meanSobelVal+stdSobel << '\n';


    //TODO filter the edges, potentially with a hard cutoff
    long long filteredSobelCount = 0;
    long long iAvg {};
    long long jAvg {};
    long long iStdAcc {};
    long long jStdAcc {};

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
            iStdAcc+=i*i;
            jStdAcc+=j*j;
            }
    }

    iAvg = iAvg/filteredSobelCount;
    jAvg = jAvg/filteredSobelCount;
    float iStd = std::sqrt(static_cast<float>(iStdAcc/filteredSobelCount - iAvg)); // could be a new float var?
    float jStd = std::sqrt(static_cast<float>(jStdAcc/filteredSobelCount - jAvg));
    std::cout << "Zero values count after filtering: " <<filteredSobelCount <<'\n';
    std::cout << "Average i = " << iAvg <<"\nSTD i = "<< iStd <<'\n';
    std::cout << "Average j = " << jAvg <<"\nSTD j = "<< jStd <<'\n';

    // The question is: which metric to use for hough transform - points being +=1 if non-zero while could also add up intensities 
    // (high intensity targets could produce anomalous signals that will overwhelm the transform)
    // long long test = iAvg - iStd*0.1;



    // const double searchBoxBoundsSTDCoef = 1.0/100.0;
    // long long iSearchSize = iStd*searchBoxBoundsSTDCoef;
    // long long jSearchSize = jStd*searchBoxBoundsSTDCoef;
    long long iSearchSize = 256; // has to be even otherwise not all values will be accesible
    long long jSearchSize = 256;
    float halfDiagLen = std::sqrt(static_cast<float>(iSearchSize*iSearchSize + jSearchSize*jSearchSize));
    std::vector<int>radii {1205,1206,1207,1208,1209}; //1195,1200,1205,1210,1215,1220
    std::vector<uint64_t>searchBox(iSearchSize*jSearchSize*radii.size(),0);

    std::cout << "centre search box coordinates (i): " << iAvg - iSearchSize/2 << " to " << iAvg + iSearchSize/2 << "\niSearchSize (Window) = " << iSearchSize <<'\n';
    std::cout << "centre search box coordinates (j): " << jAvg - jSearchSize/2 << " to " << jAvg + jSearchSize/2 << "\njSearchSize (Window) = " << jSearchSize <<'\n';

    // int radius_min = 1200, radius_max = 1205;
    // int radius = 1192;

    // int radii[] = {1192,1193,1194};
    
    



    for (int i=filter_offset+1;i<row_size-(filter_offset+1);++i) // for every pixel of the tomogram
    for (int j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        if (edges_arr[i*row_size+j]!=0){ // real edge, compute distance from every centre point

            // check if the edge is not within 0.5 (0.7?) std from the avg centre

            float distanceToAvgCentre = std::sqrt(std::pow(static_cast<float>(i - iAvg), 2.0F) + std::pow(static_cast<float>(j - jAvg), 2.0F));
            // for radius in radii should be here

            if (distanceToAvgCentre < iStd/2.0F) // iStd is ~= radius hence 0.5iStd is a half of that which is where we do not expect the edges of the cyllinder to be
                continue;
            for (size_t rIdx =0; rIdx<radii.size();rIdx++){
                if (static_cast<int>(distanceToAvgCentre-halfDiagLen) <= radii[rIdx] && static_cast<int>(distanceToAvgCentre+halfDiagLen) >= radii[rIdx]){ // the edge could theoretically be part of the circle with the centre in the search box so we can continue
                    for(auto iCentre = iAvg - iSearchSize/2; iCentre < iAvg +iSearchSize/2;++iCentre) 
                    for(auto jCentre = jAvg - jSearchSize/2; jCentre < jAvg + jSearchSize/2;++jCentre){

                        int distance = static_cast<int>(std::sqrt( std::pow(static_cast<float>(i - iCentre),2) + std::pow(static_cast<float>(j - jCentre),2)));
                    long long indexI = 0;
                    long long indexJ = 0;
                        if(distance == radii[rIdx]){
                            indexI = iCentre-iAvg + iSearchSize/2;
                        indexJ = jCentre-jAvg + jSearchSize/2;
                            searchBox.at(rIdx*(iSearchSize*jSearchSize) + indexI * jSearchSize + indexJ ) += 1; // searchBox[0]  // may be too much write accesses, could be better to revise the loop so a tmp var could be formed
                    }


                }

            }

            }

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