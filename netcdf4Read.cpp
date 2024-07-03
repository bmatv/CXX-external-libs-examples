//
// Created by bogdan on 11/06/24.
//
#include <iostream>
#include <fstream>

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

    const int32_t row_size = 2722;
    const int32_t col_size = 2722;

    const std::vector<size_t> &start{0, 0, 0};           // 3*4, 12 bytes
    const std::vector<size_t> &count{1, row_size, col_size}; // 12 bytes

    std::vector<uint16_t> arr(row_size * col_size);

    // uint16_t data_out[8364144/2] {}; //check if 8192, 8388608 the actual pass is 8364032, but the limit is inconsistent
    // e.g. uint8_t data[8377999] {}; might pass and might throw a Segmentation fault

    tomodata.getVar(start, count, arr.data()); // the data is in, it's a 2D array
    //writing mask value as it is missing in the file
    for(auto& val:arr){
        val = val!=0?val:static_cast<uint16_t>(65535); 
    }

    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            std::cout << arr[i * row_size + j] << ' ';
        }
        std::cout << '\n';
    }

    // long long accumulator = 0;
    // for (int i = 0; i < row_size; ++i)
    //     for (int j = 0; j < col_size; ++j)
    //     {
    //         accumulator += arr[i * row_size + j];
    //     }
    // std::cout << "The total sum is: " << accumulator << '\n';
    // std::cout << "The average is: " << accumulator / (row_size * col_size) << '\n';

    const int32_t avg_filter_size = 5;
    const int32_t filter_offset = avg_filter_size / 2;
    // std::vector<uint16_t> smooth_arr((row_size - filter_offset * 2) * (col_size - filter_offset * 2),0);
    std::vector<int32_t> smooth_arr(row_size * col_size,0);
    // uint16_t smooth_arr[row_size-filter_offset*2][col_size-filter_offset*2] {};
    std::cout << "---Average---\n";
    // smoothing via averaging box, might want to make smooth float, so sobel could be faster?
    for (auto i = 0 + filter_offset; i < row_size - filter_offset; ++i)
        for (auto j = 0 + filter_offset; j < col_size - filter_offset; ++j)
        {
            int32_t tmp = 0;
            for (auto x = i-filter_offset; x <= i+filter_offset; ++x){
                for (auto y = j-filter_offset; y <= j+filter_offset; ++y){
                    tmp += arr[x * row_size + y];
                }}
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

    for (auto i=filter_offset+1;i<row_size-(filter_offset+1);++i){
    for (auto j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        auto gx = static_cast<float>(
            smooth_arr[(i + 1) * row_size + j - 1] + 2 * smooth_arr[i * row_size + j - 1] +
            smooth_arr[(i - 1) * row_size + j - 1] + -smooth_arr[(i + 1) * row_size + j + 1] +
            -2 * smooth_arr[i * row_size + j + 1] + -smooth_arr[(i - 1) * row_size + j + 1]);

        auto gy = static_cast<float>(
            smooth_arr[(i - 1) * row_size + j + 1] + 2 * smooth_arr[(i - 1) * row_size + j] +
            smooth_arr[(i - 1) * row_size + j - 1] + -smooth_arr[(i + 1) * row_size + j + 1] +
            -2 * smooth_arr[(i + 1) * row_size + j] + -smooth_arr[(i + 1) * row_size + j - 1]);

        auto sobelVal = std::sqrt(gx*gx + gy*gy);

        sobelVal = sobelVal<50000.0F?sobelVal:0; // then we could ignore the mask but need to check whether the value is Nan in case of float tomo

        if (sobelVal != 0){
            if (sobelVal>maxSobel){ // just for testing?
                maxSobel = sobelVal;
            }
            countNonZero++;
            nonZeroSum+=sobelVal;
            nonZeroSquaredSum+=sobelVal*sobelVal;

            edges_arr[i*row_size+j] = sobelVal;

        }


        
    }
    }

    for (int i=380;i<400;++i){
    for (int j=380;j<400;++j){
        std::cout << edges_arr[i*row_size + j] <<' ';
    }
    std::cout << '\n';}
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

    int32_t jmin = col_size;
    int32_t jmax = 0;

    int32_t imin = row_size;
    int32_t imax = 0;

    for (int32_t i=filter_offset+1;i<row_size-(filter_offset+1);++i)
    for (int32_t j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        if (edges_arr[i*row_size+j]>(meanSobelVal+stdSobel) || edges_arr[i*row_size+j] < meanSobelVal){
                edges_arr[i*row_size+j] = 0;
        }
        else { // edge found
            filteredSobelCount++;
            iAvg+=i;
            jAvg+=j;
            iStdAcc+=i*i;
            jStdAcc+=j*j;

            jmin = j<jmin?j:jmin; // getting min j index
            jmax = j>jmax?j:jmax; // getting max j index 

            imin = i<imin?i:imin;
            imax = i>imax?i:imax;
            }
    }



    iAvg = iAvg/filteredSobelCount;
    jAvg = jAvg/filteredSobelCount;
    float iStd = std::sqrt(static_cast<float>(iStdAcc/filteredSobelCount - iAvg)); // could be a new float var?
    float jStd = std::sqrt(static_cast<float>(jStdAcc/filteredSobelCount - jAvg));
    std::cout << "Zero values count after filtering: " <<filteredSobelCount <<'\n';
    std::cout << "Average i = " << iAvg <<"\nSTD i = "<< iStd <<'\n';
    std::cout << "Average j = " << jAvg <<"\nSTD j = "<< jStd <<'\n';
    std::cout << "imin = " << imin <<"; imax = "<< imax <<";icentre = " << (imax-imin)/2 <<'\n';
    std::cout << "jmin = " << jmin <<"; jmax = "<< jmax <<";jcentre = " << (jmax-jmin)/2 <<'\n';

    std::cout << "Overriding iAvg and jAvg values with edges min/max centre\n";
    iAvg = (imax+imin)/2;
    jAvg = (jmax+jmin)/2;
    std::cout << "Average i = " << iAvg <<"; STD i = "<< iStd <<'\n';
    std::cout << "Average j = " << jAvg <<"; STD j = "<< jStd <<'\n';
    std::cout << "Overriding iAvg and jAvg values with averaged 1361 \n";
    iAvg = (iAvg+1361)/2;
    jAvg = (jAvg+1361)/2;

    // The question is: which metric to use for hough transform - points being +=1 if non-zero while could also add up intensities 
    // (high intensity targets could produce anomalous signals that will overwhelm the transform)
    // long long test = iAvg - iStd*0.1;



    // const double searchBoxBoundsSTDCoef = 1.0/100.0;
    // long long iSearchSize = iStd*searchBoxBoundsSTDCoef;
    // long long jSearchSize = jStd*searchBoxBoundsSTDCoef;
    long long iSearchSize = 64; // has to be even otherwise not all values will be accesible
    long long jSearchSize = 64;
    float halfDiagLen = std::sqrt(static_cast<float>(iSearchSize*iSearchSize + jSearchSize*jSearchSize));
    std::vector<int>radii {1215,1216,1217,1218,1219,1267,1268,1269,1270,1271,1272,1273}; //1195,1200,1205,1210,1215,1220
    std::vector<uint64_t>searchBox(iSearchSize*jSearchSize*radii.size(),0);

    std::cout << "centre search box coordinates (i): " << iAvg - iSearchSize/2 << " to " << iAvg + iSearchSize/2 << "\niSearchSize (Window) = " << iSearchSize <<'\n';
    std::cout << "centre search box coordinates (j): " << jAvg - jSearchSize/2 << " to " << jAvg + jSearchSize/2 << "\njSearchSize (Window) = " << jSearchSize <<'\n';


    // for every pixel of the tomogram
    for (int i = filter_offset + 1; i < row_size - (filter_offset + 1); ++i)
        for (int j = filter_offset + 1; j < col_size - (filter_offset + 1); ++j) {
            // real edge, compute distance from every centre point

            if (edges_arr[i * row_size + j] != 0) {
                // check if the edge is not within 0.5 (0.7?) std from the avg centre

                float distanceToAvgCentre = std::sqrt(std::pow(static_cast<float>(i - iAvg), 2.0F) +
                                                      std::pow(static_cast<float>(j - jAvg), 2.0F));
                // for radius in radii should be here
                // iStd is ~= radius hence 0.5iStd is a half of that which is where
                // we do not expect the edges of the cyllinder to be
                if (distanceToAvgCentre < (iStd / 2.0F))
                    continue;
                for (size_t rIdx = 0; rIdx < radii.size(); rIdx++) {
                    if (static_cast<int>(distanceToAvgCentre - halfDiagLen) <= radii[rIdx] &&
                        static_cast<int>(distanceToAvgCentre + halfDiagLen) >= radii[rIdx]) {
                        // the edge could theoretically be part of the circle
                        // with the centre in the search box so we can continue
                        for (auto iCentre = iAvg - iSearchSize / 2;
                             iCentre < iAvg + iSearchSize / 2; ++iCentre)
                            for (auto jCentre = jAvg - jSearchSize / 2;
                                 jCentre < jAvg + jSearchSize / 2; ++jCentre) {

                                int distance = static_cast<int>(
                                    std::round(std::sqrt(std::pow(static_cast<float>(i - iCentre), 2.0F) +
                                              std::pow(static_cast<float>(j - jCentre), 2.0F))));
                                long long indexI = 0;
                                long long indexJ = 0;
                                int32_t edgeStd = 2;
                                if (distance >= (radii[rIdx]-edgeStd) && distance <= (radii[rIdx]+edgeStd)) {
                                    indexI = iCentre - iAvg + iSearchSize / 2;
                                    indexJ = jCentre - jAvg + jSearchSize / 2;
                                    searchBox.at(rIdx * (iSearchSize * jSearchSize) +
                                                 indexI * jSearchSize + indexJ) +=
                                        1; // searchBox[0]  // may be too much write accesses, could
                                           // be better to revise the loop so a tmp var could be
                                           // formed
                                }
                            }
                    }
                }
            }
        }


    // find the maximum pixel value
    size_t idxMax = 0;
    size_t valMax = 0;
    for(size_t i=1+iSearchSize;i<searchBox.size()-iSearchSize;++i){
        auto crossSum = searchBox.at(i-iSearchSize) + searchBox.at(i-1) + searchBox.at(i) + searchBox.at(i+1) + searchBox.at(i+iSearchSize);
        if (valMax < crossSum){
            valMax = crossSum;
            idxMax = i;
        }
    }
    std::cout << "MaxVal: " << valMax << "; MaxValIdx: " << idxMax <<'\n';

    std::ofstream myfile;
    myfile.open ("example.txt");
    
    std::cout << "Outputting the accumulator\n\n";
    for (size_t rIdx =0; rIdx<radii.size();rIdx++) {
        for(auto i = 0; i < iSearchSize;i++){
        for(auto j = 0; j < jSearchSize;j++){
            myfile << searchBox[(rIdx*(iSearchSize*jSearchSize) + i*jSearchSize + j)] <<' ';
        }
        myfile << '\n';
    }
    }
    myfile.close();

    // write out the resultant netcdf after sobel filtering?




    //TODO run hough transform for the edges left (non-zero) 
}