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

    // auto filepath = "/home/bogdanm/data/containerSamples/RSES_Wood_PigTeeth_3rdMolars/tomoSliceZ-2__RMG.nc"; // inner: 1208, outer: 1271(?)
    auto filepath = "/home/bogdanm/data/containerSamples/RSES_Wood_Teeth_123_8mm/tomoSliceZ-7__R.nc"; // inner: 1097
    // auto filepath = "/home/bogdanm/data/containerSamples/Whiting_5640_5mm_031114_preserved/tomoSliceZ-13__R.nc"; // inner: 1088
    // auto filepath = "/home/bogdanm/data/containerSamples/Whiting_5640_5mm_031114_Xe_dec/tomoSliceZ-12__R.nc";
    // std::vector<int>radii {1215,1216,1217,1218,1219,1267,1268,1269,1270,1271,1272,1273}; //1195,1200,1205,1210,1215,1220
    // std::vector<int>radii {1207,1208,1209,};

    std::vector<int>radii {1097};

    // std::vector<int>radii {868,869,870,};



    netCDF::NcFile dataFile(filepath, netCDF::NcFile::read);

    netCDF::NcVar tomodata = dataFile.getVar("tomo", netCDF::NcGroup::Current);

    std::cout << "there are " << tomodata.getDimCount() << " dimensions in the tomo dataset\n";
    std::cout << "tomodata.getDim(0).getSize() " << tomodata.getDim(0).getSize() << "\n";
    std::cout << "tomodata.getDim(1).getSize() " << tomodata.getDim(1).getSize() << "\n";
    std::cout << "tomodata.getDim(2).getSize() " << tomodata.getDim(2).getSize() << "\n";



    size_t z_size   = tomodata.getDim(0).getSize();
    size_t row_size = tomodata.getDim(1).getSize();
    size_t col_size = tomodata.getDim(2).getSize();


    const std::vector<size_t> &start{0, 0, 0};           // 3*4, 12 bytes
    const std::vector<size_t> &count{z_size, row_size, col_size}; // 12 bytes

    std::vector<int16_t> vec_int16(row_size * col_size);

    // uint16_t data_out[8364144/2] {}; //check if 8192, 8388608 the actual pass is 8364032, but the limit is inconsistent
    // e.g. uint8_t data[8377999] {}; might pass and might throw a Segmentation fault

    tomodata.getVar(start, count, vec_int16.data()); // the data is in, it's a 2D array

    auto arr = *reinterpret_cast<std::vector<uint16_t>*>(&vec_int16);
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

    const size_t avg_filter_size = 5;
    const size_t filter_offset = avg_filter_size / 2;
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

    auto meanVal = static_cast<float>(nonZeroSum/countNonZero);
    auto meanSumSquared = static_cast<float>(nonZeroSquaredSum/countNonZero);

    auto stdVal = std::sqrt(static_cast<float>(meanSumSquared - meanVal*meanVal)); // approx std. good enough

    std::cout << "STD value: " << stdVal << '\n' << "Filter bounds:" << meanVal-stdVal << " to " << meanVal+stdVal << '\n';

    //TODO filter the edges, potentially with a hard cutoff
    long long filteredSobelCount = 0;
    long long iAvg {};
    long long jAvg {};
    long long iStdAcc {};
    long long jStdAcc {};

    // compute average edge i and j, and respective standard deviations. Those will be our center coordinate and search box dimensions.

    size_t jmin = col_size;
    size_t jmax = 0;

    size_t imin = row_size;
    size_t imax = 0;

    for (size_t i=filter_offset+1;i<row_size-(filter_offset+1);++i)
    for (size_t j=filter_offset+1;j<col_size-(filter_offset+1);++j){
        if (edges_arr[i*row_size+j]>(meanVal+stdVal) || edges_arr[i*row_size+j] < meanVal){
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


        //-----------------------------------------------------------------
    // doing subgrid filtering, let's make 20x20 window a standard
    // edges_arr is an input

    // row_size
    size_t gridSide = 20;
    size_t threshold = 100;
    size_t nonZeroValsCount = 0;
    nonZeroSum = 0.0; // already defined above
    nonZeroSquaredSum = 0.0; // already defined above

    std::cout << "Additional Filtering\n\n\n "<< row_size/gridSide << ' ' << col_size/gridSide << ' ' <<'\n';
    for(size_t i=0;i<row_size;i+=gridSide)
        for(size_t j=0;j<col_size;j+=gridSide){

            size_t nonZeroValsCountGrid = 0;
            float nonZeroSumGrid = 0.0F;
            float nonZeroSquaredSumGrid = 0.0F;

            for(size_t k=0;k<gridSide;k++)
            for(size_t l=0;l<gridSide;l++){
                auto idx = (i+k)*row_size+j + l;
                
                if(edges_arr[idx]!= 0){
                    nonZeroValsCountGrid ++;
                    nonZeroSumGrid+=edges_arr[idx];
                    nonZeroSquaredSumGrid+=edges_arr[idx]*edges_arr[idx];
                }
                
            }
            if (nonZeroValsCountGrid >= threshold){
                for(size_t k=0;k<gridSide;k++)
                for(size_t l=0;l<gridSide;l++){
                    edges_arr.at(i*row_size+j + k*row_size + l) = 0;
                }
            }else{
                nonZeroValsCount += nonZeroValsCountGrid;
                nonZeroSum+=nonZeroSumGrid;
                nonZeroSquaredSum+=nonZeroSquaredSumGrid;
            }

        }
    std::cout << "2. NonZero values count: " << nonZeroValsCount << '\n';
    std::cout << "2. Zero values count: " << (row_size * col_size) - nonZeroValsCount << '\n';
    std::cout << "2. Average value: " << nonZeroSum/static_cast<double>(nonZeroValsCount) << '\n';

    meanVal = static_cast<float>(nonZeroSum/nonZeroValsCount);
    meanSumSquared = static_cast<float>(nonZeroSquaredSum/nonZeroValsCount);

    stdVal = std::sqrt(static_cast<float>(meanSumSquared - meanVal*meanVal)); // approx std. good enough

    std::cout << "2. STD value: " << stdVal << '\n' << "Filter bounds:" << meanVal-stdVal << " to " << meanVal+stdVal << '\n';

    for(auto& val:edges_arr){
        val = (val<(meanVal) || val>(meanVal+stdVal))?0:val;
    }
        //     if (edges_arr[i*row_size+j]>(meanVal+stdVal) || edges_arr[i*row_size+j] < meanVal){
        //         edges_arr[i*row_size+j] = 0;
        // }


    // need to filter again! but before filtering, one needs mean and std

    std::ofstream myedges;
    myedges.open ("output_edges.txt");
    
    std::cout << "Outputting the edges\n\n";
        for(size_t i = 0; i < row_size;i++){
        for(size_t j = 0; j < col_size;j++){
            myedges << edges_arr[i*row_size+j] <<' ';
        }
        myedges << '\n';
    }
    myedges.close();

//  0 .. 20 .. 40 ..... rowsize//20
// .
// .
// 20
// .
// .
// colsize//20
    //-----------------------------------------------------------------




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
    // std::cout << "Overriding iAvg and jAvg values with averaged 1361 \n";
    // iAvg = (iAvg+static_cast<long long> (row_size)/2)/2;
    // jAvg = (jAvg+static_cast<long long>(col_size)/2)/2;

    // The question is: which metric to use for hough transform - points being +=1 if non-zero while could also add up intensities 
    // (high intensity targets could produce anomalous signals that will overwhelm the transform)
    // long long test = iAvg - iStd*0.1;



    // const double searchBoxBoundsSTDCoef = 1.0/100.0;
    // long long iSearchSize = iStd*searchBoxBoundsSTDCoef;
    // long long jSearchSize = jStd*searchBoxBoundsSTDCoef;
    long long iSearchSize = 64; // has to be even otherwise not all values will be accesible
    long long jSearchSize = 64;
    float halfDiagLen = std::sqrt(static_cast<float>(iSearchSize*iSearchSize + jSearchSize*jSearchSize));

    std::vector<uint64_t>searchBox(iSearchSize*jSearchSize*radii.size(),0);

    std::cout << "centre search box coordinates (i): " << iAvg - iSearchSize/2 << " to " << iAvg + iSearchSize/2 << "\niSearchSize (Window) = " << iSearchSize <<'\n';
    std::cout << "centre search box coordinates (j): " << jAvg - jSearchSize/2 << " to " << jAvg + jSearchSize/2 << "\njSearchSize (Window) = " << jSearchSize <<'\n';


    // for every pixel of the tomogram
    for (int i = filter_offset + 1; i < static_cast<int>(row_size - (filter_offset + 1)); ++i)
        for (int j = filter_offset + 1; j < static_cast<int>(col_size - (filter_offset + 1)); ++j) {
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
                                int32_t edgeStd = 0;
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
    // for(size_t i=1+iSearchSize;i<searchBox.size()-iSearchSize;++i){
    //     auto crossSum = searchBox.at(i-iSearchSize) + searchBox.at(i-1) + searchBox.at(i) + searchBox.at(i+1) + searchBox.at(i+iSearchSize);
    //     if (valMax < crossSum){
    //         valMax = crossSum;
    //         idxMax = i;
    //     }
    // }
    for(size_t i=1+iSearchSize;i<searchBox.size()-iSearchSize;++i){
        // auto crossSum = searchBox.at(i-iSearchSize) + searchBox.at(i-1) + searchBox.at(i) + searchBox.at(i+1) + searchBox.at(i+iSearchSize);
        if (valMax < searchBox.at(i)){
            valMax = searchBox.at(i);
            idxMax = i;
        }
    }


    std::cout << "MaxVal: " << valMax << "; MaxValIdx: " << idxMax <<'\n';
    std::cout << "Radius: " << radii.at(idxMax/(iSearchSize*iSearchSize)) << "; iIdx: " << (idxMax%(iSearchSize*iSearchSize))/iSearchSize
    << "; jIdx: " << idxMax%iSearchSize<<'\n';
    //get the actual coordinates
    std::cout << "xAbs: " << idxMax%iSearchSize - iSearchSize/2 + jAvg << 
    "; yAbs: " << (idxMax%(iSearchSize*iSearchSize))/iSearchSize - iSearchSize/2 + iAvg <<'\n';
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