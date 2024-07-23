#pragma once
#include <cstddef>  // for size_t
#include <cstdint>  // for uint32_t, uint16_t, int32_t, uint64_t
#include <utility>  // for pair
#include <vector>   // for vector

using std::size_t;

auto filterEdgesGrid(std::vector<float> &edges_arr, uint32_t row_size,
                     uint32_t col_size, uint32_t gridSide = 20,
                     uint32_t threshold = 100) -> std::pair<float, float>;

void smoothAvg(const std::vector<uint16_t> &arr, std::vector<uint16_t> &smooth_arr,
               const size_t row_size, const size_t col_size, const int32_t filterSide);

std::vector<uint64_t> edgesAccumulate(std::vector<float>& edges_arr, size_t row_size,
                                      size_t col_size, std::vector<float>& radii, long long iCentreApprox,
                                      long long jCentreApprox, float distThreshold, long long iSearchSize,
                                      long long jSearchSize, long long scaleCoeff);