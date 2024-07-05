#pragma once
#include <cstdint>
#include <vector>

auto filterEdgesGrid(std::vector<float> &edges_arr, uint32_t row_size,
                     uint32_t col_size, uint32_t gridSide = 20,
                     uint32_t threshold = 100) -> std::pair<float, float>;

void smoothAvg(const std::vector<uint16_t> &arr, std::vector<uint16_t> &smooth_arr,
               const size_t row_size, const size_t col_size, const int32_t filterSide);