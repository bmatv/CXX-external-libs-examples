#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

int main() {
    std::vector<int> a{1, 2, 3, 4, 5, 6, 1, 2, 3, 4};
    for (auto& element : a) {
        std::cout << element << ' ';
    }
    std::cout << '\n';

    auto b = std::max_element(a.begin(), a.end());

    std::cout << "Max value is: " << *b << '\n';
    std::cout << "Max value idx is: " << std::distance(a.begin(), b) << '\n';
}
