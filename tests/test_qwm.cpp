#include <cstdint>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <map>
#include <cassert>
#include <bit>

#include "essentials.hpp"
#include "qvector.hpp"
#include "qwm.hpp"


#define PRECOMPUTED_PATHS true

int main(int argc, char** argv)  {

    if (argc != 2) {
        std::cout << "USAGE: ./executable PATH_TO_DATASET" << std::endl;
        return -1;
    }

    const char* path_to_dataset = argv[1];

    if (!std::filesystem::exists(path_to_dataset)) {
        std::cout << "ERROR: \'" << path_to_dataset << "\' does not exists!" << std::endl;
        return -1;
    }

    std::cout << "[1] Loading the dataset..." << std::endl;
    essentials::dataset<uint16_t> ds;
    essentials::load(ds, path_to_dataset);

    std::cout << "[2] Creating the quaternary wavelet matrices..." << std::endl;
    wavelet_matrix<qvector<512>, PRECOMPUTED_PATHS> wm_512(ds.m_sequence, ds.m_sigma);
    wavelet_matrix<qvector<256>, PRECOMPUTED_PATHS> wm_256(ds.m_sequence, ds.m_sigma);


    std::cout << "[3] testing access queries..." << std::endl;
    for (uint64_t i = 0; i < ds.m_sequence.size(); ++i) {
        if (ds.m_sequence[i] != wm_512[i]) {
            std::cout << "[ERROR (WM 512 symbols basic block)] access(" << i << ") = " << wm_512[i] << " but it should be: " << uint64_t(ds.m_sequence[i]) << std::endl;
            return -1;
        }
        if (ds.m_sequence[i] != wm_256[i]) {
            std::cout << "[ERROR (WM 256 symbols basic block)] access(" << i << ") = " << wm_256[i] << " but it should be: " << uint64_t(ds.m_sequence[i]) << std::endl;
            return -1;
        }
    }

    std::cout << "[4] testing rank queries..." << std::endl;
    std::vector<uint64_t> ranks(ds.m_sigma + 1, 0);
    for(uint64_t i = 0; i <= wm_256.size(); ++i) {
        for (uint64_t value = 0; value <= ds.m_sigma; ++value) {
            if (ranks[value] != wm_512.rank(i, value)) {
                std::cout << "[ERROR (WM 512 symbols basic block)] rank(" << i << ", " << value << ") = " << wm_512.rank(i, value) << " but it should be: " << ranks[value] << std::endl;
                return -1;
            }
            if (ranks[value] != wm_256.rank(i, value)) {
                std::cout << "[ERROR (WM 256 symbols basic block)] rank(" << i << ", " << value << ") = " << wm_256.rank(i, value) << " but it should be: " << ranks[value] << std::endl;
                return -1;
            }
        }
        if (i < wm_256.size()) {
            ranks[wm_256[i]]++;
        }
    }

    std::cout << "[5] testing select queries..." << std::endl;
    std::vector<uint64_t> select_results[ds.m_sigma + 1];
    // Save the result of all legal select queries for each value
    for (uint64_t i = 0; i < ds.m_sequence.size(); ++i) {
        uint64_t value = ds.m_sequence[i];
        select_results[value].push_back(i); // rank = i is at position i-1
    }

    for (uint64_t value = 0; value <= ds.m_sigma; ++value) {
        for (uint64_t r = 1; r <= wm_256.rank(wm_256.size(), value); ++r) {
            if (wm_512.select<true>(r, value) != select_results[value][r - 1]) {
                std::cout << "[ERROR (WM 512 symbols basic block), linear search] select(" << r << ", " << value << ") = " << wm_512.select<true>(r, value) << " but it should be: " << select_results[value][r - 1] << std::endl;
                return -1;
            }
            if (wm_256.select<false>(r, value) != select_results[value][r - 1]) {
                std::cout << "ERROR (WM 256 symbols basic block), binary search] select(" << r << ", " << value << ") = " << wm_256.select<false>(r, value) << " but it should be: " << select_results[value][r - 1] << std::endl;
                return -1;
            }
        }
    }
    return 0;
}