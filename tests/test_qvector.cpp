#include <iostream>
#include <random>
#include <cstdint>
#include <algorithm>

#include "essentials.hpp"
#include "qvector.hpp"

std::vector<uint8_t> generate_random_text(uint64_t text_size) {
    std::vector<uint8_t> text(text_size);
    // Create a random uniform distribution in range [0, 3].
    std::random_device device;
    std::mt19937 rng(device());
    std::uniform_int_distribution<std::mt19937::result_type> distr(0, 3);
    // Generate the random text.
    for (uint64_t i = 0; i < text_size; ++i) {
        text[i] = distr(rng);
    }
    return text;
}

template <typename qvector_type> 
bool test_access(qvector_type& qv, std::vector<uint8_t>& text) {
    for (uint64_t i = 0; i < text.size(); ++i) {
        uint64_t symbol = text[i];
        if (qv[i] != symbol) {
            std::cout << "[ERROR QVECTOR] access(" << i << ") = " << qv[i] << " but it should be: " << symbol << std::endl;
            return false;
        }
    }
    return true;
}

template <typename qvector_type> 
bool test_rank(qvector_type& qv, std::vector<uint8_t>& text) {
    uint64_t ranks[4]{0, 0, 0, 0};
    // For each position of the qvector, check if the rank queries for that position are correct.
    for (uint64_t i = 0; i <= qv.size(); ++i) {
        for (uint64_t symbol = 0; symbol < 4; ++symbol) {
            if (qv.rank(i, symbol) != ranks[symbol]) {
                std::cout << "[ERROR QVECTOR] rank(" << i << ", " << symbol << ") = " << qv.rank(i, symbol) << " but it should be: " << ranks[symbol] << std::endl;
                return false;
            }       
        }
        if (i < qv.size()) {
            uint64_t symbol = qv[i];
            ranks[symbol]++;
        }
    }
    return true;
}

template <typename qvector_type> 
bool test_select(qvector_type& qv, std::vector<uint8_t>& text) {
    std::vector<uint64_t> select_results[4];
    // Save the result of all legal select queries for each symbol
    for (uint64_t i = 0; i < qv.size(); ++i) {
        uint64_t symbol = qv[i];
        select_results[symbol].push_back(i); // rank = i is at position i-1
    }

    for (uint64_t symbol = 0; symbol < 4; ++symbol) {
        for (uint64_t r = 1; r <= qv.rank(qv.size(), symbol); ++r) {
            if (qv.template select<false>(r, symbol) != select_results[symbol][r - 1]) {
                std::cout << "[ERROR SELECT (binary search)] " << "select(" << r << ", " << symbol << ") = " << qv.template select<false>(r, symbol) << " vs " << select_results[symbol][r - 1] << std::endl;
                return false;
            }
            if (qv.template select<true>(r, symbol) != select_results[symbol][r - 1]) {
                std::cout << "[ERROR SELECT (linear search)] " << "select(" << r << ", " << symbol << ") = " << qv.template select<true>(r, symbol) << " vs " << select_results[symbol][r - 1] << std::endl;
                return false;
            }
        }
    }
    return true;
}


#define BLOCK_SIZE 256 // or 512

int main() {   

    uint64_t text_size = 1e9;

    std::cout << "[1] generating random text... " << std::endl; 
    auto text = generate_random_text(text_size);

    std::cout << "[2] creating the qvector... " << std::endl;
    // Create the qvector.
    qvector<BLOCK_SIZE> qv(text.size());

    for (uint64_t i = 0; i < text.size(); ++i) {
        uint64_t symbol = text[i];
        qv.put(i, symbol);
    }

    qv.init_rank_select_support();
    
    std::cout << "[3] testing access queries... " << std::endl;
    if (test_access(qv, text)) {
        std::cout << "    all queries are correct!" << std::endl;
    } else {
        return -1;
    }
    
    std::cout << "[4] testing rank queries... " << std::endl;
    if (test_rank(qv, text)) {
        std::cout << "    all queries are correct!" << std::endl;
    } else {
        return -1;
    }
    
    
    std::cout << "[4] testing select queries... " << std::endl;
    if (test_select(qv, text)) {
        std::cout << "    all queries are correct!" << std::endl;
    } else {
        return -1;
    }
    return 0;
}