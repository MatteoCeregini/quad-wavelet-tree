#pragma once

#include <iostream>
#include <cstdint>
#include <vector>
#include <bit>
#include <array>
#include <cassert>
#include <algorithm>

#include "select.hpp"
#include "cache_aligned_allocator.hpp"


class counters {
    private:
        alignas(64) __uint128_t m_counters[4];
    public:
        inline uint64_t get_block_counter(uint64_t symbol, uint64_t block_number) const {
            __uint128_t c = m_counters[symbol] << 12;
            return (c >> (block_number * 12)) & 0b111111111111;
        }

        inline uint64_t get_superblock_counter(uint64_t symbol) const {
            return m_counters[symbol] >> 84;
        }

        inline void set_block_counter(uint64_t symbol, uint64_t block_number, uint64_t counter) {
            if (block_number > 0) {
                block_number--;
                __uint128_t mask = ~(__uint128_t(4095) << (block_number * 12));
                m_counters[symbol] &= mask;
                m_counters[symbol] |= __uint128_t(counter) << (block_number * 12);
            }
        }

        inline void set_superblock_counter(uint64_t symbol, uint64_t counter) {
            // Clear the last 44 bits
            __uint128_t mask = __uint128_t(-1) >> 44;
            m_counters[symbol] &= mask;
            // Set the counter
            m_counters[symbol] |= __uint128_t(counter) << 84 ;
        }

        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(m_counters);
        }
};


template <uint64_t block_size = 256, uint64_t select_sampling_rate = 8192>
class qvector {
    private:
        static const uint64_t blocks_in_superblock = 8;
        static const uint64_t words_in_block = block_size / 32;
        static const uint64_t words_in_superblock = blocks_in_superblock * words_in_block;
        static const uint64_t superblock_size = blocks_in_superblock * block_size;

        uint64_t m_size;
        std::vector<uint64_t, cache_aligned_allocator<uint64_t>> m_bits;
        std::vector<counters> m_counters;
        std::vector<std::vector<uint32_t>> m_select_samples;

        inline uint64_t block_number(uint64_t i) const {
            return (i / block_size) % blocks_in_superblock;
        }

        inline uint64_t superblock_index(uint64_t i) const {
            return i / superblock_size;
        }

        template <uint64_t symbol>
        inline uint64_t block_rank(uint64_t i) const {
            static_assert(symbol <= 3);
            
            uint64_t result = 0;
            uint64_t begin_block = (i / block_size) * words_in_block;

            i %= block_size;
            uint64_t k = 0;
            for (; k < i / 32; ++k) {
                uint64_t bits = m_bits[begin_block + k];
                uint64_t occ;
                if constexpr (symbol == 0) occ = ~bits & (~bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 1) occ =  bits & (~bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 2) occ = ~bits & ( bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 3) occ =  bits & ( bits >> 1) & 0x5555555555555555ULL;
                result += std::popcount(occ);
            }
            i %= 32;
            if (i != 0) {
                uint64_t bits = m_bits[begin_block + k];
                uint64_t occ;
                uint64_t mask = 0x5555555555555555ULL >> (64 - (i * 2));
                if constexpr (symbol == 0) occ = ~bits & (~bits >> 1) & mask;
                if constexpr (symbol == 1) occ =  bits & (~bits >> 1) & mask;
                if constexpr (symbol == 2) occ = ~bits & ( bits >> 1) & mask;
                if constexpr (symbol == 3) occ =  bits & ( bits >> 1) & mask;
                result += std::popcount(occ);
            }
            return result;
        }

        template <uint64_t symbol>
        inline uint64_t block_select(uint64_t r, uint64_t begin_block) const {
            static_assert(symbol <= 3);
            assert(r >=1 && r <= block_size);

            uint64_t cnt = 0;
            uint64_t prev_cnt = 0;
            uint64_t occ = 0;
            int k = 0;

            for (; k < words_in_block && cnt < r; ++k) {
                prev_cnt = cnt;
                uint64_t bits = m_bits[begin_block + k];
                if constexpr (symbol == 0) occ = ~bits & (~bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 1) occ =  bits & (~bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 2) occ = ~bits & ( bits >> 1) & 0x5555555555555555ULL;
                if constexpr (symbol == 3) occ =  bits & ( bits >> 1) & 0x5555555555555555ULL;
                cnt += std::popcount(occ); 
            }
            k--;
            r -= prev_cnt;
            uint64_t result = k * 32;
            if (r != 0) result += (select_uint64_t(occ, r - 1) / 2);
            return result;
        }

    public:
        qvector() = default;
        
        qvector(const qvector&) = default;
        qvector(qvector&&) = default;

        qvector& operator=(const qvector&) = default;
        qvector& operator=(qvector&&) = default;

        qvector(uint64_t size) : m_size(size) {
            uint64_t n = ((m_size / block_size) + 1) * words_in_block;
            m_bits.assign(n, 0ULL);
        }

        inline void put(uint64_t i, uint64_t symbol) {
            assert(symbol <= 3);
            assert(i < m_size);

            uint64_t word_index = (i % 32) * 2;
            uint64_t vector_index = i / 32;
            // Clear the position of the symbol
            uint64_t mask = ~(uint64_t(0b11) << word_index);
            m_bits[vector_index] &= mask;
            // Set the symbol
            m_bits[vector_index] |= (symbol << word_index);
        }

        inline uint64_t size() const {
            return m_size;
        }

        void init_rank_select_support() {
            // Create an empty vector of samples for each symbol.
            m_select_samples.assign(4, std::vector<uint32_t>());

            uint64_t n_superblocks = (m_size / superblock_size) + 1;
            m_counters.assign(n_superblocks, counters());

            std::array<uint64_t, 4> superblock_counters{0, 0, 0, 0};
            std::array<uint64_t, 4> block_counters{0, 0, 0, 0};

            for (uint64_t i = 0; i <= m_size; ++i) {
                if (i % superblock_size == 0) {
                    for (uint64_t symbol = 0; symbol < 4; ++symbol) {
                        m_counters[superblock_index(i)].set_superblock_counter(symbol, superblock_counters[symbol]);
                    }
                    block_counters.fill(0);
                }
                if (i % block_size == 0) {
                    for (uint64_t symbol = 0; symbol < 4; ++symbol) {
                        m_counters[superblock_index(i)].set_block_counter(symbol, block_number(i), block_counters[symbol]);
                    }
                }
                if (i < m_size) {
                    uint64_t symbol = this->operator[](i);
                    if (superblock_counters[symbol] % select_sampling_rate == 0) {
                        m_select_samples[symbol].push_back(superblock_index(i));
                    }
                    superblock_counters[symbol]++;
                    block_counters[symbol]++;
                }
            }
            // We set the unused counters in the last superblock to 4095 (all 12 bits to one), in this way the select operation is easier.
            for (uint64_t i = block_number(m_size) + 1; i < blocks_in_superblock; ++i) {
                for (uint64_t symbol = 0; symbol < 4; ++symbol) {
                    m_counters[superblock_index(m_size)].set_block_counter(symbol, i, uint64_t(4095));
                }
            }
            for (uint64_t symbol = 0; symbol < 4; ++symbol) {
                m_select_samples[symbol].shrink_to_fit();
            }
        }

        inline uint64_t operator[](uint64_t i) const {
            assert(i < m_size);

            uint64_t word_offset = (i % 32) * 2;
            uint64_t word_index = i / 32;
            return (m_bits[word_index] >> word_offset) & 0b11;
        }

        inline uint64_t rank(uint64_t i, uint64_t symbol) const {
            assert(symbol <= 3);
            assert(i <= m_size);

            uint64_t sb_index = superblock_index(i);
            __builtin_prefetch(&m_counters[sb_index]);

            uint64_t result = 0;
            switch(symbol) {
                case 0:
                    result = block_rank<0>(i);
                    break;
                case 1:
                    result = block_rank<1>(i);
                    break;
                case 2:
                    result = block_rank<2>(i);
                    break;
                case 3:
                    result = block_rank<3>(i);
                    break;
            }
            result += m_counters[sb_index].get_superblock_counter(symbol);
            result += m_counters[sb_index].get_block_counter(symbol, block_number(i));
            return result;
        }

        template <bool linear_search>
        inline uint64_t select(uint64_t r, uint64_t symbol) const {
            assert(symbol <= 3);
            assert(r >= 1 && r <= rank(m_size, symbol));

            uint64_t result = 0;

            uint64_t sample_position = (r - 1) / select_sampling_rate;
            uint64_t sb_index_sample = m_select_samples[symbol][sample_position];
            
            auto start = m_counters.begin() + sb_index_sample;
            auto finish = m_counters.end();
            auto p = start;

            if constexpr (linear_search) {
                sample_position++;
                if (sample_position < m_select_samples[symbol].size()) { 
                    // The '+1' is needed because we want to include the next sample's superblock
                    finish = m_counters.begin() + m_select_samples[symbol][sample_position] + 1;
                }

                auto distance = std::distance(start, finish);
                uint64_t step = std::sqrt(distance);

                while (distance >= step) {
                    if (p->get_superblock_counter(symbol) >= r) break;
                    std::advance(p, step);
                    distance -= step;
                }
                if (std::distance(m_counters.begin(), p) >= step) {
                    p -= step;
                } else {
                    p = m_counters.begin();
                }           
                for (; p != finish; ++p) {
                    if (p->get_superblock_counter(symbol) >= r) break;
                }
            } else {
                sample_position++;
                if (sample_position < m_select_samples[symbol].size()) {
                    // The '+1' is needed because we want to include the next sample's superblock
                    finish = m_counters.begin() + m_select_samples[symbol][sample_position] + 1;
                }
                p = std::lower_bound(start, finish, r, [symbol](const counters& c, const uint64_t value) {
                    return c.get_superblock_counter(symbol) < value;
                });
            }
            if (p != m_counters.begin()) --p; 

            uint64_t sb_index = std::distance(m_counters.begin(), p);
            result += (sb_index * superblock_size);
            r -= p->get_superblock_counter(symbol);

            uint64_t b_number = 1;
            for (; b_number < blocks_in_superblock; ++b_number) {
                if (p->get_block_counter(symbol, b_number) >= r) break; 
            }
            b_number--;
            
            r -= p->get_block_counter(symbol, b_number);
            result += (b_number * block_size);
            if (r == 0) return result;

            uint64_t begin_block = (sb_index * words_in_superblock) + (b_number * words_in_block);

            if (symbol == 0) return result + block_select<0>(r, begin_block);
            if (symbol == 1) return result + block_select<1>(r, begin_block);
            if (symbol == 2) return result + block_select<2>(r, begin_block);
            return result + block_select<3>(r, begin_block);
        }

        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(m_size);
            visitor.visit(m_bits);
            visitor.visit(m_counters);
            visitor.visit(m_select_samples);
        }
};