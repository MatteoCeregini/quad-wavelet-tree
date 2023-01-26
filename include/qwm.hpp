#pragma once

#include <vector>
#include <cstdint>
#include <cassert>
#include <bit>
#include <algorithm>
#include <iterator>


template<typename qvector, bool precomputed_paths = false>
class wavelet_matrix {
    private:
        uint64_t m_size;   // size of the text represented by the wavelet matrix.
        uint64_t m_max_level; // number of levels of the wavelet matrix
        qvector m_qvector; // qvector used to store the entire wavelet matrix, one level after the other.
        std::vector<std::vector<uint64_t>> m_rank; // for each level l, for each symbol s = 0..3, we store rank(s, l * m_size)
        std::vector<std::vector<uint64_t>> m_count; // for each level l, for each symbol s = 0..3, we store how many symbols < s there are at level l.
        // Used only when 'precomputed_paths = true'
        std::vector<std::vector<uint64_t>> m_path_off; // for each level l, for each value v = 0..sigma, we store ...
        std::vector<std::vector<uint64_t>> m_rank_path_off; // for each level l, for each value v = 0..sigma, we store ...
    
    private:
        // Counting sort used for the costruction of the wavelet matrix. 
        template <typename unsigned_int>
        std::vector<unsigned_int> counting_sort(std::vector<unsigned_int>& sequence, uint64_t shift) {
            
            std::vector<unsigned_int> sorted_sequence(sequence.size(), 0ULL);
            std::vector<uint64_t> count(4, 0ULL);

            for (auto value : sequence) {
                uint64_t symbol = (value >> shift) & 0b11;
                count[symbol]++;
            }
            for (uint64_t i = 1; i < 4; ++i) {
                count[i] += count[i - 1];
            }
            for (uint64_t i = sequence.size(); i > 0; --i) {
                uint64_t value = sequence[i - 1];
                uint64_t symbol = (value >> shift) & 0b11;
                count[symbol]--;
                sorted_sequence[count[symbol]] = value;
            }
            return sorted_sequence;
        }

    public:
        wavelet_matrix() = default;
        
        // Compute the wavelet matrix of a sequence of integers with alphabet = [0, sigma]
        template <typename unsigned_int>
        wavelet_matrix(std::vector<unsigned_int> sequence, unsigned_int sigma) {
            assert(sequence.size() != 0);

            m_size = sequence.size();
            // How many bits are needed to represent the largest value in the alphabet?
            uint64_t bits_sigma = std::bit_width(sigma);
            // If bits_sigma is odd, add one so it becomes even.
            if (bits_sigma % 2 != 0) bits_sigma++;
            // Since each symbol in the qvector is 2-bit long, we have bits_sigma / 2 levels.
            m_max_level = bits_sigma / 2;
            // Create a qvector containing all the levels of the matrix.
            m_qvector = qvector(m_size * m_max_level);
            uint64_t shift = 2 * (m_max_level - 1);
            
            uint64_t i = 0;
            for (uint64_t level = 0; level < m_max_level; ++level) {
                for (uint64_t value : sequence) {
                    uint64_t symbol = (value >> shift) & 0b11;
                    m_qvector.put(i, symbol);
                    ++i;
                }
                sequence = counting_sort(sequence, shift);
                /*
                std::stable_sort(sequence.begin(), sequence.end(), [&shift](uint64_t a, uint64_t b) {
                    return ((a >> shift) & 0b11) < ((b >> shift) & 0b11);
                });
                */
                shift -= 2;
            }
            // Initialize support for rank and select queries.
            m_qvector.init_rank_select_support();
            // Initialize precomputed info...
            m_rank.assign(m_max_level, std::vector<uint64_t>(4, 0ULL));
            m_count.assign(m_max_level, std::vector<uint64_t>(4, 0ULL));

            for (uint64_t level = 0; level < m_max_level; ++level) {
                for (uint64_t symbol = 0; symbol < 4; ++symbol) {
                    m_rank[level][symbol] = m_qvector.rank(level * m_size, symbol);
                    uint64_t cnt = 0;
                    for (uint64_t s = 0; s < symbol; ++s) {
                        uint64_t rank_level = m_qvector.rank((level + 1) * m_size, s) - m_qvector.rank(level * m_size, s);
                        cnt += rank_level;
                    }
                    m_count[level][symbol] = cnt;     
                }       
            }

            if constexpr (precomputed_paths) {
                m_path_off.assign(m_max_level + 1, std::vector<uint64_t>(sigma + 1, 0ULL));
                m_rank_path_off.assign(m_max_level + 1, std::vector<uint64_t>(sigma + 1, 0ULL));
                
                for (uint64_t value = 0; value <= sigma; ++value) {
                    uint64_t b = 0;
                    uint64_t shift = 2 * (m_max_level - 1);
                    for (uint64_t level = 0; level < m_max_level; ++level) {
                        uint64_t symbol = (value >> shift) & 0b11;
                        uint64_t rank_b = m_qvector.rank(b, symbol);
                        b = (level + 1) * m_size + (rank_b - m_rank[level][symbol]) + m_count[level][symbol];
                        shift -= 2;
                        m_path_off[level + 1][value] = b;
                        m_rank_path_off[level][value] = rank_b;
                    }
                }
            }
        }

        inline uint64_t operator[](uint64_t i) const {
            assert(i < m_size);

            uint64_t result = 0ULL;
            for (uint64_t level = 0; level < m_max_level - 1; ++level) {
                uint64_t symbol = m_qvector[i];
                result <<= 2;
                result |= symbol;
                i = (level + 1) * m_size + (m_qvector.rank(i, symbol) - m_rank[level][symbol]) + m_count[level][symbol];
            }
            // unrolling the last iteration: the compiler does not understand that updating 'i' is useless
            uint64_t symbol = m_qvector[i];
            result <<= 2;
            result |= symbol;
            return result;
        }
        /*
        inline uint64_t rank(uint64_t i, uint64_t value) const {
            assert(i <= m_size);

            uint64_t b = 0; // starting position of the interval
            uint64_t shift = 2 * (m_max_level - 1);
            for (uint64_t level = 0; level < m_max_level; ++level) {
                uint64_t symbol = (value >> shift) & 0b11;
                uint64_t rank_b = m_qvector.rank(b, symbol); // occurrences of 'symbol' in the interval [0,b)
                i = m_qvector.rank(b + i, symbol) - rank_b;  // occurrences of 'symbol' in the interval [b,i)
                b = (level + 1) * m_size + (rank_b - m_rank[level][symbol]) + m_count[level][symbol];
                shift -= 2;
            }
            return i;
        }
        */

        inline uint64_t rank(uint64_t i, uint64_t value) const {
            assert(i <= m_size);

            uint64_t b = 0; // starting position of the interval
            uint64_t shift = 2 * (m_max_level - 1);
            for (uint64_t level = 0; level < m_max_level - 1; ++level) {
                uint64_t symbol = (value >> shift) & 0b11;
                uint64_t rank_b = m_qvector.rank(b, symbol); // occurrences of 'symbol' in the interval [0,b)
                i = m_qvector.rank(b + i, symbol) - rank_b;  // occurrences of 'symbol' in the interval [b,i)
                b = (level + 1) * m_size + (rank_b - m_rank[level][symbol]) + m_count[level][symbol];
                shift -= 2;
            }
            // unrolling the last iteration: the compiler does not understand that updating 'b' is useless
            uint64_t symbol = (value >> shift) & 0b11;
            uint64_t rank_b = m_qvector.rank(b, symbol);
            i = m_qvector.rank(b + i, symbol) - rank_b;
            return i;
        }

        template <bool linear_search>
        inline uint64_t select(uint64_t i, uint64_t value) const {
            assert(i >= 1 && i <= rank(size(), value));

            if constexpr (precomputed_paths) {
                uint64_t shift = 0;
                for (uint64_t level = m_max_level; level > 0; --level) {
                    uint64_t symbol = (value >> shift) & 0b11;
                    i = m_qvector.template select<linear_search>(m_rank_path_off[level - 1][value] + i, symbol) - m_path_off[level - 1][value] + 1;
                    shift += 2;
                }
                return i - 1;
            } else {
                std::vector<uint64_t> path_off(m_max_level + 1, 0ULL);
                std::vector<uint64_t> rank_path_off(m_max_level + 1, 0ULL);
                uint64_t b = 0;
                uint64_t shift = 2 * (m_max_level - 1);
                for (uint64_t level = 0; level < m_max_level; ++level) {
                    uint64_t symbol = (value >> shift) & 0b11;
                    uint64_t rank_b = m_qvector.rank(b, symbol);
                    b = (level + 1) * m_size + (rank_b - m_rank[level][symbol]) + m_count[level][symbol];
                    shift -= 2;
                    path_off[level + 1] = b;
                    rank_path_off[level] = rank_b;
                }
                shift = 0; 
                for (uint64_t level = m_max_level; level > 0; --level) {
                    b = path_off[level - 1];
                    uint64_t rank_b = rank_path_off[level - 1];
                    uint64_t symbol = (value >> shift) & 0b11; 
                    i = m_qvector.template select<linear_search>(rank_b + i, symbol) - b + 1;
                    shift += 2;
                }
                return i - 1;
            }
        }

        inline uint64_t size() const {
            return m_size;
        }

        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(m_size);
            visitor.visit(m_max_level);
            visitor.visit(m_qvector);
            visitor.visit(m_rank);
            visitor.visit(m_count);
            visitor.visit(m_path_off);
            visitor.visit(m_rank_path_off);
        }
};