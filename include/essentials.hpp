#pragma once

#include <istream>
#include <fstream>
#include <ios>
#include <vector>
#include <chrono>
#include <map>
#include <cstdint>
#include <numeric>
#include <utility>
#include <random>
#include <algorithm>

namespace essentials {

// Represents a sequence of integers
template <typename unsigned_int>
struct dataset {
    std::vector<unsigned_int> m_sequence; // The text is encoded as a sequence of integers
    unsigned_int m_sigma; // Biggest value in the sequence

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_sequence);
        visitor.visit(m_sigma);
    }
};

// Reads a file and map the text into a sequence of integers
template <typename unsigned_int>
bool create_dataset(dataset<unsigned_int>& ds, const char* path) {
    ds.m_sequence.clear();
    ds.m_sigma = 0;
    
    std::ifstream file_stream;
    file_stream.open(path);
    if(!file_stream.is_open()) {
        return false;
    }
    char c;
    std::map<char, unsigned_int> encoding;

    while (file_stream.get(c)) {
        if (!encoding.contains(c)) {
            encoding.insert({c, ds.m_sigma});
            ds.m_sigma++;
        }
        unsigned_int enc = encoding.at(c);
        ds.m_sequence.push_back(enc);
    }
    file_stream.close();
    ds.m_sigma--;
    return true;
}

// Generates a random sequence, using a certain distribution
template <typename unsigned_int, typename distrib_type>
void generate_random_sequence(dataset<unsigned_int>& ds, uint64_t sequence_size, distrib_type& distrib) {
    std::random_device random_dev;
    std::mt19937 generator(random_dev());

    auto generate_random_value = [&distrib, &generator](){
        return distrib(generator);
    };

    ds.m_sequence.assign(sequence_size, 0ULL);
    std::generate(ds.m_sequence.begin(), ds.m_sequence.end(), generate_random_value);

    ds.m_sigma = *std::max_element(ds.m_sequence.begin(), ds.m_sequence.end());
}

// Generates legal access queries
template <typename unsigned_int>
std::vector<uint64_t> generate_access_queries(essentials::dataset<unsigned_int>& ds, uint64_t n_queries) {
    std::random_device random_dev;
    std::mt19937 generator(random_dev());
    std::uniform_int_distribution<uint64_t> distrib(0, ds.m_sequence.size() - 1);

    auto gen_random_query = [&distrib, &generator]() {
        return distrib(generator);
    };

    std::vector<uint64_t> queries(n_queries);
    std::generate(queries.begin(), queries.end(), gen_random_query);
    return queries;
}

// Generates legal rank queries
template <typename unsigned_int>
std::vector<std::pair<uint64_t, unsigned_int>> generate_rank_queries(essentials::dataset<unsigned_int>& ds, uint64_t n_queries) {
    std::vector<std::pair<uint64_t, unsigned_int>> queries(n_queries);

    std::random_device random_dev;
    std::mt19937 generator(random_dev());
    std::uniform_int_distribution<uint64_t> distrib(0, ds.m_sequence.size() - 1);

    for (uint64_t i = 0; i < n_queries; ++i) {
        uint64_t p = distrib(generator);
        unsigned_int v = ds.m_sequence[p];
        queries[i] = std::make_pair(p, v);
    } 
    return queries;
}

// Generates legal select queries
template <typename unsigned_int>
std::vector<std::pair<uint64_t, unsigned_int>> generate_select_queries(essentials::dataset<unsigned_int>& ds, uint64_t n_queries) {
    std::vector<uint64_t> occurrences(ds.m_sigma + 1, 0ULL);
    for (auto value : ds.m_sequence) {
        occurrences[value]++;
    }

    std::vector<std::pair<uint64_t, unsigned_int>> queries(n_queries);

    std::random_device random_dev;
    std::mt19937 generator(random_dev());
    std::uniform_int_distribution<uint64_t> access_distrib(0, ds.m_sequence.size() - 1);
    std::vector<std::uniform_int_distribution<uint64_t>> select_distributions(ds.m_sigma + 1);

    for (unsigned_int value = 0; value <= ds.m_sigma; ++value) {
        uint64_t occ_value = occurrences[value];
        std::uniform_int_distribution<uint64_t> d(1, occ_value);
        select_distributions[value] = d;
    }

    for (uint64_t i = 0; i < n_queries; ++i) {
        unsigned_int v = ds.m_sequence[access_distrib(generator)];
        uint64_t r = select_distributions[v](generator);
        queries[i] = std::make_pair(r, v);
    }
    return queries;
}

// Checks if a variable is aligned to n (in bytes)
bool is_aligned(const volatile void *p, std::size_t n) {
    return reinterpret_cast<std::uintptr_t>(p) % n == 0;
}

// Prevents the compiler to optimize away an expression
template <typename T>
inline void do_not_optimize_away(T&& value) {
    asm volatile("" : "+r"(value));
}

// Useful for measuring execution times
template <typename DurationType>
class timer {
    private:
        std::chrono::high_resolution_clock::time_point m_start_time;
        std::vector<double> m_timings;

    public:
        size_t runs() const {
            return m_timings.size();
        }

        void start() {
            m_start_time = std::chrono::high_resolution_clock::now();
        }

        double stop() {
            auto duration = std::chrono::duration_cast<DurationType>(std::chrono::high_resolution_clock::now() - m_start_time);
            m_timings.push_back(duration.count());
            return duration.count();
        }

        void reset() {
            m_timings.clear();
        }

        double elapsed() const {
            return std::accumulate(m_timings.begin(), m_timings.end(), 0.0);
        }

        double average() const {
            return elapsed() / runs();
        }  
};

// Given a vector, returns its size in bytes
template <typename T>
size_t vec_bytes(T const& vec) {
    return vec.size() * sizeof(vec.front()) + sizeof(typename T::size_type);
}

// Given a pod, returns its size in bytes
template <typename T>
size_t pod_bytes(T const& pod) {
    static_assert(std::is_standard_layout_v<T> && std::is_trivial_v<T>);
    return sizeof(pod);
}

template <typename T>
void load_pod(std::istream& is, T& val) {
    static_assert(std::is_standard_layout_v<T> && std::is_trivial_v<T>);
    is.read(reinterpret_cast<char*>(&val), sizeof(T));
}

template <typename T, typename Allocator>
void load_vector(std::istream& is, std::vector<T, Allocator>& vec) {
    size_t n;
    load_pod(is, n);
    vec.resize(n);
    is.read(reinterpret_cast<char*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

template <typename T>
void save_pod(std::ostream& os, T const& val) {
    static_assert(std::is_standard_layout_v<T> && std::is_trivial_v<T>);
    os.write(reinterpret_cast<char const*>(&val), sizeof(T));
}

template <typename T, typename Allocator>
void save_vector(std::ostream& os, std::vector<T, Allocator> const& vec) {
    static_assert(std::is_standard_layout_v<T> && std::is_trivial_v<T>);
    size_t n = vec.size();
    save_pod(os, n);
    os.write(reinterpret_cast<char const*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

template <typename T, typename Visitor>
size_t visit(T& data_structure, char const* filename) {
    Visitor visitor(filename);
    visitor.visit(data_structure);
    return visitor.bytes();
}

class loader {
    private:
        size_t m_num_bytes_pods;
        size_t m_num_bytes_vecs_of_pods;
        std::ifstream m_is;
    public:
        loader(char const* filename) : 
            m_num_bytes_pods(0),
            m_num_bytes_vecs_of_pods(0), 
            m_is(filename, std::ios::binary) {
            
            if (!m_is.good()) {
                throw std::runtime_error("Error in opening binary file!");
            }
        }

        ~loader() {
            m_is.close();
        }

        template <typename T>
        void visit(T& val) {
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                load_pod(m_is, val);
                m_num_bytes_pods += pod_bytes(val);
            } else {
                val.visit(*this);
            }
        }

        template <typename T, typename Allocator>
        void visit(std::vector<T, Allocator>& vec) {
            size_t n;
            visit(n);
            vec.resize(n);
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                m_is.read(reinterpret_cast<char*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
                m_num_bytes_vecs_of_pods += n * sizeof(T);
            } else {
                for (auto& v : vec) visit(v);
            }
        }

        size_t bytes() {
            return m_is.tellg();
        }

        size_t bytes_pods() {
            return m_num_bytes_pods;
        }

        size_t bytes_vecs_of_pods() {
            return m_num_bytes_vecs_of_pods;
        }
};


class saver {
    private:
        std::ofstream m_os;
    public:
        saver(char const* filename)
            : m_os(filename, std::ios::binary) {
            if (!m_os.good()) {
                throw std::runtime_error("Error in opening binary file!");
            }
        }

        ~saver() {
            m_os.close();
        }

        template <typename T>
        void visit(T& val) {
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                save_pod(m_os, val);
            } else {
                val.visit(*this);
            }
        }

        template <typename T, typename Allocator>
        void visit(std::vector<T, Allocator>& vec) {
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                save_vector(m_os, vec);
            } else {
                size_t n = vec.size();
                visit(n);
                for (auto& v : vec) visit(v);
            }
        }

        size_t bytes() {
            return m_os.tellp();
        }
};


class sizer {
    private:
        struct node {
            size_t bytes;
            size_t depth;
            std::vector<node> children;

            node(size_t b, size_t d) : bytes(b), depth(d) {}
        };
        node m_root;
        node* m_current;

    public:
        sizer() : m_root(0, 0) , m_current(&m_root) {}

        template <typename T>
        void visit(T& val) {
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                node n(pod_bytes(val), m_current->depth + 1);
                m_current->children.push_back(n);
                m_current->bytes += n.bytes;
            } else {
                val.visit(*this);
            }
        }

        template <typename T, typename Allocator>
        void visit(std::vector<T, Allocator>& vec) {
            if constexpr (std::is_standard_layout_v<T> && std::is_trivial_v<T>) {
                node n(vec_bytes(vec), m_current->depth + 1);
                m_current->children.push_back(n);
                m_current->bytes += n.bytes;
            } else {
                size_t n = vec.size();
                m_current->bytes += pod_bytes(n);
                node* parent = m_current;
                for (auto& v : vec) {
                    node n(0, parent->depth + 1);
                    parent->children.push_back(n);
                    m_current = &parent->children.back();
                    visit(v);
                    parent->bytes += m_current->bytes;
                }
                m_current = parent;
            }
        }

        size_t bytes() const {
            return m_root.bytes;
        }
};


template <typename T>
size_t load(T& data_structure, char const* filename) {
    return visit<T, loader>(data_structure, filename);
}

template <typename T>
size_t save(T& data_structure, char const* filename) {
    return visit<T, saver>(data_structure, filename);
}

template <typename T>
size_t bytes(T& data_structure) {
    sizer visitor;
    visitor.visit(data_structure);
    return visitor.bytes();
}

} // namespace essentials