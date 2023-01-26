#pragma once

#include <limits>
#include <new>

template <typename element_type>
struct cache_aligned_allocator {
    using value_type = element_type;
    static std::align_val_t constexpr ALIGNMENT{64};

    [[nodiscard]] element_type* allocate(std::size_t n_elements) {
        if (n_elements > std::numeric_limits<std::size_t>::max() / sizeof(element_type)) {
            throw std::bad_array_new_length();
        }
        auto const n_bytes = n_elements * sizeof(element_type);
        return reinterpret_cast<element_type*>(::operator new[](n_bytes, ALIGNMENT));
    }

    void deallocate(element_type* ptr, [[maybe_unused]] std::size_t n_bytes) {
        ::operator delete[](ptr, ALIGNMENT);
    }
};