#pragma once

#include <dhb/block.h>
#include <dhb/graph.h>
#include <dhb/hash_tools.h>

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace dhb {

template <typename EdgeIt, typename GetSourceF>
std::tuple<EdgeIt, EdgeIt> thread_batch(EdgeIt batch_begin, EdgeIt batch_end,
                                        GetSourceF&& get_source_f, unsigned int thread_count,
                                        unsigned int thread_id) {
    size_t const batch_size = std::distance(batch_begin, batch_end);
    size_t const elements = batch_size / thread_count;
    if (elements == 0 || elements == batch_size) {
        bool const first_thread = thread_id == 0;
        if (first_thread) {
            return {batch_begin, batch_end};
        } else {
            return {batch_end, batch_end};
        }
    }

    size_t const position = thread_id * elements;

    EdgeIt start = std::min(batch_begin + position, batch_end);

    bool const last_thread = thread_id == (thread_count - 1);
    EdgeIt end = last_thread ? batch_end : std::min(start + elements, batch_end);

    if (start != end) {
        Vertex const predecessor =
            (start == batch_begin) ? invalidVertex() : get_source_f(*std::prev(start, 1));

        while (start != end && predecessor == get_source_f(*start)) {
            std::advance(start, 1);
        }

        if (start != end) {
            for (Vertex successor = (end == batch_end) ? invalidVertex() : get_source_f(*end);
                 end != batch_end && successor == get_source_f(*(end - 1)) &&
                 get_source_f(*end) != predecessor;
                 successor = get_source_f(*end)) {
                std::advance(end, 1);
            }
        }
    }

    return {start, end};
}

template <typename key_t> key_t key_to_thread(key_t k, uint32_t t_count) {
    if constexpr (std::is_same<key_t, uint32_t>()) {
        return reduce(hash32(k), t_count);
    } else {
        return hash64(k) % t_count;
    }
}

template <typename T> class BatchParallelizer {
  public:
    template <typename Iterator, typename GetSourceF, typename Cmp, typename F>
    void operator()(Iterator begin, Iterator end, GetSourceF&& get_source_f, Cmp cmp, F func) {
        int const t_count = omp_get_max_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            for (auto it = begin; it != end; ++it)
                func(*it);
            return;
        }

        std::sort(begin, end, cmp);
#pragma omp parallel shared(begin, end)
        {
            std::tuple<Iterator, Iterator> local_batch =
                thread_batch(begin, end, std::move(get_source_f), t_count, omp_get_thread_num());
            for (auto it = std::get<0>(local_batch); it != std::get<1>(local_batch); ++it) {
                func(*it);
            }
        }
    }

    template <typename Iterator, typename K, typename Cmp, typename F>
    void apply(Iterator begin, Iterator end, K key, Cmp cmp, F func) {
        int const t_count = omp_get_max_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            for (auto it = begin; it != end; ++it)
                func(*it);
            return;
        }

        m_batch_counts.resize((t_count + 1) * t_count);
        m_batch_slots.resize(n);

#pragma omp parallel num_threads(t_count)
        {
            auto t = omp_get_thread_num();
            assert(omp_get_num_threads() == t_count);

            auto counts_of_thread = [&](int ct) -> uint64_t* {
                return &m_batch_counts[ct * (t_count + 1)];
            };

            auto n_per_thread = n / t_count;
            auto i_begin = t * n_per_thread;
            auto i_end = i_begin + n_per_thread;
            if (t == t_count - 1)
                i_end = n;

            // First, perform a local counting sort to sort updates according to associated threads.
            auto t_counts = counts_of_thread(t);
            for (size_t at = 0; at < t_count; ++at) {
                t_counts[at] = 0;
            }

            for (size_t i = i_begin; i < i_end; ++i) {
                auto k = key(*(begin + i));
                auto at = key_to_thread(k, t_count);
                ++t_counts[at];
            }

            uint64_t psum = 0;
            for (int at = 0; at < t_count; ++at) {
                psum += t_counts[at];
                t_counts[at] = i_begin + psum;
            }
            assert(i_begin + psum == i_end);
            t_counts[t_count] = i_end;

            for (size_t irev = i_end; irev > i_begin; --irev) {
                // Iterating in reverse ensures that the sort is stable;
                // this yields a better memory access pattern when performing random access later.
                auto i = irev - 1;
                auto k = key(*(begin + i));
                auto at = key_to_thread(k, t_count);
                m_batch_slots[--t_counts[at]] = i;
            }

            // Now, let each thread collect its updates.

#pragma omp barrier

            uint64_t local_count = 0;
            for (int ot = 0; ot < t_count; ++ot) {
                auto ot_counts = counts_of_thread(ot);
                auto j_begin = ot_counts[t];
                auto j_end = ot_counts[t + 1];
                for (size_t j = j_begin; j < j_end; ++j) {
                    auto i = m_batch_slots[j];
                    auto edge = *(begin + i);
                    func(*(begin + i));
                }
                local_count += j_end - j_begin;
            }
        }
    }

    template <typename Iterator, typename K> void distribute(Iterator begin, Iterator end, K key) {
        int const t_count = omp_get_num_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            return;
        }

#pragma omp single
        {
            m_batch_counts.resize((t_count + 1) * t_count);
            m_batch_slots.resize(n);
        }

        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> uint64_t* {
            return &m_batch_counts[ct * (t_count + 1)];
        };

        auto n_per_thread = n / t_count;
        auto i_begin = t * n_per_thread;
        auto i_end = i_begin + n_per_thread;
        if (t == t_count - 1)
            i_end = n;

        // First, perform a local counting sort to sort updates according to associated threads.

        auto t_counts = counts_of_thread(t);
        for (int at = 0; at < t_count; ++at)
            t_counts[at] = 0;

        for (size_t i = i_begin; i < i_end; ++i) {
            auto k = key(*(begin + i));
            auto at = key_to_thread(k, t_count);
            ++t_counts[at];
        }

        uint64_t psum = 0;
        for (int at = 0; at < t_count; ++at) {
            psum += t_counts[at];
            t_counts[at] = i_begin + psum;
        }
        assert(i_begin + psum == i_end);
        t_counts[t_count] = i_end;

        for (size_t irev = i_end; irev > i_begin; --irev) {
            // Iterating in reverse ensures that the sort is stable;
            // this yields a better memory access pattern when performing random access later.
            auto i = irev - 1;
            auto k = key(*(begin + i));
            auto at = key_to_thread(k, t_count);
            m_batch_slots[--t_counts[at]] = i;
        }

        // Now, let each thread collect its updates.
#pragma omp barrier
    }

    template <typename Iterator, typename F> void map(Iterator begin, Iterator end, F func) {
        int const t_count = omp_get_num_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
#pragma omp master
            {
                for (auto it = begin; it != end; ++it) {
                    T elem = *it;
                    func(elem);
                }
            }
            return;
        }

        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> uint64_t* {
            return &m_batch_counts[ct * (t_count + 1)];
        };

        uint64_t local_count = 0;
        for (int ot = 0; ot < t_count; ++ot) {
            auto ot_counts = counts_of_thread(ot);
            auto j_begin = ot_counts[t];
            auto j_end = ot_counts[t + 1];
            for (size_t j = j_begin; j < j_end; ++j) {
                auto i = m_batch_slots[j];
                auto edge = *(begin + i);
                func(*(begin + i));
            }
            local_count += j_end - j_begin;
        }
    }

  private:
    std::vector<uint64_t> m_batch_counts;
    std::vector<uint64_t> m_batch_slots;
};

} // namespace dhb
