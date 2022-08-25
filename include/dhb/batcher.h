#pragma once

#include <dhb/graph.h>

#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace dhb {

template <typename EdgeIt>
std::tuple<EdgeIt, EdgeIt> thread_batch(EdgeIt batch_begin, EdgeIt batch_end,
                                        unsigned int thread_count, unsigned int thread_id) {
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
            (start == batch_begin) ? invalidVertex() : std::prev(start, 1)->source;

        while (start != end && predecessor == start->source) {
            std::advance(start, 1);
        }

        if (start != end) {
            for (Vertex successor = (end == batch_end) ? invalidVertex() : end->source;
                 end != batch_end && successor == (end - 1)->source && end->source != predecessor;
                 successor = end->source) {
                std::advance(end, 1);
            }
        }
    }

    return {start, end};
}

template <typename T> class BatchParallelizer {
  public:
    template <typename Iterator, typename K, typename F>
    void operator()(Iterator begin, Iterator end, K key, F func) {
        int const t_count = omp_get_max_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            for (auto it = begin; it != end; ++it)
                func(*it);
            return;
        }

#if defined(DHB_SCATTER_SORTING)
        auto cmp = [](Edge u, Edge v) { return u.source < v.source; };
        std::sort(begin, end, cmp);
#pragma omp parallel shared(begin, end)
        {
            std::tuple<Iterator, Iterator> local_batch =
                thread_batch(begin, end, t_count, omp_get_thread_num());
            for (auto it = std::get<0>(local_batch); it != std::get<1>(local_batch); ++it) {
                func(*it);
            }
        }
#elif defined(DHB_SCATTER_DARTS)
        constexpr auto empty_cell = static_cast<unsigned int>(-1);

        // Number of slots per thread.
        // Slightly larger than the batch size per thread.
        auto s = 1 << (integer_log2_ceil((n + t_count - 1) / t_count) + 2);

        m_batch_slots.resize(s * t_count);
        std::fill(m_batch_slots.begin(), m_batch_slots.end(), empty_cell);
        m_batch_dispatched.resize(n);
        std::fill(m_batch_dispatched.begin(), m_batch_dispatched.end(), 0);

        auto slots = m_batch_slots.data();
        auto dispatched = m_batch_dispatched.data();

        std::atomic<size_t> n_done{0};

#pragma omp parallel num_threads(t_count)
        {
            auto t = omp_get_thread_num();
            assert(omp_get_num_threads() == t_count);
            auto n_per_thread = n / t_count;

            std::minstd_rand prng{m_prng_seed + t};
            std::uniform_int_distribution<int> distrib{0, s - 1};

            int r;
            for (r = 1;; ++r) {
                auto i_begin = t * n_per_thread;
                auto i_end = i_begin + n_per_thread;
                if (t == t_count - 1)
                    i_end = n;
                for (size_t i = i_begin; i < i_end; ++i) {
                    if (__atomic_load_n(&dispatched[i], __ATOMIC_RELAXED))
                        continue;
                    auto k = key(*(begin + i));
                    auto d = hash_node(k) % t_count;
                    auto j = distrib(prng);
                    __atomic_store_n(&slots[d * s + j], i, __ATOMIC_RELAXED);
                }

                size_t n_now = 0;
                for (size_t j = t * s; j < (t + 1) * s; ++j) {
                    auto i = __atomic_load_n(&slots[j], __ATOMIC_RELAXED);
                    if (i == empty_cell)
                        continue;
                    if (dispatched[i])
                        continue;
                    func(*(begin + i));
                    __atomic_store_n(&dispatched[i], 1, __ATOMIC_RELAXED);
                    __atomic_store_n(&slots[j], empty_cell, __ATOMIC_RELAXED);
                    ++n_now;
                }

                if (n_now)
                    n_done.fetch_add(n_now, std::memory_order_relaxed);
                if (n_done.load(std::memory_order_relaxed) == n)
                    break;
            }
        }

        m_prng_seed += t_count;
#elif defined(DHB_SCATTER_TWOPHASE)
        auto key_to_thread = [](unsigned int k, unsigned int t_count) -> unsigned int {
            auto hash = [](unsigned int x) -> unsigned int {
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = (x >> 16) ^ x;
                return x;
            };

            // First, hash the key to get a value that is scattered evenly in [0, 2^32).
            // For such values, the multiplication + shift yields an almost fair map,
            // see https://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/.
            return (static_cast<uint64_t>(hash(k)) * static_cast<uint64_t>(t_count)) >> 32;
        };

        m_batch_counts.resize((t_count + 1) * t_count);
        m_batch_slots.resize(n);

#pragma omp parallel num_threads(t_count)
        {
            auto t = omp_get_thread_num();
            assert(omp_get_num_threads() == t_count);

            auto counts_of_thread = [&](int ct) -> unsigned int* {
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

            unsigned int psum = 0;
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

            unsigned int local_count = 0;
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
#else
        throw std::runtime_error("DHB was compiled without support for parallel updates");
#endif
    }

    template <typename Iterator, typename K> void distribute(Iterator begin, Iterator end, K key) {
        int const t_count = omp_get_num_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count)
            return;

#if defined(DHB_SCATTER_TWOPHASE)
        auto key_to_thread = [](unsigned int k, unsigned int t_count) -> unsigned int {
            auto hash = [](unsigned int x) -> unsigned int {
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = (x >> 16) ^ x;
                return x;
            };

            // First, hash the key to get a value that is scattered evenly in [0, 2^32).
            // For such values, the multiplication + shift yields an almost fair map,
            // see https://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/.
            return (static_cast<uint64_t>(hash(k)) * static_cast<uint64_t>(t_count)) >> 32;
        };

#pragma omp single
        {
            m_batch_counts.resize((t_count + 1) * t_count);
            m_batch_slots.resize(n);
        }

        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> unsigned int* {
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

        unsigned int psum = 0;
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

#elif defined(DHB_SCATTER_COUNTING)
        auto key_to_thread = [](unsigned int k, unsigned int t_count) -> unsigned int {
            // More on this xor shift hash function can be found at:
            // https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
            auto hash = [](uint32_t x) -> uint32_t {
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = ((x >> 16) ^ x) * 0x45d9f3b;
                x = (x >> 16) ^ x;
                return x;
            };

            // We are using a fast alternative to the modulo reduction from
            // Daniel Lemire's blog. See:
            // https://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction First,
            // hash the key to get a value that is scattered evenly in [0, 2^32). For such values,
            // the multiplication + shift yields an almost fair map.
            return (static_cast<uint64_t>(hash(k)) * static_cast<uint64_t>(t_count)) >> 32;
        };

#pragma omp single
        // We start with an input buffer which is the batch that is to be
        // operated on (distribute workload). In order to operate in fair balance
        // we need to distribute the edges of the batch to our thread pool.
        // Therefore, each thread first identifies the edges that have to be
        // moved around using a hash function for the edges key (e.g. the edge
        // source). Second, an additional out buffer is allocated, offsets and
        // edges to be moved computed. Finally, all edges are moved to that out
        // buffer then ready to be used by the map() function.

        {
            m_batch_counts.resize((t_count + 1) * t_count);
            m_out.resize(t_count);
            m_wp.resize(t_count);
        }

        int const t_id = omp_get_thread_num();

        auto counts_of_thread = [&](int t_id) -> unsigned int* {
            return &m_batch_counts[t_id * (t_count + 1)];
        };

        auto n_per_thread = n / t_count;
        ptrdiff_t i_begin = t_id * n_per_thread;
        ptrdiff_t i_end = i_begin + n_per_thread;
        if (t_id == t_count - 1) {
            i_end = n;
        }

        // Determine send counts.
        auto t_counts = counts_of_thread(t_id);
        std::fill_n(std::begin(t_counts), t_count, 0);

        for (ptrdiff_t i = i_begin; i < i_end; ++i) {
            auto it = begin + i;
            auto k = key(*it);
            auto at = key_to_thread(k, t_count);
            ++t_counts[at];
        }

#pragma omp barrier
        // Do a prefix sum over the send count *to* the current thread.

        unsigned int psum = 0;
        for (int rt = 0; rt < t_count; ++rt) {
            auto rt_counts = counts_of_thread(rt);
            auto c = rt_counts[t_id];
            rt_counts[t_id] = psum;
            psum += c;
        }

        m_out[t_id].resize(psum);

#pragma omp barrier
        // We have determined which edges have to be moved to which thread buffers.
        // Now it is time to move the edges to the thread buffers.

        m_wp[t_id].resize(t_count);
        for (int at = 0; at < t_count; ++at)
            m_wp[t_id][at] = m_out[at].data() + t_counts[at];
        auto wp_ptr = m_wp[t_id].data();

        for (ptrdiff_t i = i_begin; i < i_end; ++i) {
            auto it = begin + i;
            auto k = key(*it);
            auto at = key_to_thread(k, t_count);
            *(wp_ptr[at]++) = *it;
        }

#pragma omp barrier
#else
        throw std::runtime_error("DHB was compiled without support for parallel updates");
#endif
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
#if defined(DHB_SCATTER_TWOPHASE)
        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> unsigned int* {
            return &m_batch_counts[ct * (t_count + 1)];
        };

        unsigned int local_count = 0;
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
#elif defined(DHB_SCATTER_COUNTING)
        auto t = omp_get_thread_num();

        for (auto& elem : m_out[t])
            func(elem);
#else
        throw std::runtime_error("DHB was compiled without support for parallel updates");
#endif
    }

  private:
#if defined(DHB_SCATTER_DARTS)
    unsigned int m_prng_seed = 0;
    std::vector<unsigned int> m_batch_slots;
    std::vector<char> m_batch_dispatched;
#elif defined(DHB_SCATTER_TWOPHASE)
    std::vector<unsigned int> m_batch_counts;
    std::vector<unsigned int> m_batch_slots;
#elif defined(DHB_SCATTER_COUNTING)
    std::vector<unsigned int> m_batch_counts;
    std::vector<std::vector<T>> m_out;
    std::vector<std::vector<T*>> m_wp;
#endif
};

} // namespace dhb
