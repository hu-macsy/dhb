#pragma once

#include <dhb/graph.h>
#include <dhb/integer_log2.h>

#include <algorithm>
#include <cassert>
#include <mutex>
#include <sys/mman.h>
#include <vector>

namespace dhb {

class BlockHandle;
class BlockArray;

using index_type = size_t;
inline index_type illegalIndex() { return static_cast<index_type>(-1); };
inline index_type tombstoneIndex() { return static_cast<index_type>(-2); };

// TODO: This could take the entry size into account.
inline bool uses_htab(size_t bsize) {
    const size_t cache_line = 64;
    // Each entry is typically 16 bytes or so.
    // Preliminary experiments how that using the hash index is only useful if the block becomes
    // larger than several cache lines (approx. 16 or so).
    // Hence, the following is a reasonably good heuristic:
    return 16 * bsize > 16 * cache_line;
}

// Taken from https://stackoverflow.com/a/12996028.
inline unsigned int hash_node(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

#ifndef DHB_SYSTEM_ALLOCATOR
class BlockHandle {
  public:
    BlockHandle() : m_ba{nullptr}, m_index{0} {}

    BlockHandle(BlockArray* ba, unsigned int index) : m_ba{ba}, m_index{index} {}

    explicit operator bool() { return m_ba; }

    unsigned int bsize();
    void* access_entries();
    index_type* access_htab();

    BlockArray* block_array() { return m_ba; }
    unsigned int index() { return m_index; }

  private:
    BlockArray* m_ba;
    unsigned int m_index;
};
#else
class BlockHandle {
  public:
    BlockHandle() = default;

    BlockHandle(unsigned int bsize, void* entries, index_type* htab)
        : m_bsize{bsize}, m_entries{entries}, m_htab{htab} {}

    explicit operator bool() { return m_entries; }

    unsigned int bsize() { return m_bsize; }
    void* access_entries() { return m_entries; }
    index_type* access_htab() { return m_htab; }

  private:
    unsigned int m_bsize = 0;
    void* m_entries = nullptr;
    index_type* m_htab = nullptr;
};
#endif

#ifndef DHB_SYSTEM_ALLOCATOR

struct BlockCache;

class BlockArray {
  public:
    BlockArray(unsigned int bytes_per_entry, unsigned int bsize, BlockCache* cache,
               unsigned int age)
        : m_cache(cache), m_age{age}, m_bytes_per_entry(bytes_per_entry), m_bsize(bsize) {
        const size_t min_capacity = 0x200000; // 2 MiB.
        const size_t min_blocks = 10;
        size_t allocation_per_block;
        if (uses_htab(m_bsize)) {
            allocation_per_block = block_bytes() + bsize * 2 * sizeof(index_type);
        } else {
            allocation_per_block = block_bytes();
        }
        // Determine the capacity such that:
        // * we allocate a minimum number of bytes,
        // * we can hold a minimum number of blocks.
        auto order = std::max(integer_log2_ceil(min_capacity),
                              integer_log2_ceil(min_blocks * allocation_per_block));
        m_capacity = size_t(1) << order;
        m_nblocks = m_capacity / allocation_per_block;
        assert(m_nblocks);

        // Perform the actual allocation.
        void* vma =
            mmap(nullptr, m_capacity, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (vma == MAP_FAILED)
            throw std::bad_alloc{};

        // Place the hash tables at the end of the memory range (this also ensures alignment).
        auto va = reinterpret_cast<uintptr_t>(vma);
        m_data = vma;
        if (uses_htab(m_bsize)) {
            auto bytes_of_all_htabs = m_nblocks * bsize * 2 * sizeof(index_type);
            m_htab = reinterpret_cast<index_type*>(va + (m_capacity - bytes_of_all_htabs));
        }

        // Mark all blocks as available.
        m_alloc_stack.reserve(m_nblocks);
        for (size_t i = 0; i < m_nblocks; ++i)
            m_alloc_stack.push_back(m_nblocks - 1 - i);
    }

    BlockArray(const BlockArray&) = delete;

    ~BlockArray() {
        if (munmap(m_data, m_capacity))
            std::terminate();
    }

    BlockArray& operator=(const BlockArray&) = delete;

    unsigned int age() { return m_age; }

    unsigned int bytes_per_entry() const { return m_bytes_per_entry; }
    unsigned int bsize() const { return m_bsize; }

    size_t number_of_blocks() const { return m_nblocks; }
    size_t block_bytes() const { return static_cast<size_t>(m_bytes_per_entry) * m_bsize; }

    // returns true if at least one of all associated blocks are empty
    bool occupied() const;

    bool mostly_free() const { return m_free_stack.size() > m_nblocks / 2; }

    // takes a free block from the block array, or returns nullptr.
    BlockHandle take();

    void reclaim(BlockHandle);

    void recycle();

    void* access_storage_of(unsigned int index) {
        return static_cast<char*>(m_data) + index * block_bytes();
    }

    index_type* access_htab_of(unsigned int index) {
        if (uses_htab(m_bsize))
            return m_htab + index * 2 * m_bsize;
        return nullptr;
    }

    // TODO: the following is quite hacky:
  public:
    BlockCache* m_cache;

  private:
    unsigned int m_age;
    unsigned int const m_bytes_per_entry;
    unsigned int const m_bsize;
    void* m_data;
    index_type* m_htab;
    std::vector<unsigned int> m_alloc_stack;
    std::vector<unsigned int> m_free_stack;
    size_t m_capacity;
    size_t m_nblocks;
};

inline unsigned int BlockHandle::bsize() { return m_ba->bsize(); }
inline void* BlockHandle::access_entries() { return m_ba->access_storage_of(m_index); }
inline index_type* BlockHandle::access_htab() { return m_ba->access_htab_of(m_index); }
#endif

//---------------------------------------------------------------------------------------
// BlockState: mutation and lookup functions.
//---------------------------------------------------------------------------------------

inline index_type index_hash_table(unsigned int hash, index_type current_index,
                                   unsigned int bsize) {
    return (hash + current_index) & (2 * bsize - 1);
}

template <typename E> class BlockState {
  public:
    struct Entry {
        Vertex vertex;
        E data;
    };

    static constexpr size_t entry_size = sizeof(Entry);

    struct proxy {
        explicit proxy(Entry* ptr) : m_ptr(ptr) {}

        Vertex vertex() const { return m_ptr->vertex; }
        E& data() const { return m_ptr->data; }

      private:
        Entry* m_ptr;
    };

    struct const_proxy {
        explicit const_proxy(const Entry* ptr) : m_ptr(ptr) {}

        Vertex vertex() const { return m_ptr->vertex; }
        const E& data() const { return m_ptr->data; }

      private:
        const Entry* m_ptr;
    };

    // Helper class to implement operator-> of iterator.
    struct arrow {
        arrow(Entry* ptr) : m_p(ptr) {}

        proxy* operator->() { return &m_p; }

      private:
        proxy m_p;
    };

    // Helper class to implement operator-> of const_iterator.
    struct const_arrow {
        const_arrow(const Entry* ptr) : m_p(ptr) {}

        const_proxy* operator->() { return &m_p; }

      private:
        const_proxy m_p;
    };

    struct iterator {
        using value_type = proxy;
        using reference = proxy; // decltype(operator*)
        using pointer = arrow;   // decltype(operator->)
        using difference_type = ptrdiff_t;
        // In the legacy sense, this is only an input iterator since the reference type
        // is not a true reference. In the modern sense, this is a random access iterator.
        using iterator_category = std::input_iterator_tag;
        using iterator_concept = std::random_access_iterator_tag;

        friend bool operator==(iterator lhs, iterator rhs) { return lhs.m_ptr == rhs.m_ptr; }
        friend bool operator!=(iterator lhs, iterator rhs) { return lhs.m_ptr != rhs.m_ptr; }

        friend ptrdiff_t operator-(iterator lhs, iterator rhs) { return lhs.m_ptr - rhs.m_ptr; }

        iterator() : m_ptr(nullptr) {}

        explicit iterator(Entry* ptr) : m_ptr(ptr) {}

        proxy operator*() const { return proxy(m_ptr); }
        arrow operator->() const { return arrow(m_ptr); }

        iterator operator+(ptrdiff_t delta) const {
            iterator copy = *this;
            copy += delta;
            return copy;
        }

        iterator& operator++() {
            ++m_ptr;
            return *this;
        }
        iterator& operator--() {
            --m_ptr;
            return *this;
        }
        iterator operator++(int) {
            iterator copy = *this;
            ++(*this);
            return copy;
        }
        iterator operator--(int) {
            iterator copy = *this;
            --(*this);
            return copy;
        }
        iterator& operator+=(ptrdiff_t delta) {
            m_ptr += delta;
            return *this;
        }
        iterator& operator-=(ptrdiff_t delta) {
            m_ptr -= delta;
            return *this;
        }

        Entry* get_ptr() const { return m_ptr; }

      private:
        Entry* m_ptr;
    };

    struct const_iterator {
        using value_type = const_proxy;
        using reference = const_proxy; // decltype(operator*)
        using pointer = const_arrow;   // decltype(operator->)
        using difference_type = ptrdiff_t;
        // In the legacy sense, this is only an input iterator since the reference type
        // is not a true reference. In the modern sense, this is a random access iterator.
        using iterator_category = std::input_iterator_tag;
        using iterator_concept = std::random_access_iterator_tag;

        friend bool operator==(const_iterator lhs, const_iterator rhs) {
            return lhs.m_ptr == rhs.m_ptr;
        }
        friend bool operator!=(const_iterator lhs, const_iterator rhs) {
            return lhs.m_ptr != rhs.m_ptr;
        }

        friend ptrdiff_t operator-(const_iterator lhs, const_iterator rhs) {
            return lhs.m_ptr - rhs.m_ptr;
        }

        const_iterator() : m_ptr(nullptr) {}

        explicit const_iterator(const Entry* ptr) : m_ptr(ptr) {}

        // Conversion of iterator to const_iterator.
        const_iterator(iterator other) : m_ptr(other.get_ptr()) {}

        const_proxy operator*() const { return const_proxy(m_ptr); }
        const_arrow operator->() const { return const_arrow(m_ptr); }

        const_iterator operator+(ptrdiff_t delta) const {
            const_iterator copy = *this;
            copy += delta;
            return copy;
        }

        const_iterator& operator++() {
            ++m_ptr;
            return *this;
        }
        const_iterator& operator--() {
            --m_ptr;
            return *this;
        }
        const_iterator operator++(int) {
            const_iterator copy = *this;
            ++(*this);
            return copy;
        }
        const_iterator operator--(int) {
            const_iterator copy = *this;
            --(*this);
            return copy;
        }
        const_iterator& operator+=(ptrdiff_t delta) {
            m_ptr += delta;
            return *this;
        }
        const_iterator& operator-=(ptrdiff_t delta) {
            m_ptr -= delta;
            return *this;
        }

        const Entry* get_ptr() const { return m_ptr; }

      private:
        const Entry* m_ptr;
    };

    BlockState() : m_bsize{0}, m_degree{0}, m_entries{nullptr}, m_htab{nullptr} {}

    BlockState(BlockHandle bhandle)
        : m_bsize{bhandle.bsize()}, m_degree{0},
          m_entries{static_cast<Entry*>(bhandle.access_entries())}, m_htab{bhandle.access_htab()} {
        // Clear the hash table.
        if (uses_htab(m_bsize)) {
            index_type* const hash_table = m_htab;
            for (size_t j = 0; j < 2 * m_bsize; ++j) {
                hash_table[j] = illegalIndex();
            }
        }
    }

    BlockState(BlockHandle bhandle, const BlockState& other)
        : m_bsize{bhandle.bsize()}, m_degree{other.m_degree},
          m_entries{static_cast<Entry*>(bhandle.access_entries())}, m_htab{bhandle.access_htab()} {
        assert(other.m_degree <= m_bsize);
        for (size_t i = 0; i < other.m_degree; ++i)
            m_entries[i] = other.m_entries[i];
        if (uses_htab(m_bsize))
            rebuild_ht();
    }

    BlockState(const BlockState&) = delete;

    BlockState(BlockState&& other) = default;

    BlockState& operator=(const BlockState&) = delete;

    BlockState& operator=(BlockState&&) = default;

    ptrdiff_t bsize() const { return m_bsize; }

    ptrdiff_t degree() const { return m_degree; }

    // Returns true if the neighbor was inserted.
    // Returns false if the neighbor was found.
    std::tuple<iterator, bool> insert(Vertex v, E ed) {
        assert(m_degree < m_bsize);

        if (uses_htab(m_bsize)) {
            index_type* htab = m_htab;
            auto h = hash_node(v);
            index_type j;
            index_type ts = illegalIndex();
            for (index_type i = 0; true; ++i) {
                if (i == m_bsize) {
                    assert(ts != illegalIndex());
                    break;
                }
                j = index_hash_table(h, i, m_bsize);

                if (htab[j] == illegalIndex())
                    break;
                if (htab[j] == tombstoneIndex()) {
                    ts = j;
                    continue;
                }
                auto ptr = &m_entries[htab[j]];
                if (ptr->vertex == v)
                    return {iterator(ptr), false};
            }

            // If we did hit a tombstone, insert at the tombstone.
            if (ts != illegalIndex())
                j = ts;

            // Insert into the adjacency list.
            auto i = m_degree++;
            htab[j] = i;
            m_entries[i] = {v, ed};
            return {iterator(&m_entries[i]), true};
        } else {
            // Find the node with the adjacency list.
            for (index_type i = 0; i < m_degree; ++i)
                if (m_entries[i].vertex == v)
                    return {iterator(&m_entries[i]), false};

            // Insert into the adjacency list.
            auto i = m_degree++;
            m_entries[i] = {v, ed};
            return {iterator(&m_entries[i]), true};
        }
    }

    // true = node was removed.
    bool remove(Vertex const v) {
        if (uses_htab(m_bsize)) {
            // Find the node within the hash table.
            index_type* const hash_table = m_htab;

            unsigned int const hv = hash_node(v);
            size_t jv = 0; // Hash table index.
            for (index_type i = 0; true; ++i) {
                if (i == m_bsize) {
                    return false;
                }

                jv = index_hash_table(hv, i, m_bsize);

                if (hash_table[jv] == illegalIndex()) {
                    return false;
                }

                if (hash_table[jv] == tombstoneIndex()) {
                    continue;
                }

                if ((m_entries + hash_table[jv])->vertex == v) {
                    break;
                }
            }

            auto back_it = std::prev(m_entries + m_degree);
            if (m_entries + hash_table[jv] != back_it) {
                auto w = back_it->vertex;

                unsigned int const hw = hash_node(w);
                size_t jw = 0; // Hash table index.
                for (index_type i = 0; true; ++i) {
                    // Otherwise, w is not in the HT.
                    assert(i != m_bsize);

                    jw = index_hash_table(hw, i, m_bsize);

                    // Otherwise, w is not in the HT.
                    assert(hash_table[jw] != illegalIndex());

                    if (hash_table[jw] == tombstoneIndex()) {
                        continue;
                    }

                    if ((m_entries + hash_table[jw])->vertex == w) {
                        break;
                    }
                }

                *(m_entries + hash_table[jv]) = *back_it;
                hash_table[jw] = hash_table[jv];
            }

            hash_table[jv] = tombstoneIndex();
            --m_degree;

            return true;
        } else {
            // Find the index.
            index_type iv;
            for (iv = 0; iv < m_degree; ++iv)
                if (m_entries[iv].vertex == v)
                    break;

            // Otherwise, v is not in the block.
            if (iv == m_degree)
                return false;

            // Swap with the last entry.
            if (iv + 1 != m_degree)
                m_entries[iv] = m_entries[m_degree - 1];

            // Remove the last entry.
            --m_degree;

            return true;
        }
    }

    void clear() {
        m_degree = 0;
        if (uses_htab(m_bsize)) {
            for (size_t j = 0; j < 2 * m_bsize; ++j)
                m_htab[j] = illegalIndex();
        }
    }

    template <typename L> void sort(L less) {
        std::sort(m_entries, m_entries + m_degree, std::move(less));
        if (uses_htab(m_bsize))
            rebuild_ht();
    }

    iterator iterator_to(Vertex v) const {
        if (uses_htab(m_bsize)) {
            // Find the node within the hash table.
            auto h = hash_node(v);
            size_t j = 0; // Hash table index.
            index_type* const hash_table = m_htab;
            for (index_type i = 0; true; ++i) {
                if (i == m_bsize)
                    return valid_end();
                j = index_hash_table(h, i, m_bsize);

                if (hash_table[j] == illegalIndex())
                    return valid_end();
                if (hash_table[j] == tombstoneIndex())
                    continue;
                if ((m_entries + hash_table[j])->vertex == v)
                    break;
            }

            return iterator(m_entries + hash_table[j]);
        } else {
            // Find the node within the adjacency list.
            for (index_type i = 0; i < m_degree; ++i)
                if (m_entries[i].vertex == v)
                    return iterator(&m_entries[i]);

            return valid_end();
        }
    }

    bool empty() const { return !m_degree; }

    bool full() const { return m_degree == m_bsize; }

    iterator valid_end() const { return iterator(m_entries + m_degree); }

    iterator begin() { return iterator(m_entries); }
    const_iterator begin() const { return const_iterator(m_entries); }

    const_iterator cbegin() const { return const_iterator(m_entries); }

    const_iterator end() const { return const_iterator(m_entries + m_bsize); }

  private:
    void rebuild_ht() {
        assert(uses_htab(m_bsize));

        index_type* const hash_table = m_htab;

        // Clear the hash table.
        for (size_t j = 0; j < 2 * m_bsize; ++j) {
            hash_table[j] = illegalIndex();
        }

        // Re-insert all vertices.
        for (size_t i = 0; i < m_degree; ++i) {
            auto v = m_entries[i].vertex;

            // Find a free slot within the hash table.
            unsigned int const h = hash_node(v);
            size_t j = 0; // Hash table index.
            for (index_type p = 0; true; ++p) {
                assert(p < m_bsize);
                j = index_hash_table(h, p, m_bsize);

                if (hash_table[j] == illegalIndex() || hash_table[j] == tombstoneIndex()) {
                    break;
                }

                // There cannot be any duplicates.
                assert((m_entries + hash_table[j])->vertex != v);
            }

            // Insert into the hash table.
            hash_table[j] = i;
        }
    }

    unsigned int m_bsize;
    unsigned int m_degree;
    Entry* m_entries;
    index_type* m_htab;
};

template <typename E> constexpr size_t bytes_per_entry_for() { return BlockState<E>::entry_size; }

//---------------------------------------------------------------------------------------
// BlockArray: allocation functions.
//---------------------------------------------------------------------------------------

#ifndef DHB_SYSTEM_ALLOCATOR
inline bool BlockArray::occupied() const { return m_alloc_stack.empty(); }

inline BlockHandle BlockArray::take() {
    assert(!m_alloc_stack.empty());

    auto index = m_alloc_stack.back();
    m_alloc_stack.pop_back();
    return {this, index};
}

inline void BlockArray::reclaim(BlockHandle bhandle) {
    assert(bhandle.block_array() == this);

    m_free_stack.push_back(bhandle.index());
}

inline void BlockArray::recycle() {
    for (unsigned int index : m_free_stack)
        m_alloc_stack.push_back(index);
    m_free_stack.clear();
}
#endif

} // namespace dhb
