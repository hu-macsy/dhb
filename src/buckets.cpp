#include <dhb/buckets.h>
#include <dhb/integer_log2.h>

#include <cassert>

namespace dhb {

dhb::BlockManager::BlockManager(unsigned int bytes_per_entry) : m_bytes_per_entry(bytes_per_entry) {
#ifndef DHB_SYSTEM_ALLOCATOR
    pthread_key_create(&m_cache, nullptr);
#endif
}

BlockManager::~BlockManager() {
#ifndef DHB_SYSTEM_ALLOCATOR
    while (!m_all_caches.empty()) {
        auto cache = &m_all_caches.front();
        for (size_t i = 0; i < cache->buckets.size(); ++i) {
            auto bucket = &cache->buckets[i];
            for (BlockArray* ba : bucket->all)
                delete ba;
        }
        m_all_caches.pop_front();
    }

    pthread_key_delete(m_cache);
#endif
}

BlockHandle BlockManager::allocate_block(Degree degree) {
    // the max() operation will make sure that we return at least BA of bsize 4 / index 2,
    // no matter the degree of the vertex
    size_t const index = std::max(2u, integer_log2_ceil(degree));
    size_t const bsize_requested = 1 << index;
    assert(bsize_requested >= m_min_bsize);

#ifndef DHB_SYSTEM_ALLOCATOR
    auto cache = get_own_cache();
    if (index >= cache->buckets.size())
        throw std::runtime_error("bucket size is too large");
    auto bucket = &cache->buckets[index];

    std::unique_lock<std::mutex> lock{bucket->mutex};
    [&] {
        // Fast path if a BA is available.
        unsigned int age;
        if (!bucket->available.empty())
            return;
        age = bucket->next_age++;

        // Allocate the BA with locks dropped.
        lock.unlock();
        auto ba = new BlockArray(m_bytes_per_entry, bsize_requested, cache, age);
        lock.lock();

        bucket->all.insert(ba);
        bucket->available.insert(ba);
    }();
    assert(!bucket->available.empty());

    auto ba = *bucket->available.begin();
    assert(!ba->occupied());

    auto bhandle = ba->take();

    if (ba->occupied()) {
        if (ba->mostly_free()) {
            ba->recycle();
        } else {
            bucket->available.erase(bucket->available.begin());
        }
    }

    return bhandle;
#else
    void* ptr_entries;
    void* ptr_htab = nullptr;
    ptr_entries = operator new(m_bytes_per_entry* bsize_requested);
    if (uses_htab(bsize_requested))
        ptr_htab = operator new(sizeof(index_type) * 2 * bsize_requested);
    return BlockHandle(bsize_requested, ptr_entries, reinterpret_cast<index_type*>(ptr_htab));
#endif
}

void BlockManager::free_block(BlockHandle bhandle) {
#ifndef DHB_SYSTEM_ALLOCATOR
    if (!bhandle)
        return;

    auto ba = bhandle.block_array();
    auto index = integer_log2_ceil(ba->bsize());
    auto cache = ba->m_cache;
    auto bucket = &cache->buckets[index];

    std::unique_lock<std::mutex> lock{bucket->mutex};

    ba->reclaim(bhandle);

    if (ba->occupied()) {
        if (ba->mostly_free()) {
            ba->recycle();
            bucket->available.insert(ba);
        }
    }
#else
    operator delete(bhandle.access_entries(), m_bytes_per_entry* bhandle.bsize());
    operator delete(bhandle.access_htab(), sizeof(index_type) * 2 * bhandle.bsize());
#endif
}

#ifndef DHB_SYSTEM_ALLOCATOR
BlockCache* BlockManager::get_own_cache() {
    void* p = pthread_getspecific(m_cache);
    if (!p) {
        std::unique_lock<std::mutex> meta_lock{m_meta_mutex};

        m_all_caches.emplace_back();
        p = &m_all_caches.back();
        if (pthread_setspecific(m_cache, p))
            std::terminate();
    }
    return static_cast<BlockCache*>(p);
}
#endif

} // namespace dhb
