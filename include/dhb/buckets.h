#pragma once

#include <dhb/block.h>

#include <array>
#include <list>
#include <set>
#include <unordered_set>
#include <vector>

namespace dhb {

#ifndef DHB_SYSTEM_ALLOCATOR

struct BlockBucket {
    struct AgeLess {
        bool operator()(BlockArray* lhs, BlockArray* rhs) const { return lhs->age() < rhs->age(); }
    };

    std::mutex mutex;
    unsigned int next_age = 1;
    std::unordered_set<BlockArray*> all;
    std::set<BlockArray*, AgeLess> available;
};

struct BlockCache {
    std::array<BlockBucket, 48> buckets;
};

#endif // DHB_SYSTEM_ALLOCATOR

class BlockManager {
  public:
    BlockManager(unsigned int bytes_per_entry);

    BlockManager(const BlockManager&) = delete;

    BlockManager(BlockManager&& other) = delete;

    ~BlockManager();

    BlockManager& operator=(const BlockManager& other) = delete;

    BlockHandle allocate_block(Degree degree);

    void free_block(BlockHandle bptr);

  private:
    unsigned int m_bytes_per_entry;

#ifndef DHB_SYSTEM_ALLOCATOR
    BlockCache* get_own_cache();

    // Stores the block cache of each thread.
    pthread_key_t m_cache;

    // Protectes access to m_all_caches.
    std::mutex m_meta_mutex;

    // Stores all thread caches.
    // TODO: Add some mechanism (std::shared_ptr?) to deal with exiting threads.
    // We need caches with stable addresses; simply use a std::list here.
    std::list<BlockCache> m_all_caches;
#endif

    size_t m_min_bsize{4};
};

} // namespace dhb
