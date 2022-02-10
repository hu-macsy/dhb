#pragma once

#include <dhb/block.h>
#include <dhb/buckets.h>
#include <dhb/graph.h>

#include <cassert>
#include <memory>
#include <vector>

namespace dhb {

template <typename E> struct Vec {
  public:
    using iterator = typename BlockState<E>::iterator;
    using const_iterator = typename BlockState<E>::const_iterator;

    friend void swap(Vec& p, Vec& q) noexcept {
        using std::swap;
        swap(p.m_state, q.m_state);
        swap(p.m_handle, q.m_handle);
        swap(p.m_manager, q.m_manager);
    }

    Vec() : Vec(std::make_shared<BlockManager>(sizeof(Target))) {}

    Vec(std::shared_ptr<BlockManager> manager) : m_manager(std::move(manager)) {}

    Vec(const Vec& other) : Vec{} {
        for (auto ent : other)
            insert(ent.vertex(), ent.data());
    }

    Vec(Vec&& other) noexcept : Vec{} { swap(*this, other); }

    ~Vec() { m_manager->free_block(m_handle); }

    Vec& operator=(Vec other) {
        swap(*this, other);
        return *this;
    }

    iterator begin() { return m_state.begin(); }
    const_iterator begin() const { return m_state.begin(); }
    iterator end() { return m_state.valid_end(); }
    const_iterator end() const { return m_state.valid_end(); }

    iterator find(Vertex u) { return m_state.iterator_to(u); }
    const_iterator find(Vertex u) const { return m_state.iterator_to(u); }

    std::tuple<iterator, bool> insert(Vertex u, E ed) {
        auto& state = m_state;

        if (!state.full())
            return state.insert(u, ed);

        // We need to reallocate the adjacency list.
        auto new_bhandle = m_manager->allocate_block(state.bsize() + 1);
        BlockState<E> new_block{new_bhandle, state};
        auto result = new_block.insert(u, ed);

        auto old_block = std::move(m_state);
        auto old_bhandle = m_handle;
        m_state = std::move(new_block);
        m_handle = new_bhandle;
        m_manager->free_block(old_bhandle);

        return result;
    }

    void erase(const_iterator it) {
        auto success = m_state.remove(it->vertex());
        assert(success);
    }

    void clear() { m_state.clear(); }

  private:
    BlockState<E> m_state;
    BlockHandle m_handle;
    std::shared_ptr<BlockManager> m_manager;
};

} // namespace dhb
