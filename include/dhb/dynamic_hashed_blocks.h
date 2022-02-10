#pragma once

#include <dhb/batcher.h>
#include <dhb/block.h>
#include <dhb/buckets.h>
#include <dhb/graph.h>

#include <cassert>
#include <functional>
#include <memory>
#include <mutex>
#include <omp.h>
#include <vector>

namespace dhb {

std::vector<Degree> degrees_from(Edges const&);

template <typename E> struct Matrix {
  public:
    struct NeighborView {
        using iterator = typename BlockState<E>::iterator;
        using proxy = typename BlockState<E>::proxy;

        NeighborView(Matrix* g, Vertex u) : m_g{g}, m_u{u} {}

        auto begin() { return m_g->m_vertices[m_u].begin(); }

        auto end() { return m_g->m_vertices[m_u].valid_end(); }

        auto iterator_to(Vertex v) { return m_g->m_vertices[m_u].iterator_to(v); }

        bool exists(Vertex v) { return iterator_to(v) != end(); }

        void clear() { m_g->m_vertices[m_u].clear(); }

        size_t degree() { return end() - begin(); }

        proxy operator[](size_t i) { return *(begin() + i); }

        std::tuple<iterator, bool> insert(Vertex v, E ed) {
            auto& state = m_g->m_vertices[m_u];

            if (!state.full())
                return state.insert(v, ed);

            // We need to reallocate the adjacency block.
            auto new_bhandle = m_g->m_manager->allocate_block(state.bsize() + 1);
            BlockState<E> new_block{new_bhandle, state};
            auto result = new_block.insert(v, ed);

            auto old_block = std::move(m_g->m_vertices[m_u]);
            auto old_bhandle = m_g->m_handles[m_u];
            m_g->m_vertices[m_u] = std::move(new_block);
            m_g->m_handles[m_u] = new_bhandle;
            m_g->m_manager->free_block(old_bhandle);

            return result;
        }

      private:
        Matrix* m_g;
        Vertex m_u;
    };

    struct ConstNeighborView {
        using iterator = typename BlockState<E>::const_iterator;
        using proxy = typename BlockState<E>::const_proxy;

        ConstNeighborView(const Matrix* g, Vertex u) : m_g{g}, m_u{u} {}

        auto begin() { return m_g->m_vertices[m_u].begin(); }

        auto end() { return m_g->m_vertices[m_u].valid_end(); }

        auto iterator_to(Vertex v) { return m_g->m_vertices[m_u].iterator_to(v); }

        bool exists(Vertex v) { return iterator_to(v) != end(); }

        size_t degree() { return end() - begin(); }

        proxy operator[](size_t i) { return *(begin() + i); }

      private:
        const Matrix* m_g;
        Vertex m_u;
    };

    friend void swap(Matrix& p, Matrix& q) noexcept {
        using std::swap;
        swap(p.m_vertices, q.m_vertices);
        swap(p.m_handles, q.m_handles);
        swap(p.m_manager, q.m_manager);
    }

    Matrix() : Matrix(std::make_shared<BlockManager>(bytes_per_entry_for<E>())) {}

    Matrix(std::shared_ptr<BlockManager> manager) : m_manager(std::move(manager)) {}

    Matrix(Edges&& edges) : Matrix{} {
        std::vector<Degree> const extracted_degrees = degrees_from(edges);
        Degree const max_degree =
            *std::max_element(std::cbegin(extracted_degrees), std::end(extracted_degrees));

        unsigned int count = graph::vertex_count(edges);
        m_vertices.reserve(count);
        m_handles.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            auto degree = extracted_degrees[i];
            dhb::BlockHandle bhandle = m_manager->allocate_block(degree);
            BlockState<E> block{bhandle};
            m_vertices.push_back(std::move(block));
            m_handles.push_back(bhandle);
        }

        // insert vertex and it's neighbours to respective block arrays and increase size
        // of block array if necessary
        for (Edges::iterator edge = std::begin(edges); edge != std::end(edges);) {
            Vertex source = edge->source;
            Vertex next_source = source;
            Degree degree = extracted_degrees[source];

            auto& block = m_vertices[source];

            while (source == next_source && edge != std::end(edges)) {
                block.insert(edge->target.vertex, edge->target.data);
                ++edge;
                if (edge != std::end(edges)) {
                    next_source = edge->source;
                }
            }
        }
    }

    Matrix(Vertex vertex_count) : Matrix() { resize(vertex_count); }

    Matrix(const Matrix& other) : Matrix{} {
        if (other.vertices_count())
            resize(other.vertices_count());
        other.for_edges([&](Vertex u, Vertex v, E ed) { neighbors(u).insert(v, ed); });
    }

    Matrix(Matrix&& other) noexcept : Matrix{} { swap(*this, other); }

    ~Matrix() {
        for (auto bhandle : m_handles)
            m_manager->free_block(bhandle);
    }

    Matrix& operator=(Matrix other) {
        swap(*this, other);
        return *this;
    }

    void resize(Vertex id) {
        m_vertices.resize(id);
        m_handles.resize(id);
    }

    void preallocate(Vertex u, unsigned int degree) {
        auto& associated_block = m_vertices[u];

        if (associated_block.bsize() < degree) {
            auto new_bhandle = m_manager->allocate_block(degree);
            BlockState<E> new_block{new_bhandle, associated_block};

            auto old_block = std::move(m_vertices[u]);
            auto old_bhandle = m_handles[u];
            m_vertices[u] = std::move(new_block);
            m_handles[u] = new_bhandle;
            m_manager->free_block(old_bhandle);
        }
    }

    bool insert(Vertex u, Vertex v, E ed, bool do_update = true) {
        auto result = neighbors(u).insert(v, ed);
        if (std::get<1>(result))
            return true;
        if (do_update) {
            std::get<0>(result)->data() = ed;
            return true;
        }
        return false;
    }

    template <typename EdgeIt> void insert(EdgeIt begin, EdgeIt end, bool do_update = true) {
        for (Edges::const_iterator edge = begin; edge != end; ++edge) {
            insert(edge->source, edge->target.vertex, edge->target.data, do_update);
        }
    }

    bool removeEdge(Vertex source, Vertex target) {
        if (source >= m_vertices.size()) {
            return false;
        }

        auto& associated_block = m_vertices[source];
        if (associated_block.remove(target)) {
            return true;
        }

        return false;
    }

    template <typename L> void sort(Vertex u, L less) { m_vertices[u].sort(std::move(less)); }

    NeighborView neighbors(Vertex u) { return {this, u}; }

    ConstNeighborView neighbors(Vertex u) const { return {this, u}; }

    template <typename F> void for_nodes(F&& handle) const {
        for (Vertex node = 0; node < m_vertices.size(); ++node) {
            handle(node);
        }
    }

    template <typename F> void for_edges(F&& handle) const {
        for (Vertex i = 0; i < m_vertices.size(); ++i) {
            if (m_vertices[i].degree() > 0) {
                auto& block = m_vertices[i];
                for (auto it = block.cbegin(); it != block.valid_end(); it++) {
                    handle(i, it->vertex(), it->data());
                }
            }
        }
    }

    template <typename F> void for_neighbors_node(Vertex const node, F&& handle) const {
        assert(node < m_vertices.size());
        auto& block = m_vertices[node];

        for (auto it = block.cbegin(); it != block.valid_end(); it++) {
            handle(it->vertex());
        }
    }

    template <typename F> void for_neighbors_target(Vertex const node, F&& handle) const {
        assert(node < m_vertices.size());
        auto& block = m_vertices[node];

        for (auto it = block.cbegin(); it != block.valid_end(); it++) {
            handle(it->vertex(), it->data());
        }
    }

    Vertex degree(Vertex node) const {
        assert(node < m_vertices.size());
        return m_vertices[node].degree();
    }

    Vertex vertices_count() const { return m_vertices.size(); }

    Vertex edges_count() const {
        return std::accumulate(
            std::cbegin(m_vertices), std::end(m_vertices), 0u,
            [](Vertex sum, const BlockState<E>& block) { return sum + block.degree(); });
    }

  private:
    std::vector<BlockState<E>> m_vertices;
    std::vector<BlockHandle> m_handles;
    std::shared_ptr<BlockManager> m_manager;
};

template <typename E> inline std::tuple<Degree, Vertex> max_degree(const Matrix<E>& dhb) {
    Vertex v = 0;
    for (Vertex u = 1; u < dhb.vertices_count(); ++u) {
        if (dhb.degree(u) > dhb.degree(v))
            v = u;
    }
    assert(v < dhb.vertices_count());
    return {dhb.degree(v), v};
}

} // namespace dhb
