#pragma once

#include <dhb/dynamic_hashed_blocks.h>

#include <algorithm>

namespace dhb {

using SubmatrixSet = std::vector<Vertex>;

template <typename E> class Submatrix {
  public:
    friend void swap(Submatrix<E>& sub_a, Submatrix<E>& sub_b) noexcept {
        using std::swap;
        swap(sub_a.m_matrix, sub_b.m_matrix);
        swap(sub_a.m_set, sub_b.m_set);
    }

    Submatrix() = delete;

    ~Submatrix() = default;

    Submatrix(Matrix<E>& matrix) : m_matrix(matrix) {}

    Submatrix(Submatrix<E> const& other) : m_matrix(other.m_matrix), m_set(other.m_set) {}

    Submatrix<E>(Submatrix<E>&& other) noexcept : m_matrix(other.m_matrix), m_set() {
        m_set = other.m_set;

        other.m_set.clear();
    }

    Submatrix<E>& operator=(Submatrix<E>& other) {
        if (this != &other) {
            m_set.resize(other.vertices_count());
            std::copy(std::begin(other.m_set), std::end(other.m_set), m_set);

            if (m_matrix != other.m_matrix) {
                m_matrix = other.m_matrix;
            }
        }

        return *this;
    }

    Submatrix<E>& operator=(Submatrix<E>&& other) {
        if (this != &other) {
            m_set = std::move(other.m_set);

            if (&m_matrix != &other.m_matrix) {
                m_matrix = other.m_matrix;
            }
        }

        return *this;
    }

    Submatrix(Matrix<E>& matrix, Vertex const begin, Vertex const end) : m_matrix(matrix) {
        extend_set(begin, end);
    }

    Submatrix(Matrix<E>& matrix, SubmatrixSet const& set) : m_matrix(matrix), m_set(set) {}

    Submatrix(Matrix<E>& matrix, SubmatrixSet&& set) : m_matrix(matrix), m_set(std::move(set)) {}

    void extend_set(Vertex const begin, Vertex const end) {
        Vertex const actual_begin = [&](Vertex const begin) {
            Vertex new_begin = begin;
            std::for_each(std::begin(m_set), std::end(m_set), [&](Vertex v) {
                if (v > new_begin) {
                    new_begin = v;
                    return;
                }
            });
            return new_begin;
        }(begin);

        for (Vertex vertex = actual_begin; vertex < end; ++vertex) {
            m_set.push_back(vertex);
        }
    }

    template <typename F> void for_neighbors_node(Vertex const v, F&& handle) const {
        m_matrix.for_neighbors_node(v, std::move(handle));
    }

    template <typename F> void for_nodes(F&& handle) const {
        std::for_each(std::begin(m_set), std::end(m_set), handle);
    }

    typename Matrix<E>::NeighborView neighbors(Vertex const u) { return m_matrix.neighbors(u); }

    typename Matrix<E>::ConstNeighborView neighbors(Vertex const u) const {
        return typename Matrix<E>::ConstNeighborView{&m_matrix, u};
    }

    Vertex vertices_count() const { return m_set.size(); }

    unsigned long edges_count() const {
        return std::accumulate(
            std::cbegin(m_set), std::end(m_set), 0u,
            [this](unsigned long sum, Vertex const v) { return sum + m_matrix.degree(v); });
    }

    Vertex degree(Vertex const v) const { return m_matrix.degree(v); }

    SubmatrixSet const& set() const { return m_set; }

    Matrix<E>& matrix() { return m_matrix; }

    Matrix<E> const& matrix() const { return m_matrix; }

  private:
    Matrix<E>& m_matrix;
    SubmatrixSet m_set;
};

template <typename E> bool operator==(Submatrix<E> const& a, Submatrix<E> const& b) {
    if (a.set() != b.set() || &a.matrix() != &b.matrix()) {
        return false;
    }

    return true;
}

template <typename E> bool operator!=(Submatrix<E> const& a, Submatrix<E> const& b) {
    return !(a == b);
}

} // namespace dhb
