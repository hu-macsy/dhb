#include "graph_io.h"

#include <dhb/dynamic_hashed_blocks.h>

#include <omp.h>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>

using namespace dhb;

// https://networkrepository.com/bio-celegans.php
static std::string graph_name{"bio-celegans.mtx"};
static std::string graph_path =
#ifdef DHB_TEST_GRAPH_DIR
    std::string(DHB_TEST_GRAPH_DIR) + "/"
#else
    "test/graphs/"
#endif
    + graph_name;

TEST_CASE("Matrix") {
    Edges edges = read_graph_unweighted(graph_path);

    // we start at 0 so we must have 453 + 1 vertices and not 453 as described on
    // https://networkrepository.com/bio-celegans.php
    unsigned int const bio_celegans_vertex_count = 453 + 1;
    unsigned int const bio_celegans_edge_count = 2025;

    // we consider our graphs to be directed thus we find a maximum degree of 94
    // and not 237 (which it would be in the non directed case) at stated on
    // https://networkrepository.com/bio-celegans.php
    Degree const bio_celegans_max_degree = 94;
    std::vector<Degree> extracted_degrees = degrees_from(edges);

    SECTION("read graph unweighted") {
        CHECK(edges.size() == bio_celegans_edge_count);
        CHECK(graph::vertex_count(edges) == bio_celegans_vertex_count);
        Degree const max_degree =
            *std::max_element(std::begin(extracted_degrees), std::end(extracted_degrees));

        CHECK(max_degree == bio_celegans_max_degree);
    }

    SECTION("sanity check") {
        Matrix<EdgeData> m(std::move(edges));
        CHECK(m.vertices_count() == bio_celegans_vertex_count);

        size_t const max_degree_vertex = std::distance(
            std::begin(extracted_degrees),
            std::max_element(std::begin(extracted_degrees), std::end(extracted_degrees)));

        CHECK(m.degree(max_degree_vertex) == 94);
        CHECK(m.edges_count() == bio_celegans_edge_count);
        CHECK(m.vertices_count() == bio_celegans_vertex_count);
    }

    SECTION("NeighborView") {
        Matrix<EdgeData> m(std::move(edges));
        // N(89) = 13, 31
        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.source() == 89);
        CHECK(nv.degree() == 3);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        auto n89 = nv.cbegin();
        CHECK(n89 != nv.end());

        // yes, we assume an order here
        CHECK(n89->vertex() == 13);
        CHECK(++n89 != nv.cend());
        CHECK(n89->vertex() == 31);
        CHECK(++n89 != nv.cend());
        CHECK(n89->vertex() == 86);

        auto n31 = nv.iterator_to(31);
        CHECK(n31->vertex() == 31);

        CHECK(nv[0].vertex() == 13);
        std::tuple<Matrix<EdgeData>::NeighborView::iterator, bool> r =
            nv.insert(14, EdgeData{10.f, 450});

        CHECK(std::get<0>(r)->vertex() == 14);
        CHECK(nv.degree() == 4);
        CHECK(nv.exists(14));
    }

    SECTION("NeighborView, copy constructor") {
        Matrix<EdgeData> m(std::move(edges));
        // N(89) = 13, 31
        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);

        Matrix<EdgeData>::NeighborView nv_copy{nv};

        CHECK(nv_copy.source() == nv.source());
        CHECK(nv_copy.degree() == nv.degree());
        CHECK(nv_copy.exists(13) == nv.exists(13));
        CHECK(nv_copy.exists(31) == nv.exists(31));
        CHECK(nv_copy.exists(86) == nv.exists(86));
    }

    SECTION("NeighborView, copy assignment") {
        Matrix<EdgeData> m(std::move(edges));
        // N(89) = 13, 31
        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        Matrix<EdgeData>::NeighborView nv_other = m.neighbors(1);
        REQUIRE(nv_other.degree() > 1);

        nv_other = nv;

        CHECK(nv_other.source() == nv.source());
        CHECK(nv_other.degree() == nv.degree());
        CHECK(nv_other.exists(13) == nv.exists(13));
        CHECK(nv_other.exists(31) == nv.exists(31));
        CHECK(nv_other.exists(86) == nv.exists(86));
    }

    SECTION("change target vertex of edge") {
        Matrix<EdgeData> m(std::move(edges));
        constexpr Vertex source = 89;
        constexpr Vertex target = 31;
        constexpr Vertex target_update = 32;
        Matrix<EdgeData>::NeighborView nv = m.neighbors(source);

        REQUIRE(nv.exists(target));
        Matrix<EdgeData>::NeighborView::iterator edge = nv.iterator_to(target);
        edge.get_ptr()->vertex = target_update;

        CHECK(nv.exists(target_update));
        CHECK(!nv.exists(target));
    }

    SECTION("ConstNeighborView") {
        Matrix<EdgeData> const m(std::move(edges));
        // N(89) = 13, 31
        Matrix<EdgeData>::ConstNeighborView nv = m.neighbors(89);
        CHECK(nv.source() == 89);
        CHECK(nv.degree() == 3);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        auto n89 = nv.begin();
        CHECK(n89 != nv.end());

        // yes, we assume an order here
        CHECK(n89->vertex() == 13);
        REQUIRE(++n89 != nv.end());
        CHECK(n89->vertex() == 31);
        REQUIRE(++n89 != nv.end());
        CHECK(n89->vertex() == 86);

        auto n31 = nv.iterator_to(31);
        CHECK(n31->vertex() == 31);

        CHECK(nv[0].vertex() == 13);
    }

    SECTION("sort") {
        Matrix<EdgeData> m(std::move(edges));
        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        nv.insert(e89_14.target.vertex, e89_14.target.data);
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};
        nv.insert(e89_8.target.vertex, e89_8.target.data);

        m.sort(89, [](BlockState<EdgeData>::Entry& a, BlockState<EdgeData>::Entry& b) {
            return a.vertex < b.vertex;
        });

        CHECK(nv.degree() == 5);

        auto n = nv.begin();
        CHECK(n->vertex() == 8);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 13);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 14);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 31);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 86);
        CHECK(++n == nv.end());
    }

    SECTION("removeEdge") {
        Matrix<EdgeData> m(std::move(edges));
        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);

        CHECK(nv.degree() == 3);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(m.removeEdge(89, 13));

        CHECK(!nv.exists(13));

        CHECK(!m.removeEdge(89, 13));
    }

    SECTION("add loop") {
        constexpr Vertex vertex_count = 3;
        Matrix<Weight> m(vertex_count);

        REQUIRE(3 == m.vertices_count());
        REQUIRE(m.insert(0, 1, 1.f, false));
        REQUIRE(1 == m.edges_count());

        CHECK(m.insert(0, 0, 1.f, false));
        REQUIRE(m.neighbors(0).exists(0));
        CHECK(2 == m.edges_count());

        REQUIRE(m.insert(0, 0, 2.f, true));
        REQUIRE(!m.insert(0, 0, 2.f, false));
    }

    SECTION("remove loop") {
        constexpr Vertex vertex_count = 3;
        Matrix<Weight> m(vertex_count);

        REQUIRE(3 == m.vertices_count());
        REQUIRE(m.insert(0, 1, 1.f, false));
        REQUIRE(1 == m.edges_count());

        CHECK(m.insert(0, 0, 1.f, false));
        CHECK(2 == m.edges_count());

        REQUIRE(m.removeEdge(0, 0));

        REQUIRE(!m.neighbors(0).exists(0));
        REQUIRE(1 == m.edges_count());
    }

    SECTION("parallel insert no update single threaded") {
        Matrix<EdgeData> m(std::move(edges));

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto fun = [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); };
        auto get_edge_f = [](Edge const& e) { return e.source; };
        BatchParallelizer<Edge> par;
        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp),
            std::move(fun));

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("parallel insert with update single threaded") {
        Matrix<EdgeData> m(std::move(edges));

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edge e89_13_update{89, Target{13, EdgeData{11.f, 455}}};

        Edges new_edges{e89_14, e89_8, e89_13_update};

        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key = [](Edge e) { return e.source; };
        auto fun = [&](Edge e) {
            std::tuple<dhb::BlockState<dhb::EdgeData>::iterator, bool> insertion_result =
                m.neighbors(e.source).insert(e.target.vertex, e.target.data);
            if (!std::get<1>(insertion_result)) {
                std::get<0>(insertion_result)->data() = e.target.data;
            }
        };
        auto get_edge_f = [](Edge const& e) { return e.source; };

        BatchParallelizer<Edge> par;
        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp),
            std::move(fun));

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(nv.iterator_to(e89_13_update.target.vertex)->data().weight ==
              e89_13_update.target.data.weight);
    }

    SECTION("parallel insert no update dual threaded") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        Matrix<EdgeData> m(graph::vertex_count(edges));

        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key = [](Edge e) { return e.source; };
        auto fun = [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); };
        auto get_edge_f = [](Edge const& e) { return e.source; };
        BatchParallelizer<Edge> par;
        par(edges.begin(), edges.end(), std::move(get_edge_f), std::move(cmp), fun);

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto cmp2 = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key2 = [](Edge e) { return e.source; };
        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp2),
            [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); });

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("parallel insert no update dual threaded, two phase, apply") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        Matrix<EdgeData> m(graph::vertex_count(edges));

        auto key = [](Edge e) { return e.source; };
        auto fun = [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); };
        auto get_edge_f = [](Edge const& e) { return e.source; };
        BatchParallelizer<Edge> par;

        par.apply(edges.begin(), edges.end(), std::move(key), fun);

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto cmp2 = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key2 = [](Edge e) { return e.source; };

        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp2),
            [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); });

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("parallel insert no update dual threaded, two phase, map, distribute") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        Matrix<EdgeData> m(graph::vertex_count(edges));

        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key = [](Edge e) { return e.source; };
        auto fun = [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); };
        auto get_edge_f = [](Edge const& e) { return e.source; };
        BatchParallelizer<Edge> par;

        std::vector<unsigned int> degrees(m.vertices_count());
#pragma omp parallel
        {
            par.distribute(edges.begin(), edges.end(), key);

            par.map(edges.begin(), edges.end(), [&](const auto& e) { ++degrees[e.source]; });

            par.map(edges.begin(), edges.end(), [&](auto const& e) {
                if (degrees[e.source]) {
                    m.preallocate(e.source, degrees[e.source]);
                    degrees[e.source] = 0;
                }

                m.neighbors(e.source).insert(e.target.vertex, e.target.data);
            });
        }

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto cmp2 = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key2 = [](Edge e) { return e.source; };

        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp2),
            [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); });

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("parallel insert with update dual threaded") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        Matrix<EdgeData> m(graph::vertex_count(edges));

        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key = [](Edge e) { return e.source; };
        auto get_edge_f = [](Edge const& e) { return e.source; };

        BatchParallelizer<Edge> par;
        par(edges.begin(), edges.end(), std::move(get_edge_f), std::move(cmp),
            [&](Edge e) { m.neighbors(e.source).insert(e.target.vertex, e.target.data); });

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edge e89_13_update{89, Target{13, EdgeData{11.f, 455}}};

        Edges new_edges{e89_14, e89_8, e89_13_update};

        auto cmp2 = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key2 = [](Edge e) { return e.source; };

        par(new_edges.begin(), new_edges.end(), std::move(get_edge_f), std::move(cmp), [&](Edge e) {
            std::tuple<dhb::BlockState<dhb::EdgeData>::iterator, bool> insertion_result =
                m.neighbors(e.source).insert(e.target.vertex, e.target.data);
            if (!std::get<1>(insertion_result)) {
                std::get<0>(insertion_result)->data() = e.target.data;
            }
        });

        Matrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 5);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(nv.iterator_to(e89_13_update.target.vertex)->data().weight ==
              e89_13_update.target.data.weight);
    }
}
