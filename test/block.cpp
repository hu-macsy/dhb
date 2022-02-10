#include <catch2/catch_test_macros.hpp>
#include <dhb/block.h>
#include <dhb/buckets.h>
#include <dhb/graph.h>

using namespace dhb;

#ifndef DHB_SYSTEM_ALLOCATOR
TEST_CASE("Block") {
    unsigned int const bytes_per_entry = sizeof(Target);
    unsigned int const degree = 2;

    BlockManager bm(bytes_per_entry);

    SECTION("BlockState Sanity Check") {
        BlockHandle bh = bm.allocate_block(degree);
        BlockState<EdgeData> bs(bh);

        CHECK(bs.bsize() >= degree);
        CHECK(bs.degree() == 0);
    }

    Edge edge(0, Target{1, EdgeData{10.f, 0}});

    SECTION("BlockState insert") {
        BlockHandle bh = bm.allocate_block(degree);
        BlockState<EdgeData> bs(bh);

        std::tuple<BlockState<EdgeData>::iterator, bool> insertion_result =
            bs.insert(edge.target.vertex, edge.target.data);

        CHECK(std::get<1>(insertion_result));
        CHECK(bs.degree() == 1);

        BlockState<EdgeData>::iterator it = std::get<0>(insertion_result);
        CHECK(it->vertex() == edge.target.vertex);
        CHECK(it->data().weight == edge.target.data.weight);
        CHECK(it->data().id == edge.target.data.id);
    }

    SECTION("BlockState remove") {
        BlockHandle bh = bm.allocate_block(degree);
        BlockState<EdgeData> bs(bh);

        std::tuple<BlockState<EdgeData>::iterator, bool> insertion_result =
            bs.insert(edge.target.vertex, edge.target.data);

        REQUIRE(std::get<1>(insertion_result));
        CHECK(bs.degree() == 1);
        CHECK(bs.remove(edge.target.vertex));
        CHECK(bs.degree() == 0);

        Edge edge_2(1, Target{2, EdgeData{20.f, 1}});

        insertion_result = bs.insert(edge.target.vertex, edge.target.data);
        REQUIRE(std::get<1>(insertion_result));
        CHECK(bs.degree() == 1);
        CHECK(!bs.full());

        insertion_result = bs.insert(edge_2.target.vertex, edge_2.target.data);

        REQUIRE(std::get<1>(insertion_result));
        CHECK(bs.degree() == 2);
        CHECK(bs.remove(edge.target.vertex));
        CHECK(bs.degree() == 1);

        CHECK(bs.begin()->vertex() == edge_2.target.vertex);
    }
}
#endif
