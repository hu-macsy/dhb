#include <catch2/catch_test_macros.hpp>
#include <dhb/vec.h>

TEST_CASE("Vec") {
    dhb::Vec<int> v;

    REQUIRE(std::get<1>(v.insert(2, 1)));
    REQUIRE(std::get<1>(v.insert(3, 1)));
    REQUIRE(std::get<1>(v.insert(5, 1)));
    REQUIRE(std::get<1>(v.insert(7, 1)));
    REQUIRE(std::get<1>(v.insert(11, 1)));

    SECTION("iterate") {
        std::array<int, 5> indices{2, 3, 5, 7, 11};

        bool entries_are_equal =
            std::equal(v.begin(), v.end(), indices.begin(), indices.end(),
                       [](auto ent, int idx) -> bool { return ent.vertex() == idx; });
        CHECK(entries_are_equal);

        bool all_edge_data_is_1 =
            std::all_of(v.begin(), v.end(), [](auto ent) -> bool { return ent.data() == 1; });
        CHECK(all_edge_data_is_1);
    }

    SECTION("lookup") {
        REQUIRE(v.find(1) == v.end());
        REQUIRE(v.find(2) != v.end());
        REQUIRE(v.find(3) != v.end());
        REQUIRE(v.find(4) == v.end());
        REQUIRE(v.find(5) != v.end());
    }

    SECTION("insert") {
        auto r1 = v.insert(1, 2);
        auto r2 = v.insert(2, 2);
        auto r4 = v.insert(4, 2);

        // Validate the returned bool (whether the value was inserted or not).
        CHECK(std::get<1>(r1));
        CHECK(!std::get<1>(r2));
        CHECK(std::get<1>(r4));

        // Validate the returned iterator.
        CHECK(std::get<0>(r1)->data() == 2);
        CHECK(std::get<0>(r2)->data() == 1);
    }

    SECTION("copy") {
        auto w = v;

        bool vertices_are_copied =
            std::equal(v.begin(), v.end(), w.begin(), w.end(),
                       [](auto ve, auto we) -> bool { return ve.vertex() == we.vertex(); });
        CHECK(vertices_are_copied);

        bool data_is_copied =
            std::equal(v.begin(), v.end(), w.begin(), w.end(),
                       [](auto ve, auto we) -> bool { return ve.data() == we.data(); });
        CHECK(data_is_copied);
    }

    SECTION("erase") {
        v.erase(v.find(11));
        v.erase(v.find(5));
        std::array<int, 3> indices{2, 3, 7};

        bool erase_took_place =
            std::equal(v.begin(), v.end(), indices.begin(), indices.end(),
                       [](auto ent, int idx) -> bool { return ent.vertex() == idx; });
        CHECK(erase_took_place);
    }

    SECTION("clear") {
        v.clear();
        CHECK(v.begin() == v.end());
    }
}
