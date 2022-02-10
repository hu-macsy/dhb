#pragma once

#include <limits>
#include <random>
#include <tuple>
#include <vector>

namespace dhb {

using Weight = double;
using Degree = unsigned int;
#ifdef DHB_64BIT_IDS
using Vertex = uint64_t;
using EdgeID = uint64_t;
#else
using Vertex = uint32_t;
using EdgeID = uint32_t;
#endif

struct EdgeData {
    Weight weight;
    EdgeID id;
};

struct Target {
    Vertex vertex;
    EdgeData data;
};

struct Edge {
    Edge() = default;

    Edge(Vertex _source, Target _target) : source(_source), target(_target) {}

    Vertex source;
    Target target;
};

using Edges = std::vector<Edge>;
using Targets = std::vector<Target>;

Vertex invalidVertex();
Weight invalidWeight();
EdgeID invalidEdgeID();
Edge invalidEdge();
Target invalidTarget();

namespace graph {
Vertex vertex_count(Edges const& edges);
Edge into(Vertex source, Target const&);
Target into(Edge const&);
} // namespace graph

} // namespace dhb
