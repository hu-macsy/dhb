#include <dhb/graph.h>

#include <algorithm>

namespace dhb {

Vertex invalidVertex() { return std::numeric_limits<Vertex>::max(); }
Weight invalidWeight() { return std::numeric_limits<float>::infinity(); }
EdgeID invalidEdgeID() { return std::numeric_limits<EdgeID>::max(); }

Edge invalidEdge() { return Edge(invalidVertex(), invalidTarget()); }

Target invalidTarget() { return {invalidVertex(), {invalidWeight(), invalidEdgeID()}}; }

Vertex graph::vertex_count(Edges const& edges) {
    // Determine n as the maximal node ID.
    Vertex n = 0;
    for (Edge const& edge : edges) {
        Vertex const vertex = std::max(edge.source, edge.target.vertex);
        n = std::max(n, vertex);
    }

    n += 1;
    return n;
}

} // namespace dhb
