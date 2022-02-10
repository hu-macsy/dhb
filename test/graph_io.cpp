#include "graph_io.h"

#include <experimental/filesystem>

dhb::Edges read_graph_unweighted(std::string path) {
    namespace fs = std::experimental::filesystem;

    fs::path graph_path(std::move(path));

    if (!fs::exists(graph_path)) {
        throw std::runtime_error("Path to graph does not exist!");
    }

    std::ifstream graph_input(graph_path);

    dhb::Edges edges;
    read_graph_unweighted(graph_input, [&](unsigned int u, unsigned int v) {
        edges.emplace_back(dhb::Edge{u, {v, 1., 0u}});
    });
    return edges;
}

std::tuple<dhb::Edges, Timestamps> read_temporal_graph(std::string path, bool const weighted) {
    namespace fs = std::experimental::filesystem;

    fs::path graph_path(std::move(path));

    if (!fs::exists(graph_path)) {
        throw std::runtime_error("Path to graph does not exist!");
    }

    std::ifstream graph_input(graph_path);

    std::tuple<dhb::Edges, Timestamps> temporal_edges;

    read_temporal_graph(graph_input, weighted, [&](unsigned int u, unsigned int v, unsigned int t) {
        dhb::Target target{v, {1., 0u}};
        std::get<0>(temporal_edges).emplace_back(dhb::Edge(u, std::move(target)));
        std::get<1>(temporal_edges).push_back(t);
    });

    return temporal_edges;
}
