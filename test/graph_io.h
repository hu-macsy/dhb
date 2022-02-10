#pragma once

#include <dhb/graph.h>

#include <fstream>
#include <iostream>
#include <sstream>

using Timestamp = unsigned int;
using Timestamps = std::vector<Timestamp>;

dhb::Edges read_graph_unweighted(std::string path);

template <typename F> void read_graph_unweighted(std::istream& ins, F fn) {
    std::string line;
    bool seen_header = false;
    while (std::getline(ins, line)) {
        if (line.empty())
            continue;
        if (line.front() == '%')
            continue;
        if (line.front() == '#')
            continue;

        std::istringstream ls(line);
        dhb::Vertex u = 0;
        dhb::Vertex v = 0;
        if (!(ls >> u >> v))
            throw std::runtime_error("Parse error while reading input graph");

        if (!seen_header) {
            seen_header = true;
            continue;
        }

        fn(u, v);
    }
}

std::tuple<dhb::Edges, Timestamps> read_temporal_graph(std::string path, bool weighted);

template <typename F> void read_temporal_graph(std::istream& ins, bool const weighted, F fn) {
    std::string line;
    bool seen_header = false;
    while (std::getline(ins, line)) {
        if (line.empty())
            continue;
        if (line.front() == '%')
            continue;

        std::istringstream ls(line);
        dhb::Vertex u = 0;
        dhb::Vertex v = 0;
        unsigned int t = 0;
        float w = 0.;

        if (weighted) {
            if (!(ls >> u >> v >> w >> t)) {
                throw std::runtime_error("Parse error while reading weighted input graph");
            }
        } else {
            if (!(ls >> u >> v >> t)) {
                throw std::runtime_error("Parse error while reading unweighted input graph");
            }
        }

        if (!seen_header) {
            seen_header = true;
            continue;
        }

        fn(u, v, t);
    }
}
