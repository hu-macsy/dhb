#include <dhb/dynamic_hashed_blocks.h>
#include <dhb/integer_log2.h>

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>
#include <mutex>
#include <numeric>
#include <tuple>

namespace dhb {

std::vector<Degree> degrees_from(Edges const& edges) {
    unsigned int count = graph::vertex_count(edges);
    std::vector<Degree> associated_degree(count, 0u);

    for (Edge const& edge : edges) {
        associated_degree[edge.source]++;
    }

    return associated_degree;
}

} // namespace dhb
