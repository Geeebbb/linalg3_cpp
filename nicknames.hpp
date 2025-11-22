#include <set>
#include "graph.hpp"
#include "Matrix.hpp"
#include <limits>
#include <vector>
using namespace linalg;
using namespace graph;
using namespace std;
using node_name_t = unsigned int;
using weight_t = double;
using matrix_t = linalg::Matrix<weight_t>;
template <typename T>
using graph_t = graph::Graph<node_name_t, T, weight_t>;
using components_t = vector<vector<node_name_t>>;
using route_t = vector<node_name_t>;

#define INF std::numeric_limits<double>::infinity()
#define EPS std::numeric_limits<double>::epsilon()
