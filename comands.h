#include <iostream>
#include <unordered_map>
#include <fstream>
#include <stdlib.h>
#include <regex>
#include <cmath>
#include <deque>
#include <ctime>
#include <thread>
#include <chrono>
#include <string>
#include <tuple>
#include <stdexcept>
#include <vector>
#include "nicknames.hpp"
#ifndef COMMANDS_H
#define COMMANDS_H
using namespace std;
tuple<const char*, node_name_t, node_name_t> parse_args(int arg_count, char* arg_vars[]);
matrix_t load_matrix(const char* filename);
template<typename Node>
graph_t<Node> create_graph(const matrix_t& matr) noexcept;
components_t compute_components(const matrix_t& matr) noexcept;

pair<weight_t, route_t> dijkstra(const matrix_t& matr, node_name_t key_from, node_name_t key_to) noexcept;
struct NodeDijkstra {
public:
    double distance;
    bool visited;
    node_name_t predecessor;

    NodeDijkstra() : distance(INF), visited(false), predecessor(0) {}
};

class NodesToBeVisited : public std::vector<graph_t<NodeDijkstra>::iterator>
{
public:
    NodesToBeVisited(graph_t<NodeDijkstra>& graph);

    graph_t<NodeDijkstra>::iterator pop_min_weight();
};


pair<weight_t, route_t> spfa(const matrix_t& matr, node_name_t key_from, node_name_t key_to) noexcept;
struct NodeSPFA {
    weight_t distance;
    unsigned int  update_count;
    bool visited;
    node_name_t predecessor;

    NodeSPFA() : distance(INF), update_count(0), visited(false), predecessor(0) {}
};
class NodeQueue : public std::deque<graph_t<NodeSPFA>::iterator>
{
public:
    void push(const graph_t<NodeSPFA>::iterator& it);

    void pop();
};

#endif COMMANDS_H
