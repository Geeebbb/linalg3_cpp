#include "comands.h"
#include <vector>

std::tuple<const char*, node_name_t, node_name_t> parse_args(int arg_count, char* arg_vars[])
{
    string data = " ", tmp[3], flag = "";
    for (int i = 1; i < arg_count; i++)
    {
        data += (arg_vars[i] + (string)(" "));
    }

    const regex regexp1(R"(--file (.\w*\.\w*.))", regex_constants::icase);
    const regex regexp2(R"(--from (\w*))", regex_constants::icase);
    const regex regexp3(R"(--to (\w*))", regex_constants::icase);
    const regex regexp[3] = { regexp1, regexp2, regexp3 };
    for (int i = 0; i < 3; i++)
    {
        std::sregex_token_iterator iterator{ data.begin(), data.end(), regexp[i], 1 };
        std::sregex_token_iterator end;
        if ((++iterator) != end)
        {
            if (i == 0)
                flag = "file";
            else if (i == 1)
                flag = "from";
            else if (i == 2)
                flag = "to";
            throw runtime_error("Duplicate argument: " + flag);
        }
    }

    const regex file_arg_regexp(R"( --file )", regex_constants::icase);
    const regex from_arg_regexp(R"( --from )", regex_constants::icase);
    const regex to_arg_regexp(R"( --to )", regex_constants::icase);
    const regex from_node_regexp(R"( --from \d+ )", regex_constants::icase);
    const regex to_node_regexp(R"( --to \d+ )", regex_constants::icase);

    smatch file_arg_match, from_arg_match, to_arg_match, from_node_match, to_node_match;
    std::regex_search(data, file_arg_match, file_arg_regexp);
    std::regex_search(data, from_arg_match, from_arg_regexp);
    std::regex_search(data, to_arg_match, to_arg_regexp);
    std::regex_search(data, from_node_match, from_node_regexp);
    std::regex_search(data, to_node_match, to_node_regexp);

    if (file_arg_match.empty())
    {
        const regex file_err(R"( -*f*i*l*e* )");
        std::sregex_token_iterator iterator{ data.cbegin(), data.cend(), file_err };
        throw runtime_error("Invalid argument:" + string(*iterator));
    }
    else if (from_arg_match.empty())
    {
        const regex from_err(R"( -*f*r*o*m* )");
        std::sregex_token_iterator iterator{ data.cbegin(), data.cend(), from_err };
        throw runtime_error("Invalid argument:" + string(*iterator));
    }
    else if (to_arg_match.empty())
    {
        const regex to_err(R"( -*t*o* )");
        std::sregex_token_iterator iterator{ data.cbegin(), data.cend(), to_err };
        throw runtime_error("Invalid argument:" + string(*iterator));
    }
    else if (from_node_match.empty())
    {
        const regex from_node_err(R"( --from (\w*))");
        std::sregex_token_iterator iterator{ data.cbegin(), data.cend(), from_node_err, 1 };

        if (string(*iterator) == "")
            throw runtime_error("Invalid number of arguments. No argument to \"from\"");
        throw runtime_error("Invalid argument: " + string(*iterator));
    }
    else if (to_node_match.empty())
    {
        const regex to_node_err(R"( --to (\w*))");
        std::sregex_token_iterator iterator{ data.cbegin(), data.cend(), to_node_err, 1 };
        if (string(*iterator) == "")
            throw runtime_error("Invalid number of arguments. No argument to \"to\"");
        throw runtime_error("Invalid argument: " + string(*iterator));
    }

    for (int i = 0; i < 3; i++)
    {
        std::sregex_token_iterator iterator{ data.begin(), data.end(), regexp[i], 1 };
        std::sregex_token_iterator end;
        tmp[i] = string(*iterator);
    }

    return std::make_tuple(tmp[0].c_str(), node_name_t(atoi(tmp[1].c_str())), node_name_t(atoi(tmp[2].c_str())));
}


matrix_t load_matrix(const char* filename)
{
    {
        std::ifstream fin("../" + std::string(filename));
        if (!fin.is_open())
        {
            throw std::runtime_error("Can't find file with this name");
        }

        int rows = 0, columns = 0;
        std::string str;
        bool first_row = true;

        while (std::getline(fin, str, '|'))
        {
            std::stringstream ss(str);
            std::string word;
            int count = 0;

            while (ss >> word)
            {
                count++;
            }

            if (!str.empty() && str != "\n")
            {
                if (first_row)
                {
                    columns = count;
                    first_row = false;
                }

                if (count != columns)
                {
                    throw std::runtime_error("Invalid elements number in row: " + std::to_string(rows + 1) +
                        " (expected: " + std::to_string(columns) +
                        ", actual: " + std::to_string(count) + ")");
                }

                rows++;
            }
        }

        if (rows != columns)
        {
            throw std::runtime_error("Matrix is not square");
        }

        fin.clear();
        fin.seekg(0, std::ios::beg);

        linalg::Matrix<double> mat(rows, columns);
        int count = 0;

        while (std::getline(fin, str, '|'))
        {
            std::stringstream ss(str);
            std::string word;

            while (ss >> word)
            {
                try
                {
                    if (count >= rows * columns)
                    {
                        throw std::runtime_error("More elements found than expected in file");
                    }
                    mat(count / columns, count % columns) = std::stod(word);
                    count++;
                }
                catch (const std::exception& e)
                {
                    throw std::runtime_error("Incorrect file contents");
                }
            }
        }

        return mat;
    }
}

template <typename Node>
graph_t<Node> create_graph(const matrix_t& matr) noexcept
{
    graph_t<Node> gr;
    for (node_name_t i = 0; i < matr.rows(); i++)
    {
        gr.insert(i, Node{});
    }
    for (node_name_t i = 0; i < matr.rows(); i++)
    {
        for (node_name_t j = 0; j < matr.columns(); j++)
        {
            if (matr(i, j) != 0)
            {
                gr.insert_edge({ i, j }, matr(i, j));
            }
        }
    }
    return gr;
}

/*-------------------Algoritms on graph--------------------*/

std::vector<node_name_t> topology_sort(graph_t<bool> graph_inverted)
{
    std::vector<node_name_t> vec;
    int flag = 0, check = 0, x = 1;
    for (int i = 0; i < graph_inverted.size(); i++)
    {
        graph_t<bool>::iterator p = graph_inverted.find(i);
        for (int j = 0; j < vec.size(); j++)
        {
            if (vec[j] == i)
            {
                flag = 1;
                break;
            }
        }
        if (flag)
        {
            flag = 0;
            continue;
        }
        std::vector<node_name_t> tmp;

        while (x == 1)
        {

            for (graph_t<bool>::Node::iterator it = p->second.begin(); it != p->second.end(); ++it)
            {
                if (!graph_inverted.at(it->first))
                {
                    p->second.value() = true;
                    tmp.push_back(p->first);
                    p = graph_inverted.find(it->first);
                    check = 1;
                    break;
                }
            }

            if (!check)
            {
                p->second.value() = true;
                vec.insert(vec.begin(), p->first);
                if (tmp.empty())
                {
                    break;
                }
                p = graph_inverted.find(tmp[tmp.size() - 1]);
                tmp.pop_back();
            }
            check = 0;
        }
        flag = 0;
    }
    return vec;
}

components_t compute_components(const matrix_t& matr) noexcept
{
    // graph_t<bool> graph_initial = create_graph<bool>(matr);

    graph_t<bool> graph_initial;
    thread th([&graph_initial, &matr]()
        { graph_initial = create_graph<bool>(matr); });

    graph_t<bool> graph_inverted = create_graph<bool>(transpose(matr));
    std::vector<node_name_t> sorted_nodes = topology_sort(graph_inverted);

    th.join();

    components_t res;

    int flag = 0, x = 1, check = 0, counter = -1;
    for (int i = 0; i < sorted_nodes.size(); i++)
    {
        for (int j = 0; j < res.size(); j++)
        {
            for (int k = 0; k < res[j].size(); k++)
            {
                if (sorted_nodes[i] == res[j][k])
                {
                    flag = 1;
                    break;
                }
            }
            if (flag == 1)
                break;
        }
        if (flag == 1)
        {
            flag = 0;
            continue;
        }
        counter++;
        res.push_back(vector<node_name_t>());
        graph_t<bool>::iterator p = graph_initial.find(sorted_nodes[i]);
        std::vector<node_name_t> tmp;
        while (x == 1)
        {

            for (graph_t<bool>::Node::iterator it = p->second.begin(); it != p->second.end(); ++it)
            {
                if (!graph_initial.at(it->first))
                {
                    p->second.value() = true;
                    tmp.push_back(p->first);
                    p = graph_initial.find(it->first);
                    check = 1;
                    break;
                }
            }

            if (!check)
            {
                p->second.value() = true;
                if (!res[counter].empty() && res[counter].front() == p->first)
                    res[counter].erase(res[counter].begin());

                res[counter].insert(res[counter].begin(), p->first);

                if (tmp.empty())
                {
                    break;
                }
                p = graph_initial.find(tmp[tmp.size() - 1]);
                tmp.pop_back();
            }
            check = 0;
        }
        flag = 0;
    }
    return res;
}

NodesToBeVisited::NodesToBeVisited(graph_t<NodeDijkstra>& graph)
{
    for (graph_t<NodeDijkstra>::iterator it = graph.begin(); it != graph.end(); ++it)
    {
        this->push_back(it);
    }
}

graph_t<NodeDijkstra>::iterator NodesToBeVisited::pop_min_weight()
{
    double min = this->front()->second.value().distance;
    vector<graph_t<NodeDijkstra>::iterator>::iterator pos = this->begin();

    for (vector<graph_t<NodeDijkstra>::iterator>::iterator it = this->begin(); it != this->end(); ++it)
    {
        if ((*it)->second.value().distance < min)
        {

            min = (*it)->second.value().distance;
            pos = it;
        }
    }
    graph_t<NodeDijkstra>::iterator tmp = *pos;
    (*pos)->second.value().visited = true;
    this->erase(pos);
    return tmp;
}

void dijkstra_step(graph_t<NodeDijkstra>& graph, NodesToBeVisited& nodes_to_be_visited)
{
    graph_t<NodeDijkstra>::iterator pos = nodes_to_be_visited.pop_min_weight();

    for (graph_t<NodeDijkstra>::Node::iterator it = pos->second.begin(); it != pos->second.end(); ++it)
    {
        if (pos->second.value().distance + it->second < graph.at(it->first).distance)
        {
            graph.at(it->first).distance = pos->second.value().distance + it->second;
            graph.at(it->first).predecessor = pos->first;
        }
    }
}

std::pair<weight_t, route_t> dijkstra(const matrix_t& matr, node_name_t key_from, node_name_t key_to) noexcept
{
    graph_t<NodeDijkstra> graph = create_graph<NodeDijkstra>(matr);

    graph.at(key_from).distance = 0.0;

    NodesToBeVisited nodes_to_be_visited(graph);

    while (!nodes_to_be_visited.empty())
        dijkstra_step(graph, nodes_to_be_visited);

    if (graph.at(key_to).distance == INF)
    {
        std::vector<node_name_t> tmp{ key_from, key_to };
        std::pair<weight_t, route_t> res(INF, tmp);
        return res;
    }

    route_t tmp;
    node_name_t x = key_to;
    tmp.insert(tmp.begin(), x);
    while (x != key_from)
    {
        x = graph.at(x).predecessor;
        tmp.insert(tmp.begin(), x);
    }

    std::pair<weight_t, route_t> res(graph.at(key_to).distance, tmp);
    return res;
}

void NodeQueue::push(const graph_t<NodeSPFA>::iterator& it)
{
    it->second.value().visited = true;
    this->push_back(it);
}

void NodeQueue::pop()
{
    this->front()->second.value().visited = false;
    this->pop_front();
}

void spfa_step(graph_t<NodeSPFA>& graph, NodeQueue& nodes_in_work)
{
    node_name_t pos = nodes_in_work.front()->first;
    if (graph.at(pos).update_count == graph.size())
    {
        graph.at(pos).distance = -INF;
    }

    for (graph_t<NodeSPFA>::Node::iterator it = graph.find(pos)->second.begin(); it != graph.find(pos)->second.end(); ++it)
    {
        if (graph.at(pos).distance + it->second < graph.at(it->first).distance)
        {
            graph.at(it->first).distance = graph.at(pos).distance + it->second;
            nodes_in_work.push(graph.find(it->first));
            graph.at(it->first).predecessor = pos;
            graph.at(it->first).update_count++;
        }
    }

    nodes_in_work.pop();
}

std::pair<weight_t, route_t> spfa(const matrix_t& matr, node_name_t key_from, node_name_t key_to) noexcept
{
    graph_t<NodeSPFA> graph = create_graph<NodeSPFA>(matr);

    auto node_start_it = graph.find(key_from);
    node_start_it->second.value().distance = 0.0;

    NodeQueue nodes_in_work;
    nodes_in_work.push(node_start_it);

    while (!nodes_in_work.empty())
        spfa_step(graph, nodes_in_work);

    if (graph.at(key_to).distance == -INF)
    {
        std::vector<node_name_t> tmp{ key_from, key_to };
        std::pair<weight_t, route_t> res(-INF, tmp);
        return res;
    }
    if (graph.at(key_to).distance == INF)
    {
        std::vector<node_name_t> tmp{ key_from, key_to };
        std::pair<weight_t, route_t> res(INF, tmp);
        return res;
    }
    std::vector<node_name_t> route;
    node_name_t x = key_to;
    route.insert(route.begin(), x);
    while (x != key_from)
    {
        x = graph.at(x).predecessor;
        route.insert(route.begin(), x);
    }

    return std::pair<weight_t, route_t>(graph.at(key_to).distance, route);
}

void print_result(components_t comp, std::pair<weight_t, route_t> route) noexcept
{
    if (route.first == INF)
    {
        cout << "There is no way from '" << route.second[0] << "' to '" << route.second[1] << "' !";
    }
    else if (route.first == -INF)
    {
        cout << "There is a negative cycle between '" << route.second[0] << "' and '" << route.second[1] << "' !";
    }

    else
    {
        cout << "Shortest route (weight " << route.first << "): ";
        for (int i = 0; i < route.second.size(); i++)
        {
            cout << route.second[i] << " ";
        }
    }
    cout << "\nStrongly connected components:" << endl;
    for (int i = 0; i < comp.size(); i++)
    {
        cout << i << ") ";
        for (int j = 0; j < comp[i].size(); j++)
        {
            cout << comp[i][j] << " ";
        }
        cout << endl;
    }
}

void create_big_matrix()
{
    std::ofstream out;
    out.open("../big_matrix.txt");

    for (int i = 0; i < 2000; i++)
    {
        out << "|";
        for (int j = 0; j < 2000; j++)
        {
            int x = rand() % (10 - 0 + 1) + 0;
            out << x << " ";
        }
        out << "|\n";
    }

    out.close();
}

