#include "graph.h"

using namespace std;

namespace graph
{


    /*-------------------------------------------*/



    template<typename Key, typename Value, typename Weight>
    Graph<Key, Value, Weight>::Node::Node(Node&& GG) noexcept : val(move(GG.val)), m_edges(move(GG.m_edges)) {}

    //      Node
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::Node Graph<Key, Value, Weight>::Node::operator=(const Node& G) {
        if (!this->empty())
        {
            this->clear();
        }
        if (this != &G)
        {
            val = G.val;
            m_edges = G.m_edges;
        }
        return *this;
    }
    //     Node
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::Node Graph<Key, Value, Weight>::Node::operator=(Node&& GG) {
        if (this != &GG) {
            val = move(GG.val);
            m_edges = move(GG.m_edges);
        }
        return *this;
    }

    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::Node::iterator Graph<Key, Value, Weight>::Node::begin() {
        return m_edges.begin();
    }

    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::Node::iterator Graph<Key, Value, Weight>::Node::end() {
        return m_edges.end();
    }




    /*------------------ GRAPH--------------------*/



    template <typename key_type, typename value_type, typename weight_type>
    void Graph<key_type, value_type, weight_type>::swap(Graph<key_type, value_type, weight_type>& G)
    {
        this->m_nodes.swap(G.m_nodes);
    }

    template <typename key_type, typename value_type, typename weight_type>
    void swap(Graph<key_type, value_type, weight_type>& first, Graph<key_type, value_type, weight_type>& second)
    {
        first.m_nodes.swap(second.m_nodes);
    }

    // =
    template <typename key_type, typename value_type, typename weight_type>
    Graph<key_type, value_type, weight_type> Graph<key_type, value_type, weight_type>::operator=(const Graph<key_type, value_type, weight_type>& G)
    {
        if (!this->empty())
        {
            this->m_nodes.clear();
        }
        this->m_nodes = G.m_nodes;
        return *this;
    }

    template <typename key_type, typename value_type, typename weight_type>
    Graph<key_type, value_type, weight_type> Graph<key_type, value_type, weight_type>::operator=(Graph<key_type, value_type, weight_type>&& GG) noexcept
    {
        this->m_nodes = std::move(GG.m_nodes);
        return *this;
    }


    //      []
    template <typename key_type, typename value_type, typename weight_type>
    value_type& Graph<key_type, value_type, weight_type>::operator[](key_type key)

    {
        return this->m_nodes[key].value();
    }




    //      
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::value_type& Graph<Key, Value, Weight>::at(key_type key) {
        return m_nodes.at(key).val;
    }


    //    
    template<typename Key, typename Value, typename Weight>
    size_t Graph<Key, Value, Weight>::degree_in(key_type key) {
        if (m_nodes.find(key) != m_nodes.end())
        {
            size_t tmp = 0;
            for (iterator it = this->m_nodes.begin(); it != this->m_nodes.end(); ++it)
            {
                tmp += it->second.m_edges.count(key);
            }
            return tmp;
        }
        return 0;
    }

    //    
    template<typename Key, typename Value, typename Weight>
    size_t Graph<Key, Value, Weight>::degree_out(key_type key) {
        if (m_nodes.find(key) != m_nodes.end()) {
            return m_nodes[key].size();
        }
        return 0;
    }
    // 
    template<typename Key, typename Value, typename Weight>
    bool Graph<Key, Value, Weight>::loop(key_type key) {
        auto it = m_nodes.find(key);
        if (it == m_nodes.end()) {
            throw runtime_error("Node with the specified key does not exist.");
        }

        for (const auto& edge_pair : it->second.m_edges) {
            if (edge_pair.first == key) {
                return true; //  
            }
        }

        return false; //   
    }



    //    
    template <typename key_type, typename value_type, typename weight_type>
    pair<typename Graph<key_type, value_type, weight_type>::iterator, bool> Graph<key_type, value_type, weight_type>::insert(key_type key, value_type value) {
        pair<typename Graph<key_type, value_type, weight_type>::iterator, bool> tmp = this->m_nodes.insert({ key, Node(value) });
        return tmp;
    }
    //     orrrrr
    template <typename key_type, typename value_type, typename weight_type>
    pair<typename Graph<key_type, value_type, weight_type>::iterator, bool> Graph<key_type, value_type, weight_type>::insert_or_assign(key_type key, value_type value)
    {
        pair <typename Graph<key_type, value_type, weight_type>::iterator, bool> result = this->m_nodes.insert_or_assign(key, value);
        return result;
    }

    template <typename key_type, typename value_type, typename weight_type>
    pair<typename Graph<key_type, value_type, weight_type>::Node::iterator, bool> Graph<key_type, value_type, weight_type>::insert_edge(pair<key_type, key_type> keys, weight_type weight)
    {
        try
        {
            this->at(keys.first);
        }
        catch (const exception& e)
        {
            throw runtime_error("Error! There is not first key in this graph");
        }

        try
        {
            this->at(keys.second);
        }
        catch (const exception& e)
        {
            throw runtime_error("Error! There is not second key in this graph");
        }

        pair<typename Graph<key_type, value_type, weight_type>::Node::iterator, bool> result = this->m_nodes.at(keys.first).m_edges.insert({ keys.second, weight });
        return result;
    }

    template <typename key_type, typename value_type, typename weight_type>
    pair<typename Graph<key_type, value_type, weight_type>::Node::iterator, bool> Graph<key_type, value_type, weight_type>::insert_or_assign_edge(pair<key_type, key_type> keys, weight_type weight)
    {
        try
        {
            this->at(keys.first);
        }
        catch (const exception& e)
        {
            throw runtime_error("Error! There is not first key in this graph");
        }

        try
        {
            this->at(keys.second);
        }
        catch (const exception& e)
        {
            throw runtime_error("Error! There is not second key in this graph");
        }
        pair<typename Graph<key_type, value_type, weight_type>::Node::iterator, bool> result = this->m_nodes.at(keys.first).m_edges.insert_or_assign(keys.second, weight);
        return result;
    }

    template <typename key_type, typename value_type, typename weight_type>
    void Graph<key_type, value_type, weight_type>::delete_edge(pair<key_type, key_type> keys)
    {
        this->m_nodes.at(keys.first).m_edges.erase(keys.second);
    }

    template <typename key_type, typename value_type, typename weight_type>
    void Graph<key_type, value_type, weight_type>::delete_node(const key_type& key)
    {
        try
        {
            this->m_nodes.at(key);
        }
        catch (const exception& e)
        {
            throw runtime_error("Error! There s not such Node in this graph");
        }

        for (iterator it = this->m_nodes.begin(); it != this->m_nodes.end(); ++it)
        {
            it->second.m_edges.erase(key);
        }
        this->m_nodes.erase(key);
    }

    //    
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::iterator Graph<Key, Value, Weight>::find(key_type key) {
        return m_nodes.find(key);
    }

    //    
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::iterator Graph<Key, Value, Weight>::begin() {
        return m_nodes.begin();
    }

    //    
    template<typename Key, typename Value, typename Weight>
    typename Graph<Key, Value, Weight>::iterator Graph<Key, Value, Weight>::end() {
        return m_nodes.end();
    }

}