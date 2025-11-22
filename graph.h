#include<iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;
namespace graph
{
    template <typename Key, typename Value, typename Weight>
    class Graph
    {


    public:
        // Вложенный класс Node
        class Node;
        using key_type = Key;
        using value_type = Value;
        using weight_type = Weight;

        //-----------КОНСТРУКТОРЫ И ОПЕРАТОРЫ-------------------//
        Graph() = default;
        Graph(const Graph<key_type, value_type, weight_type>& G) : m_nodes(G.m_nodes) {} // Конструктор копирования
        Graph(Graph<key_type, value_type, weight_type>&& GG)  noexcept : m_nodes(move(GG.m_nodes)) {} //перемещения
        Graph<key_type, value_type, weight_type> operator=(const Graph<key_type, value_type, weight_type>& G);
        Graph<key_type, value_type, weight_type> operator=(Graph<key_type, value_type, weight_type>&& moved) noexcept;
        ~Graph() = default;

        using iterator = typename unordered_map<key_type, Node>::iterator; //итератор для чтения / записи для обхода узлов графа
        using const_iterator = typename unordered_map<key_type, Node>::const_iterator;//только для чтения/онстантный итератор для обхода узлов графа

        // Операторы доступа к элементам графа
        value_type& operator[](key_type key);
        value_type& at(key_type key);// Метод доступа к значению по ключу с проверкой

        bool empty() {
            return m_nodes.empty();
        }
        size_t size() const {
            return m_nodes.size();
        }
        void clear()
        {
            m_nodes.clear();
        }

        void swap(Graph<key_type, value_type, weight_type>& G);

        template <typename key_type, typename value_type, typename weight_type>
        friend void swap(Graph<key_type, value_type, weight_type>& first, Graph<key_type, value_type, weight_type>& second);

        // Методы для работы с узлами и рёбрами графа
        size_t degree_in(key_type key);
        size_t degree_out(key_type key);
        bool loop(key_type key);// Метод, проверяющий наличие петли (рёбра, соединяющего узел с самим собой) для узла

        // Методы вставки и удаления узлов и рёбер

        pair<Graph::iterator, bool> insert(key_type, value_type);

        pair<Graph::iterator, bool> insert_or_assign(key_type, value_type);

        pair<typename Graph::Node::iterator, bool> insert_edge(std::pair<key_type, key_type>, weight_type);

        pair<typename Graph::Node::iterator, bool> insert_or_assign_edge(std::pair<key_type, key_type>, weight_type);

        void delete_node(const key_type& key); //УЗЕЛ
        void delete_edge(pair<key_type, key_type> edge);//РЕБРО


        // Методы поиска и итерации по графу
        iterator find(key_type key);
        iterator begin();
        iterator end();
    private:
        unordered_map<key_type, Node> m_nodes; //УЗЛЫ ГРАФА Хранение узлов графа в хеш-таблице







    };



    // Вложенный класс Node
    template <typename Key, typename Value, typename Weight>
    class Graph<Key, Value, Weight>::Node
    {
    private:
        // Значение узла и его рёбра
        friend Graph<key_type, value_type, weight_type>;
        value_type val;
        unordered_map<key_type, weight_type> m_edges;

    public:
        // Конструкторы и операторы для работы с узлом
        Node() = default;
        explicit Node(value_type value) : val(value) {}
        Node(const Node& G) : val(G.val), m_edges(G.m_edges) {}
        Node(Node&& GG) noexcept;
        Node operator=(const Node& G);
        Node operator=(Node&& GG);
        ~Node() = default;



        bool empty() {
            return m_edges.empty();
        }
        size_t size() {
            return m_edges.size();
        }
        value_type& value() {
            return val;
        }

        using iterator = typename std::unordered_map<key_type, weight_type>::iterator;// Итератор для обхода рёбер узла
        using const_iterator = typename std::unordered_map<key_type, weight_type>::const_iterator;

        iterator begin(); // Возвращает итератор на начало рёбер узла
        iterator end();
        const_iterator begin() const { return m_edges.begin(); }
        const_iterator end() const { return m_edges.end(); }
        void clear() {
            m_edges.clear();
        }
    };
};
