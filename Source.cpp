#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include "comands.hpp"
#include <filesystem> 
namespace fs = std::filesystem;

using namespace std;
using namespace graph;
void test_load_matrix();
void test_create_graph();
void test_algorithms();

int main(int argc, char* argv[]) {

    try {
        //test_load_matrix();
        const char* filename = "Matrix.txt"; // Замените на имя вашего тестового файла
        try
        {
            // Загружаем матрицу из файла
            linalg::Matrix<double> mat = load_matrix(filename);

            // Выводим матрицу на экран
            cout << "Loaded matrix:" << endl;
            for (int i = 0; i < mat.rows(); ++i)
            {
                for (int j = 0; j < mat.columns(); ++j)
                {
                    cout << setw(10) << mat(i, j) << " ";
                }
                cout << endl;
            }
        }
        catch (const std::exception& e)
        {
            cerr << "Error: " << e.what() << endl;
        }

        test_create_graph();
        test_algorithms();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

void test_load_matrix() {
    try {
        // Test 1: Valid matrix
        matrix_t mat1 = load_matrix("valid.txt");
        cout << "Loaded valid matrix successfully." << endl;

        // Test 2: Invalid matrix (non-square)
        try {
            matrix_t mat2 = load_matrix("invalid.txt");
        }
        catch (const std::runtime_error& e) {
            cout << "Caught expected exception for invalid matrix: " << e.what() << endl;
        }

        // Test 3: Incorrect file
        try {
            matrix_t mat3 = load_matrix("non_existent_file.txt");
        }
        catch (const std::runtime_error& e) {
            cout << "Caught expected exception for non-existent file: " << e.what() << endl;
        }
    }
    catch (const std::exception& e) {
        throw std::runtime_error(string("test_load_matrix: ") + e.what());
    }
}

void test_create_graph() {
    try {
        // Example matrix
        matrix_t mat = {
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1},
            {1, 0, 0, 0}
        };

        Graph<string, bool, double> graph;
        graph.insert("A", true);
        graph.insert("B", false);
        graph.insert("C", true);
        graph.insert("D", false);

        graph.insert_edge({ "A", "B" }, 1.0);
        graph.insert_edge({ "B", "C" }, 1.0);
        graph.insert_edge({ "C", "D" }, 1.0);
        graph.insert_edge({ "D", "A" }, 1.0);

        cout << "Graph created successfully from matrix." << endl;

        // Print the graph (optional)
        for (const auto& node : graph) {
            cout << "Node " << node.first << ": ";
            for (const auto& edge : node.second) {
                cout << " -> " << edge.first << " (weight: " << edge.second << ")";
            }
            cout << endl;
        }
    }
    catch (const std::exception& e) {
        throw std::runtime_error(string("test_create_graph: ") + e.what());
    }
}

void test_algorithms() {
    try {
        // Example matrix for algorithms
        matrix_t mat = {
            {0, 1, 4, 0},
            {0, 0, 4, 0},
            {0, 0, 0, 5},
            {0, 0, 0, 0}
        };

        // Compute components
        components_t components = compute_components(mat);
        cout << "Computed strongly connected components successfully." << endl;

        // Dijkstra's algorithm
        auto dijkstra_result = dijkstra(mat, 0, 3);
        cout << "Dijkstra result: Weight = " << dijkstra_result.first << ", Path = ";
        for (auto node : dijkstra_result.second) {
            cout << node << " ";
        }
        cout << endl;
        // SPFA algorithm
        auto spfa_result = spfa(mat, 0, 3);
        cout << "SPFA result: Weight = " << spfa_result.first << ", Path = ";
        for (auto node : spfa_result.second) {
            cout << node << " ";
        }
        cout << endl;
    }
    catch (const std::exception& e) {
        throw std::runtime_error(string("test_algorithms: ") + e.what());
    }
}