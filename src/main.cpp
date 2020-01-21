#include <fstream>
#include <sstream>
#include <iostream>

#include "Graph.h"

void write_times(std::vector<uint_pair> times, const std::string &filename) {
    std::ofstream write(filename.c_str());
    for (unsigned int i = 0; i < times.size(); i++){
        write << times[i].first << "," << times[i].second << std::endl;
    }
}

void read_instance(unsigned int& n, std::vector<uint_pair> &edges, const std::string &filename){
    // edges are numbered from 0 to n - 1, not from 1 to n !!!
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    while (std::getline(targetFile, input_line)){
        std::stringstream line(input_line);
        // ignore family_id
        std::getline(line, entity, ' ');
        if (entity == "p"){
            std::getline(line, entity, ' ');
            std::getline(line, entity, ' ');
            n = std::stoi(entity);
        } else if (entity == "e") {
            std::getline(line, entity, ' ');
            unsigned int i = std::stoi(entity) - 1;
            std::getline(line, entity, ' ');
            unsigned int j = std::stoi(entity) - 1;
            edges.push_back(uint_pair(i, j));
        }
    }
}

std::vector<std::string> read_file_index(const std::string &filename){
    // edges are numbered from 0 to n - 1, not from 1 to n !!!
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    std::vector<std::string> res(0);
    while (std::getline(targetFile, entity)){
        res.push_back(entity);
    }
    return res;
}

int graph_coloring(const std::string& filename, const clock_t &time_limit){
    // we assume absence of duplicate edges...
    // return -1 if we exceed the time limit
    unsigned int n;
    std::vector<uint_pair> edges(0);
    read_instance(n, edges, filename);
    std::cout << n << " vertexes, " << edges.size() <<"/2 edges." << std::endl;
    std::vector<unsigned int> degree(n, 0);

    // we use the maximum degree as upper bound on the nb of colors needed
    for (unsigned int i = 0; i < edges.size(); i++){
        degree[edges[i].first]++;
    }
    unsigned int max = 0;
    for (unsigned int i = 0; i < n; i++)
        if (degree[i] > max)
            max = degree[i];
    int K = max - 1;

    // decrease K until no more solution
    while (K > 0){
        Graph G(n, edges, K);
        status s = G.solve(time_limit);
//        if (s == FOUND){
//            std::vector<unsigned int> sol = G.get_solution(); std::cout << sol[0] << " " << sol[4] << " " << sol[20] << " " << std::endl;
//        }
        if (s == ABORT)
            return -1;
        else if (s == FOUND)
            K--;
        else
            break;
    }
    return K + 1;
}

void benchmark_queens(unsigned int n_max, std::vector<uint_pair> times) {
    for (unsigned int n = 2; n <= n_max; n++){
        clock_t begin = clock();
        clock_t t0 = clock();
        Graph G(n);
        status s  = G.solve(begin + 3*60*CLOCKS_PER_SEC);
        times.push_back(uint_pair(n, clock() - t0));
    }
}


void print_queens_solution(const unsigned int &n, const std::vector<unsigned int>& solution) {
    for (unsigned int i = 0; i < solution.size(); i++) {
            for (unsigned int m = 0; m < solution[i]; m++)
                std::cout << "0 ";
            std::cout << "1 ";
            for (unsigned int m = solution[i] + 1; m < n; m++)
                std::cout << "0 ";
            std::cout << std::endl;
        }
}

void full_benchmark() {
    // full benchmark
    std::vector<uint_pair> times(0);

    // n_queens
    benchmark_queens(40, times);

    // graph coloring
    std::vector<std::string> index = read_file_index("../../index.txt");
    for (unsigned int i = 0; i < index.size(); i++){
        clock_t t0 = clock();
        clock_t begin = t0;
        int a = graph_coloring("../../data/" + index[i], begin + 3*60*CLOCKS_PER_SEC);
        times.push_back(uint_pair(i, clock() - t0));
    }

    // write results
    write_times(times, "../../times_full.csv");
}

int main() {
//    full_benchmark();
//    return 0;
    clock_t begin = clock();

    // test n-queens on an example
//    unsigned int n = 5;
//    Graph G(n);
//    std::cout << "start solve" << std::endl;
//    if(G.solve(begin + 3*60*CLOCKS_PER_SEC) == FOUND) {
//        std::cout << n << "-queens has a solution:" << std::endl;
//        std::vector<unsigned int> solution = G.get_solution();
//        print_queens_solution(n, solution);
//    } else {
//        std::cout << n << "-queens has no solution." << std::endl;
//    }
//    return 0;

    // test graph coloring on an example
//    std::cout << graph_coloring("../../data/anna.col", begin + 3*60*CLOCKS_PER_SEC) << std::endl;
    std::cout << graph_coloring("../../data/queen5_5.col", begin + 3*60*CLOCKS_PER_SEC) << std::endl;
    return 0;
}