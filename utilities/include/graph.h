#ifndef UTILITIES_GRAPH_H_
#define UTILITIES_GRAPH_H_

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Graph
{
 public:
    Graph(int input_n);
    ~Graph();
    bool has_edge(int u, int v);
    void add_edge(int u, int v);
    void remove_edge(int u, int v);
    void reduce_degree(int u);
    void set_edge(int u, int v, bool value);
    int get_n();
    int get_degree(int u);
    vector<int>* get_seed_set(int k);
    void set_ordering(int index, int vertex);
    int get_ordering(int index);
    void shuffle_ordering();
    void print();

 protected:
    int n;
    bool *adjacency_matrix;
    int *degree;
    std::vector<int> ordering;
};

class Hardware : public Graph
{
 public:
    Hardware(int input_n, int input_diameter);
    ~Hardware();
  int get_diameter();

 protected:
    int diameter;
};

/* Only contains the dimensions and vertex count */
class Chimera
{
 public:
    Chimera(int input_c, int input_m, int input_n);
    int get_c();
    int get_m();
    int get_n();
    int get_num_vertices();

 protected:
    int c;
  int m;
  int n;
  int num_vertices;
};

#endif
