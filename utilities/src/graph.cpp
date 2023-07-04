#include "../include/graph.h"

using namespace std;

Graph::Graph(int input_n)
{
    n = input_n;
    adjacency_matrix = new bool[n*n]();
    degree = new int[n]();
    ordering = std::vector<int>(n, -1);
}

Graph::~Graph() {
    delete[] adjacency_matrix;
    delete[] degree;
}

bool Graph::has_edge(int u, int v)
{
    if (u < 0 || v < 0) {
        return false;
    }
    return adjacency_matrix[u * n + v];
}

void Graph::add_edge(int u, int v)
{
    degree[u] = degree[u] + 1;
    degree[v] = degree[v] + 1;
    adjacency_matrix[u * n + v] = true;
    adjacency_matrix[v * n + u] = true;
}

void Graph::remove_edge(int u, int v)
{
    degree[u] = degree[u] - 1;
    degree[v] = degree[v] - 1;
    adjacency_matrix[u * n + v] = false;
    adjacency_matrix[v * n + u] = false;
}

void Graph::reduce_degree(int u) {
    degree[u] --;
}


void Graph::set_edge(int u, int v, bool value)
{
    adjacency_matrix[u * n + v] = value;
}

int Graph::get_n()
{
    return n;
}

int Graph::get_degree(int u)
{
    return degree[u];
}

vector<int>* Graph::get_seed_set(int k) {
    int n = this->get_n();
    int i = 0;
    vector<int>* res = new vector<int>();
    for (int i = 0; i < n; ++i) {
        res->push_back(i);
    }
    int degree_P[n];
    for (int i = 0; i < n; ++i) {
        degree_P[i] = get_degree(i);
    }
    for (int i = k; i < n; ++i) {
        int minn = 10000000, saved = -1;
        for (auto it = res->begin(); it != res->end(); ++it) {
            if (degree_P[*it] < minn) {
                minn = degree_P[*it];
                saved = *it;
            }
        }
        res->erase(find(res->begin(), res->end(), saved));
        for (auto it = res->begin(); it != res->end(); ++it) {
            if (has_edge(*it, saved)) {
                degree_P[*it]--;
            }
        }
    }
    return res;
}

void Graph::set_ordering(int index, int vertex)
{
    ordering[index] = vertex;
}

int Graph::get_ordering(int index)
{
    return ordering[index];
}

void Graph::print() {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << has_edge(i,j) << " ";
        }
        std::cout << std::endl;
    }
}


Hardware::Hardware(int input_n, int input_diameter) : Graph(input_n)
{
    diameter = input_diameter;
}

Hardware::~Hardware() {}

int Hardware::get_diameter()
{
    return diameter;
}


Chimera::Chimera(int input_c, int input_m, int input_n)
{
    c = input_c;
    m = input_m;
    n = input_n;
    num_vertices = 2 * c * m * n;
}

int Chimera::get_c()
{
    return c;
}

int Chimera::get_m()
{
    return m;
}

int Chimera::get_n()
{
    return n;
}

int Chimera::get_num_vertices()
{
    return num_vertices;
}
