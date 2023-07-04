#ifndef UTILITIES_EMBEDDING_H_
#define UTILITIES_EMBEDDING_H_

#include <fstream>
#include <memory>

#include "graph.h"

using namespace std;

class Point {
    public:
        int x;
        int y;
        int k;
        Point(int x, int y, int k);
        Point();
        ~Point();
        vector<Point*> * get_neighbors(int min_r, int min_c, int max_r, int max_c);
};

class Tree_Point: public Point {
public:
    Tree_Point(Point* p, int color);
    ~Tree_Point();
    int color;
    int value = 0;
    int expected_color = -1;
    int depth = 0;
    int numOfdes = 0;
    vector<Tree_Point*>* children;
    Tree_Point* parent;
    void insert_child(Tree_Point* p);
    void insert_parent(Tree_Point* p);
};


class EB_Point: public Point {
    public:
        int color;
        EB_Point(int x, int y, int k, int color);
        EB_Point();
        ~EB_Point();
        void print();
        void modify_value(int x, int y, int k, int color);
};

class Pass {
    public:
        vector<int>* old_node;
        int root;
        Pass(int root, vector<int>* old_node);
        ~Pass();
        vector<int>* get_old_node();
        int get_root();
};

class Embedding
{
  public:
    Embedding(int rr, int cc);
    ~Embedding();
    void append(int x, int y, int k, int color);
    void append(EB_Point* p);
    int get_topo_row();
    int get_topo_column();
    void set_topo_row(int r);
    void set_topo_column(int c);
    void assign_selected(vector<EB_Point*>* new_selected);
    vector<int>* get_embedded_nodes();
    vector<Point*>* get_chain(int u);
    vector<EB_Point*>* get_selected();
    void expanding_border();
    int* generate_labelHW();
    void remove_vertex(int x, int y, int k);
    bool contains(int color);
    bool is_color_empty(int color);
    bool is_empty();
    void clear();
    int compute_qubits_used();
    void print();
    void write_to_file(std::string &filename);

 protected:
    vector<EB_Point*>* selected;
    vector<int> * embedded_nodes;
    int topo_row;
    int topo_column;
};

#endif
