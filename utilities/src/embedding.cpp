#include "../include/embedding.h"

Embedding::Embedding(int rr, int cc)
{
    selected = new vector<EB_Point*>();
    embedded_nodes = new vector<int>();
    topo_row = rr;
    topo_column = cc;
}

Embedding::~Embedding()
{
    delete[] selected;
    delete[] embedded_nodes;
}

int Embedding::get_topo_row() {
    return this->topo_row;
}

int Embedding::get_topo_column() {
    return this->topo_column;
}

void Embedding::set_topo_row(int r) {
    this->topo_row = r;
}

void Embedding::set_topo_column(int c) {
    this->topo_column = c;
}

void Embedding::assign_selected(vector<EB_Point*>* new_selected) {
    vector<EB_Point*>* former = this->selected;
    this->selected = new_selected;
    for (auto it = former->begin(); it != former->end(); ++it) {
        delete *it;
    }
    delete former;
}

void Embedding::append(int x, int y, int k, int color)
{
    selected->push_back(new EB_Point(x,y,k,color));
    vector<int>::iterator it = std::find (embedded_nodes->begin(), embedded_nodes->end(), color);
    if (it != embedded_nodes->end()) {
        embedded_nodes->push_back(color);
    }
}

void Embedding::append(EB_Point* p)
{
    selected->push_back(p);
    int color = p->color;
    vector<int>::iterator it = std::find (embedded_nodes->begin(), embedded_nodes->end(), color);
    if (it != embedded_nodes->end()) {
        embedded_nodes->push_back(color);
    }
}

vector<int>* Embedding::get_embedded_nodes() {
    return this->embedded_nodes;
}

vector<Point*>* Embedding::get_chain(int u) {
    vector<Point*>* chain = new vector<Point*>();
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        if ((*it)->color == u) {
            chain->push_back(new Point((*it)->x, (*it)->y, (*it)->k));
        }
    }
    return chain;
}

vector<EB_Point*>* Embedding::get_selected() {
    return this->selected;
}

void Embedding::expanding_border() {
    for (auto it = selected->begin(); it != selected->end(); it++) {
        (*it) -> modify_value((*it)->x + 1, (*it)->y + 1, (*it)->k, (*it)->color);
    }
    topo_row = topo_row + 2;
    topo_column = topo_column + 2;
}

int* Embedding::generate_labelHW() {
    int* arr_label_HW = new int[topo_row*topo_column*8]();
    for (int i = 0; i < topo_row*topo_column*8; ++i) {
        arr_label_HW[i] = -1;
    }
    //int* label_HW = new int[topo_row*topo_column*8]();
    for (auto it = selected->begin(); it != selected->end(); it++) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, co = (*it)->color;
        arr_label_HW[x*topo_column*8 + y*8 + k] = co;
    }
    //int* label_HW = arr_label_HW;
    return arr_label_HW;
}

void Embedding::remove_vertex(int x, int y, int k)
{
}

bool Embedding::contains(int color)
{
    return false;
}

bool Embedding::is_color_empty(int color)
{
    return false;
}

bool Embedding::is_empty()
{
    return false;
}

void Embedding::clear()
{
    
}

int Embedding::compute_qubits_used()
{
    return (int)(selected->size());
}

void Embedding::print() {
    std::cout << "Embedding:" << std::endl;
    for (auto it = selected->begin(); it != selected->end(); it++) {
        (*it) -> print();
    }
    std::cout << "Requires "
              << compute_qubits_used()
              << " qubits"
              << std::endl;
}

void Embedding::write_to_file(std::string &filename)
{
}

Point::Point(int x, int y, int k) {
    this->x = x;
    this->y = y;
    this->k = k;
}

Point::Point() {
    this->x = -1;
    this->y = -1;
    this->k = -1;
}

Point::~Point() {
}

Tree_Point::Tree_Point(Point* p, int color) {
    this->x = p->x;
    this->y = p->y;
    this->k = p->k;
    this->color = color;
    this->children = new vector<Tree_Point*>();
    this->parent = NULL;
}

Tree_Point::~Tree_Point() {
    delete this->children;
}

void Tree_Point::insert_child(Tree_Point* p) {
    this->children->push_back(p);
}

void Tree_Point::insert_parent(Tree_Point* p) {
    this->parent = p;
}

EB_Point::EB_Point(int x, int y, int k, int color) {
    this->x = x;
    this->y = y;
    this->k = k;
    this->color = color;
}

EB_Point::EB_Point() {
    this->x = -1;
    this->y = -1;
    this->k = -1;
    this->color = -1;
}

EB_Point::~EB_Point(){
    
}

void EB_Point::print(){
    std::cout << this->x << ' '<< this->y <<' ' <<this->k <<' '<< this->color << std::endl;
}

void EB_Point::modify_value(int x, int y, int k, int color){
    this->x = x;
    this->y = y;
    this->k = k;
    this->color = color;
}

vector<Point*> * Point::get_neighbors(int min_r, int min_c, int max_r, int max_c) {
    vector<Point*> * nb_list = new vector<Point*>();
    for (int i = 0; i < 4; i++) {
        int newx = this->x, newy = this->y, newk;
        if (k > 3) {
            newk = i;
        } else {
            newk = i + 4;
        }
        nb_list->push_back(new Point(newx, newy, newk));
    }
    for (int i = -1; i <= 1; i += 2) {
        int newx = this->x, newy = this->y, newk = this->k;
        if (k > 3) {
            newy = newy + i;
        } else {
            newx = newx + i;
        }
        if (newx > min_r && newx <max_r && newy > min_c && newy < max_c){
            nb_list->push_back(new Point(newx, newy, newk));
        }
    }
    return nb_list;
}

Pass::Pass(int root, vector<int>* old_node) {
    this->old_node = old_node;
    this->root = root;
}

Pass::~Pass() {
    delete[] old_node;
}

vector<int>* Pass::get_old_node() {
    return this->old_node;
}

int Pass::get_root() {
    return this->root;
}
