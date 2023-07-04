
#include <iostream>
#include <fstream>
#include <queue>
#include <stdlib.h>
#include <time.h>
#include <memory>
#include <thread>

#include "utilities/include/graph.h"
#include "utilities/include/io.h"
#include "utilities/include/embedding.h"


using namespace std;

void delete_pointer(vector<Point*>* vec) {
    for (auto it = vec->begin(); it != vec->end(); ++it) {
        delete *it;
    }
    vec->clear();
    delete vec;
}

void delete_queue(queue<Point*>* q) {
    while (!q->empty()) {
        delete q->front();
        q->pop();
    }
    delete q;
}

void delete_tree(Tree_Point* r) {
    for (auto it = r->children->begin(); it != r->children->end(); ++it) {
        delete_tree(*it);
    }
    delete r;
}

void initialize_embedding(Embedding* eb, vector<int>* seed_set, int seed_limit) {
    int c_count = 0;
    for (auto it = seed_set->begin(); it != seed_set->end(); ++it){
        if (c_count <= 2) {
            eb->append(0, 0, c_count, *it);
            eb->append(0, 0, c_count + 4, *it);
        }
        else {
            if (c_count == 3) {
                eb->append(0, 0, 3, *it);
            } else {
                eb->append(0, 0, 7, *it);
            }
        }
        c_count++;
    }
}

void extract_order(Graph* P, vector<Pass*>* list_passes, vector<int>* seed_set, int seed_limit) {
    int n = P->get_n();
    vector<int> orders;
    while (true) {
        double maxn = -10000000;
        int saved = -1;
        for (int i = 0; i < n; ++i) {
            double countt = 0;
            vector<int>::iterator it = find(seed_set->begin(), seed_set->end(), i);
            if (it == seed_set->end()) {
                for (auto seed_it = seed_set->begin(); seed_it != seed_set->end(); ++seed_it) {
                    if (P->has_edge(i, *seed_it)) {
                        countt = countt - (double)(distance(seed_set->begin(), seed_it)) + n;
                        //countt = countt + 1;
                    }
                }
                //cout<<countt<<' ';
                if (countt > maxn) {
                    maxn = countt;
                    saved = i;
                }
                /*if (countt == maxn && rand()%2 == 1) {
                    maxn = countt;
                    saved = i;
                }*/
            }
        }
        //cout<<endl;
        if (saved == -1) {
            break;
        } else {
            seed_set->push_back(saved);
            vector<int>* old_node = new vector<int>();
            for (auto seed_it = seed_set->begin(); seed_it != seed_set->end(); ++seed_it) {
                if (P->has_edge(saved, *seed_it)) {
                    old_node->push_back(*seed_it);
                }
            }
            Pass* ps = new Pass(saved, old_node);
            list_passes->push_back(ps);
        }
    }
    for (int i = 0; i < seed_limit; i++){
        seed_set->erase(seed_set->begin());
    }
}

void expanding_chain(Embedding* embedding, Graph* P, Pass* ps) {
    int n = P->get_n(), topo_row = embedding->get_topo_row(), topo_column = embedding ->get_topo_column();
    
    
    
    int* label_HW = embedding->generate_labelHW();
    //cout << *(label_HW + 37);
    //cout << *(label_HW + 0);
    
    vector<int>* old_node = ps->get_old_node();
    vector<Point*> *free_neighbors = new vector<Point*>[old_node->size()];
    int len_old = (int)(old_node->size());
    bool* check_HW = new bool[topo_row*topo_column*len_old];
    for (int i = 0; i < topo_row*topo_column*len_old; ++i) {
        check_HW[i] = false;
    }
    vector<EB_Point*> *selected = embedding->get_selected();
    
    for (auto it = selected->begin(); it != selected->end(); it++) {
        auto pos = find(old_node->begin(), old_node->end(), (*it)->color);
        int co = (*it)->color;
        
        if (pos != old_node->end()) {
            int idx = (int)(pos - old_node->begin());
            vector<Point*> *neighbors = (*it) -> get_neighbors(-1, -1, topo_row, topo_column);
            for (auto p_it = neighbors->begin(); p_it != neighbors->end(); ++p_it) {
                int x = (*p_it)->x, y = (*p_it)->y, k = (*p_it)->k;
                if (not check_HW[x*topo_column*len_old+y*len_old+idx] && label_HW[x*topo_column*8 + y*8 + k] == -1) {
                    check_HW[x*topo_column*len_old+y*len_old+idx] = true;
                    free_neighbors[idx].push_back(new Point((*p_it)->x,(*p_it)->y,(*p_it)->k));
                }
            }
            delete_pointer(neighbors);
        }
    }
    
    int idx = 0;
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        int gap = P->get_degree(*it) - (int)(free_neighbors[idx].size());
        if (gap > -3) {
            for (auto p_it = free_neighbors[idx].begin(); p_it != free_neighbors[idx].end(); ++p_it) {
                int x = (*p_it)->x, y = (*p_it)->y, k = (*p_it)->k;
                if (label_HW[x*topo_column*8 + y*8 + k] == -1){
                    embedding->append(x, y, k, *it);
                    label_HW[x*topo_column*8 + y*8 + k] = *it;
                    gap = gap - 1;
                }
            }
        }
        idx ++;
    }
    
    idx = 0;
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        for (auto p_it = free_neighbors[idx].begin(); p_it != free_neighbors[idx].end(); ++p_it) {
            delete *p_it;
        }
        free_neighbors[idx].clear();
        idx++;
    }
    free_neighbors->clear();
    delete [] free_neighbors;
    delete[] label_HW;
    delete[] check_HW;
}

void BFS(Embedding* embedding, int it, int idx, int topo_row, int topo_column, int* label_HW, int* f) {
    /*long long dem = 0;
    for (int i = 0; i < 10000000000; ++i){
        ++dem;
    }*/
    clock_t time_Start = clock();

    vector<Point*>* cardinity = embedding->get_chain(it);
    queue<Point*>* q = new queue<Point*>;
    //int visited[topo_row][topo_column][8];
    int* visited = new int[topo_row*topo_column*8]();
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                visited[i*topo_column*8+j*8+k] = -1;
            }
        }
    }
    
    for (auto it = cardinity->begin(); it != cardinity->end(); ++it) {
        int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
        visited[x*topo_column*8+y*8+k] = 0;
        q->push(new Point(x,y,k));
    }
    
    while (! q->empty()) {
        Point* head = q->front();
        q->pop();
        vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
        for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
            int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
            if (visited[x*topo_column*8+y*8+k] == -1 && label_HW[x*topo_column*8+y*8+k] == -1) {
                visited[x*topo_column*8+y*8+k] = visited[head->x*topo_column*8+head->y*8+head->k] + 1;
                f[idx*topo_row*topo_column*8 + x*topo_column*8 + y*8 + k] = visited[x*topo_column*8+y*8+k];
                q->push(new Point((*it)->x,(*it)->y,(*it)->k)) ;
            }
        }
        delete head;
        delete_pointer(neighbors);
    }

    delete [] visited;
    delete_queue(q);
    delete_pointer(cardinity);
}

pair<int,int> calculate_DP(Tree_Point* curr_root) {
    if (curr_root->children->size() == 0) {
        if (curr_root->color == -1) {
            curr_root->value = 0;
            return make_pair(0, 0);
        } else {
            curr_root->value = 1;
            curr_root->expected_color = curr_root->color;
            return make_pair(1, 1);
        }
    }
    int sum_val = 0, sum_numdes = 0;
    for (auto it = curr_root->children->begin(); it != curr_root->children->end(); ++it) {
        pair<int,int> saved = calculate_DP(*it);
        sum_val += saved.first;
        sum_numdes += saved.second;
        if (saved.second > 0) {
            curr_root->expected_color = (*it)->expected_color;
            curr_root->depth = (*it)->depth + 1;
        }
    }
    
    if (sum_val > 0) {
        sum_val = sum_val + 1;
    }
    curr_root->value = sum_val;
    curr_root->numOfdes = sum_numdes;
    return make_pair(sum_val, sum_numdes);
}

Tree_Point* find_dynamic_tree(Point* good_pos, int topo_row, int topo_column, int* label_HW, vector<int>* old_node, int n) {
    queue<Tree_Point*> q;
    int* visited = new int[topo_row*topo_column*8]();
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                visited[i*topo_column*8+j*8+k] = -1;
            }
        }
    }
    bool check_eb[n], final_check = true;
    for (int i = 0; i < n; ++i) {
        check_eb[i] = false;
    }
    Tree_Point* root = new Tree_Point(good_pos, -1);
    q.push(root);
    visited[root->x*topo_column*8+root->y*8+root->k] = 0;
    while (! q.empty()) {
        Tree_Point* head = q.front();
        q.pop();
        vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
        for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
            int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
            if (visited[x*topo_column*8+y*8+k] == -1){
                int color_tail = label_HW[x*topo_column*8+y*8+k];
                if (color_tail >= 0) {
                    if (find(old_node->begin(), old_node->end(), color_tail) != old_node->end() && !check_eb[color_tail]) {
                        check_eb[color_tail] = true;
                        Tree_Point* new_point = new Tree_Point(*it, color_tail);
                        head->insert_child(new_point);
                        new_point->insert_parent(head);
                    }
                } else {
                    visited[x*topo_column*8+y*8+k] = visited[head->x*topo_column*8+head->y*8+head->k] + 1;
                    Tree_Point* new_point = new Tree_Point(*it, -1);
                    head->insert_child(new_point);
                    new_point->insert_parent(head);
                    q.push(new_point);
                }
            }
        }
        delete_pointer(neighbors);
    }
    delete [] visited;
    pair<int,int> res = calculate_DP(root);
    return root;
}

bool is_good_pos(Point* A, Point* B) {
    return (A->x == B->x && A->y == B->y);
}

void allocate_node(Tree_Point* curr_root, int new_color, Embedding* embedding, bool check_first, bool check_alone, int curr_depth, Graph* P) {
    bool chek_alone = check_alone;
    if (curr_root->numOfdes == 0 || curr_root->color != -1){
        return;
    }
    if (curr_root->numOfdes > 1 || check_first){
        EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, new_color);
        embedding->append(new_embed);
    } else {
        int old_color = curr_root->expected_color;
        if (check_alone) {
            EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
            embedding->append(new_embed);
        } else {
            if (P->get_degree(old_color) - curr_depth < 0) {
                EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
                embedding->append(new_embed);
                chek_alone = true;
            } else {
                //cout<<P->get_degree(old_color)<<' '<<P->get_degree(new_color) <<endl;

                int rand_val = rand()% (P->get_degree(old_color) + P->get_degree(new_color));
                if (rand_val < P->get_degree(old_color) || is_good_pos(curr_root, curr_root->parent)) {
                    EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
                    embedding->append(new_embed);
                    chek_alone = true;
                } else {
                    EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, new_color);
                    embedding->append(new_embed);
                }
            }
        }
    }

    for (auto it = curr_root->children->begin(); it != curr_root->children->end(); ++it) {
        allocate_node(*it, new_color, embedding, false, chek_alone, curr_depth + 1, P);
    }
    return;
}

void allocate_node_test(Tree_Point* curr_root, int new_color, vector<EB_Point*>* new_selected, bool check_first, bool check_alone, int curr_depth, Graph* P) {
    bool chek_alone = check_alone;
    if (curr_root->numOfdes == 0 || curr_root->color != -1){
        return;
    }
    if (curr_root->numOfdes > 1 || check_first){
        EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, new_color);
        new_selected->push_back(new_embed);
    } else {
        int old_color = curr_root->expected_color;
        if (check_alone) {
            EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
            new_selected->push_back(new_embed);
        } else {
            if (P->get_degree(old_color) - curr_depth < 0) {
                EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
                new_selected->push_back(new_embed);
                chek_alone = true;
            } else {
                int rand_val = rand()% (P->get_degree(old_color) + P->get_degree(new_color) - curr_depth + 1);
                if (rand_val < P->get_degree(old_color) || is_good_pos(curr_root, curr_root->parent)) {
                    EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, old_color);
                    new_selected->push_back(new_embed);
                    chek_alone = true;
                } else {
                    EB_Point* new_embed = new EB_Point(curr_root->x, curr_root->y, curr_root->k, new_color);
                    new_selected->push_back(new_embed);
                }
            }
        }
    }

    for (auto it = curr_root->children->begin(); it != curr_root->children->end(); ++it) {
        allocate_node_test(*it, new_color, new_selected, false, chek_alone, curr_depth + 1, P);
    }
    return;
}

bool find_embedding_own(Embedding* embedding, Graph* P, Pass* ps) {
    int n = P->get_n(), topo_row = embedding->get_topo_row(), topo_column = embedding ->get_topo_column();
    int new_color = ps->get_root();
    int* label_HW = embedding->generate_labelHW();
    vector<EB_Point*> *selected = embedding->get_selected();
    vector<int>* old_node = ps->get_old_node();
    //cout<<old_node->size()<<endl;
    //cout<<"REAL topo: "<<topo_row<<' '<<topo_column<<endl;
    
    int max_size = topo_row*topo_column*8*old_node->size();
    int* f = new int[max_size];
   // int* sum_f = new int[max_size];
   // bool* check_f = new bool[max_size];
    for (int i = 0; i < max_size; ++i) {
        f[i] = -1;
 //       sum_f[i] = 0;
 //       check_f[i] = false;
    }
    
    int num_threads = 10;
    
    int idx = 0;
    auto it = old_node->begin();
    while (it != old_node->end()) {
        vector<std::thread> threads;
        clock_t time_Start = clock();
        for (int i = 0; i < num_threads; ++i) {
            if (it == old_node->end()) {
                num_threads = i;
                break;
            }
            time_Start = clock();
            std::thread t(BFS, embedding, *it, idx, topo_row, topo_column, label_HW, f);
            threads.push_back(move(t));
            idx++;
            it++;
        }
        for (int i = 0; i < num_threads; ++i) {
            threads[i].join();
        }
    }
    
    /*int idx = 0;
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        //vector<Point*>* cardinity = embedding->get_chain(*it);
        BFS(embedding, *it, idx, topo_row, topo_column, label_HW, f);
        idx ++;
        //delete_pointer(cardinity);
    }*/
    
    int maxn = 10000000;
    Point* good_pos = new Point();
    Point* good_pos_2 = new Point();
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                bool check_embedded = true;
                int sum_embedded = 0, idx = 0;
                for (auto it = old_node->begin(); it != old_node->end(); ++it) {
                    if (label_HW[i*topo_column*8+j*8+k] >= 0 || f[idx*topo_row*topo_column*8 + i*topo_column*8 + j*8 + k] == -1) {
                        check_embedded = false;
                        break;
                    } else {
                        sum_embedded += f[idx*topo_row*topo_column*8 + i*topo_column*8 + j*8 + k];
//                        sum_f[i*topo_row*8+j*8+k] += f[idx*topo_row*topo_column*8 + i*topo_column*8 + j*8 + k];
                    }
                    idx = idx + 1;
                }
//                check_f[i*topo_row*8+j*8+k] = check_embedded;
                if (check_embedded && sum_embedded < maxn) {
                    maxn = sum_embedded;
                    good_pos->x = i;
                    good_pos->y = j;
                    good_pos->k = k;
                }
            }
        }
    }
    
    //cout << "Good Pos: " <<good_pos->x<< ' ' << good_pos->y << ' '<<good_pos->k<<' '<<maxn<<endl;
    /*cout << "Point Test: "<<endl;
    int epsilon = 50;
    int min_size = 1000000;
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                int idx = i*topo_row*8+j*8+k;
                if (label_HW [idx] == -1 && check_f[idx] && sum_f[idx] >= maxn - epsilon) {
                    Point* estimate_pos = new Point();
                    estimate_pos->x = i;
                    estimate_pos->y = j;
                    estimate_pos->k = k;
                    Tree_Point* root = find_dynamic_tree(estimate_pos, topo_row, topo_column, label_HW, old_node, n);
                    vector<EB_Point*>* selected = embedding->get_selected();
                    vector<EB_Point*>* new_selected = new vector<EB_Point*>();
                    for (auto it = selected->begin(); it != selected->end(); it++) {
                        new_selected->push_back(new EB_Point((*it)->x, (*it)->y, (*it)->k, (*it)->color));
                    }
                    bool check_first = true, check_alone = false;
                    int curr_depth = 0;
                    allocate_node_test(root, new_color, new_selected, check_first, check_alone, curr_depth, P);
                    if (i == good_pos->x && j == good_pos->y && k == good_pos->k) {
                        cout << i<< ' ' <<j<<' '<<k<<' '<<new_selected->size() << endl;
                    }
                    if (new_selected->size() < min_size) {
                        min_size = new_selected->size();
                        good_pos_2->x = i;
                        good_pos_2->y = j;
                        good_pos_2->k = k;
                    }
                }
            }
            
        }
        
    }
    cout << "Min Size: "<< min_size <<' '<<maxn<<endl;
    good_pos = good_pos_2;*/
    if (good_pos->x == -1) {
        delete good_pos;
        delete[] f;
        delete[] label_HW;
        //delete[] sum_f;
        //delete [] check_f;
        return false;
    } else {
        Tree_Point* root = find_dynamic_tree(good_pos, topo_row, topo_column, label_HW, old_node, n);
        bool check_first = true, check_alone = false;
        int curr_depth = 0;
        allocate_node(root, new_color, embedding, check_first, check_alone, curr_depth, P);
        delete_tree(root);
        delete good_pos;
        delete[] label_HW;
        delete[] f;
        //delete[] sum_f;
        //delete [] check_f;
        return true;
    }
}

void refresh(queue<Point*>*q, int* visited, Point* tail, vector<Point*>* cardinity, int topo_row, int topo_column, int limit_visited, const int num_HW) {
    while (!q->empty()) {
        q->pop();
    }
    Point* head = new Point(tail->x, tail->y, tail->k);
    while (visited[head->x*topo_column*8+head->y*8+head->k] != num_HW*limit_visited) {
        cardinity->push_back(new Point(head->x, head->y, head->k));
        vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
        for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
            int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
            if (visited[x*topo_column*8+y*8+k] == visited[head->x*topo_column*8+head->y*8+head->k] - 1) {
                head = *it;
                break;
            }
        }
        //delete_pointer(neighbors);
    }
    
    for (auto it = cardinity->begin(); it != cardinity->end(); ++it) {
        int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
        visited[x*topo_column*8+y*8+k] = (limit_visited+1)*num_HW;
        q->push(new Point(x,y,k));
    }
    
}

bool BFS_TST(vector<Point*>* cardinity, int topo_row, int topo_column, int* label_HW, vector<int>* old_node, int n, int idx) {
    queue<Point*>* q = new queue<Point*>();
    int* visited = new int[topo_row*topo_column*8];
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                visited[i*topo_column*8+j*8+k] = -1;
            }
        }
    }
    
    bool dx_old[n];
    for (int i = 0; i < n; ++i) {
        dx_old[i] = false;
    }
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        dx_old[*it] = true;
    }
    
    int limit_visited = 0;
    const int num_HW = topo_row*topo_column*8;
    
    for (auto it = cardinity->begin(); it != cardinity->end(); ++it) {
        int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
        visited[x*topo_column*8+y*8+k] = limit_visited*num_HW;
        q->push(new Point(x,y,k));
    }
    dx_old[idx] = false;
    
    while (! q->empty()) {
        Point* head = q->front();
        q->pop();
        vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
        for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
            int x = (*it) -> x, y = (*it) -> y, k = (*it) -> k;
            //cout << x*topo_column*8+y*8+k << endl;
            if (label_HW[x*topo_column*8+y*8+k] >= 0 && label_HW[head->x*topo_column*8+head->y*8+head->k] >= 0) {
                continue;
            }
            if (visited[x*topo_column*8+y*8+k] < limit_visited*num_HW) {
                if (label_HW[x*topo_column*8+y*8+k] == -1) {
                    visited[x*topo_column*8+y*8+k] = visited[head->x*topo_column*8+head->y*8+head->k] + 1;
                    q->push(*it);
                } else {
                    if (dx_old[label_HW[x*topo_column*8+y*8+k]]) {
                        if (limit_visited == 0) {
                            cardinity->clear();
                        }
                        refresh(q, visited, head, cardinity, topo_row, topo_column, limit_visited, num_HW);
                        dx_old[label_HW[x*topo_column*8+y*8+k]] = false;
                        limit_visited++;
                        break;
                    }
                }
            }
        }
        //delete_pointer(neighbors);
    }
    bool check = true;
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        if (dx_old[*it]) {
            check = false;
            break;
        }
    }
    delete_queue(q);
    return check;
}

bool find_embedding_2(Embedding* embedding, Graph* P, Pass* ps) {
    int n = P->get_n(), topo_row = embedding->get_topo_row(), topo_column = embedding ->get_topo_column();
    int new_color = ps->get_root();
    int* label_HW = embedding->generate_labelHW();
    vector<EB_Point*> *selected = embedding->get_selected();
    vector<int>* old_node = ps->get_old_node();
    
    bool check = true;
    int minn = 1000000;
    vector<Point*>* cardinity;
    for (auto it = old_node->begin(); it != old_node->end(); ++it) {
        cardinity = embedding->get_chain(*it);
        check = BFS_TST(cardinity, topo_row, topo_column, label_HW, old_node, P->get_n(), *it);
        if (check && cardinity->size()< minn) {
            minn = cardinity->size();
        }
        delete_pointer(cardinity);
        if (!check) {
            break;
        }
        break;
    }
    
    cout << "CARDINITY " << ' ' <<minn<< endl;
    for (auto it = cardinity->begin(); it != cardinity->end(); ++it) {
        cout << (*it)->x << ' ' << (*it)->y << ' ' << (*it)->k << endl;
    }
    //cout << cardinity->size() << endl;
    return check;
}

bool check_exist(vector<EB_Point*>* dx, Point* p) {
    bool check = false;
    for (auto it = dx->begin(); it != dx->end(); ++it) {
        if (p->x == (*it)->x && p->y == (*it)->y && p->k == (*it)->k) {
            check = true;
            break;
        }
    }
    return check;
}

bool test_chain_connection(Embedding* embedding, Graph* P) {
    int n = P->get_n();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    int* label_HW = embedding->generate_labelHW();
    vector<EB_Point*>* selected = embedding->get_selected();
    int nums[n];
    bool dx[topo_row][topo_column][8];
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                dx[i][j][k] = false;
            }
        }
    }
    
    for (int i = 0; i < n; ++i) {
        nums[i] = 0;
    }
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        nums[(*it)->color]++;
    }
    bool total_check = true;
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if (dx[x][y][k]) {
            continue;
        } else {
            queue<Point*> q;
            q.push(new Point(x,y,k));
            dx[x][y][k] = true;
            int countt = 0;
            while (!q.empty()) {
                Point* head = q.front();
                q.pop();
                ++countt;
                vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
                for (auto it_p = neighbors->begin(); it_p != neighbors->end(); ++it_p) {
                    int xp = (*it_p)->x, yp = (*it_p)->y, kp = (*it_p)->k;
                    if (label_HW[xp*topo_column*8+yp*8+kp] == color && !dx[xp][yp][kp]) {
                        q.push(new Point(xp,yp,kp));
                        dx[xp][yp][kp] = true;
                    }
                }
                delete_pointer(neighbors);
            }
            
            if (countt != nums[color]) {
                cout << "Wrong at: " << x << ' ' <<y << ' ' << k <<' '<<color <<endl;
                total_check = false;
            }
        }
    }
    delete [] label_HW;
    return total_check;
}

bool test_global_connection(Embedding* embedding, Graph* P) {
    int n = P->get_n();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    int* label_HW = embedding->generate_labelHW();
    vector<EB_Point*>* selected = embedding->get_selected();
    bool total_check = true;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (P->has_edge(i, j)) {
                bool check = false;
                for (auto it = selected->begin(); it != selected->end(); it++) {
                    if (check) {
                        break;
                    }
                    int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
                    if (color == i || color == j) {
                        vector<Point*>* neighbors = (*it)->get_neighbors(-1, -1, topo_row, topo_column);
                        for (auto it_p = neighbors->begin(); it_p != neighbors->end(); ++it_p) {
                            int xp = (*it_p)->x, yp = (*it_p)->y, kp = (*it_p)->k, colorp = label_HW[xp*topo_column*8+yp*8+kp];
                            if ((colorp == i && color == j) || (colorp == j && color == i)) {
                                check = true;
                                break;
                            }
                        }
                        delete_pointer(neighbors);
                    }
                }
                if (!check) {
                    total_check = false;
                    cout << "Wrong at: (" <<i<<','<<j<<')'<<endl;
                }
            }
        }
    }
    delete [] label_HW;
    return total_check;
}

bool test_onetomany_connection(Embedding* embedding, Graph* P) {
    int n = P->get_n();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    int* label_HW = embedding->generate_labelHW();
    bool total_check = true;
    vector<EB_Point*>* selected = embedding->get_selected();
    int se_size = (int)(selected->size());
    for (int i = 0; i < se_size; ++i) {
        EB_Point* p1 = selected->at(i);
        for (int j = i + 1; j < se_size; ++j) {
            EB_Point* p2 = selected->at(j);
            if (p1->x == p2->x && p1->y == p2->y && p1->k == p2->k) {
                cout <<"Wrong at: ("<<p1->x<<','<<p1->y<<','<<p1->k<<','<<p1->color<<"), ("<<p2->x<<','<<p2->y<<','<<p2->k<<','<<p2->color<<')'<<endl;
                total_check = false;
            }
        }
    }
    delete [] label_HW;
    return total_check;
}


bool test_chain_connection_lightweight(int n, Embedding* embedding, EB_Point* removed_point, int* label_HW) {
    int nums = 0, topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    vector<EB_Point*>* selected = embedding->get_selected();
    EB_Point* saved_node = NULL;
    for (auto it = selected->begin(); it != selected->end(); it++) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k;
        if (label_HW[x*topo_column*8+y*8+k] == removed_point->color) {
            ++nums;
            saved_node = *it;
        }
    }
    if (nums == 0) {
        return false;
    }
    vector<EB_Point*>* dx = new vector<EB_Point*>();
    dx->push_back(new EB_Point(saved_node->x,saved_node->y,saved_node->k,saved_node->color));
    int idx = 0, color = saved_node->color, countt = 0;
    while (idx < dx->size()) {
        EB_Point* head = dx->at(idx);
        ++countt;
        vector<Point*>* neighbors = head->get_neighbors(-1, -1, topo_row, topo_column);
        for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
            int x = (*it)->x, y = (*it)->y, k = (*it)->k;
            if (label_HW[x*topo_column*8+y*8+k] == color && ! check_exist(dx, *it)) {
                dx->push_back(new EB_Point(x,y,k,color));
            }
        }
        delete_pointer(neighbors);
        ++idx;
    }
    for (auto it = dx->begin(); it != dx->end(); ++it) {
        delete *it;
    }
    delete dx;
    return (countt == nums);
}

bool test_global_connection_lightweight(Embedding* embedding, EB_Point* removed_point, int* label_HW, Graph* P) {
    vector<EB_Point*>* selected = embedding->get_selected();
    int color = removed_point->color, topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    vector<Point*>* neighbors = removed_point->get_neighbors(-1, -1, topo_row, topo_column);
    vector<int>* check_colors = new vector<int>();
    for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k;
        int neighbor_color = label_HW[x*topo_column*8+y*8+k];
        if (neighbor_color != -1 && neighbor_color != color && P->has_edge(color, neighbor_color)) {
            check_colors->push_back(neighbor_color);
        }
    }
    
    for (auto it = selected->begin(); it != selected->end(); it++) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k;
        if (label_HW[x*topo_column*8+y*8+k] == removed_point->color) {
            vector<Point*>* neighbors = (*it)->get_neighbors(-1, -1, topo_row, topo_column);
            for (auto it_p = neighbors->begin(); it_p != neighbors->end(); ++it_p) {
                int x = (*it_p)->x, y = (*it_p)->y, k = (*it_p)->k;
                if (label_HW[x*topo_column*8+y*8+k] != -1 && label_HW[x*topo_column*8+y*8+k] != color) {
                    auto it_check = find(check_colors->begin(), check_colors->end(), label_HW[x*topo_column*8+y*8+k]);
                    if (it_check != check_colors->end()) {
                        check_colors->erase(it_check);
                    }
                }
            }
        }
    }
    delete_pointer(neighbors);
    long sz = check_colors->size();
    delete check_colors;
    
    return sz == 0;
}

void remove_trash(int n, Embedding* embedding, int last_selected_size, Graph* P) {
    vector<EB_Point*>* selected = embedding->get_selected();
    auto it = selected->begin() + last_selected_size;
    int* label_HW = embedding -> generate_labelHW();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    while (it != selected->end()) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k;
        int old_val = label_HW[x*topo_column*8+y*8+k];
        label_HW[x*topo_column*8+y*8+k] = -1;
        if (! test_chain_connection_lightweight(n, embedding, *it, label_HW) || ! test_global_connection_lightweight(embedding, *it, label_HW, P)) {
            label_HW[x*topo_column*8+y*8+k] = old_val;
            it++;
        } else {
            selected->erase(it);
        }
    }
    delete [] label_HW;
}

void expanding_vertical(int pivot, vector<int>* old_node, Embedding* embedding, Graph* P) {
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    vector<EB_Point*>* selected = embedding->get_selected();
    int* label_HW = embedding->generate_labelHW();
    int end_point = (int)(selected->size());
    for (int i = 0; i < end_point; ++i) {
        auto it = selected->begin() + i;
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if ((*it)->x > pivot) {
            (*it)->modify_value(x + 1, y, k, color);
        } else {
            if ((*it)->x == pivot) {
                if ((*it)->k <= 3) {
                    int prev_co = label_HW[x*topo_column*8+y*8+k];
                    int next_co = -1;
                    if ((*it)->x <= topo_row - 2) {
                        next_co = label_HW[(x+1)*topo_column*8+y*8+k];
                    }
                    int co = -1;
                    if (P->has_edge(prev_co, next_co) || (prev_co == next_co && prev_co != -1)) {
                        if (find(old_node->begin(), old_node->end(), prev_co) != old_node->end()) {
                            co = prev_co;
                        } else {
                            if (find(old_node->begin(), old_node->end(), next_co) != old_node->end()) {
                                co = next_co;
                            } else {
                                if (rand()%(P->get_degree(prev_co) + P->get_degree(next_co) + 1) < P->get_degree(prev_co)) {
                                    co = prev_co;
                                } else {
                                    co = next_co;
                                }
                            }
                        }
                
                        embedding->append(x + 1, y, k, co);
                    }
                    
                }
            }
        }
    }
    embedding->set_topo_row(topo_row + 1);
    delete [] label_HW;
}

void expanding_horizontal(int pivot, vector<int>* old_node, Embedding* embedding, Graph* P) {
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    vector<EB_Point*>* selected = embedding->get_selected();
    int* label_HW = embedding->generate_labelHW();
    int end_point = (int)(selected->size());
    for (int i = 0; i < end_point; ++i) {
        auto it = selected->begin() + i;
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if ((*it)->y > pivot) {
            (*it)->modify_value(x, y + 1, k, color);
        } else {
            if ((*it)->y == pivot) {
                if ((*it)->k > 3) {
                    int prev_co = label_HW[x*topo_column*8+y*8+k];
                    int next_co = -1;
                    if ((*it)->y <= topo_column - 2) {
                        next_co = label_HW[x*topo_column*8+(y+1)*8+k];
                    }
                    int co = -1;
                    if (P->has_edge(prev_co, next_co) || (prev_co == next_co && prev_co != -1)) {
                        if (find(old_node->begin(), old_node->end(), prev_co) != old_node->end()) {
                            co = prev_co;
                        } else {
                            if (find(old_node->begin(), old_node->end(), next_co) != old_node->end()) {
                                co = next_co;
                            } else {
                                if (rand()%(P->get_degree(prev_co) + P->get_degree(next_co) + 1) < P->get_degree(prev_co)) {
                                    co = prev_co;
                                } else {
                                    co = next_co;
                                }
                            }
                        }
                    
                        embedding->append(x, y + 1, k, co);
                    }
                }
            }
        }
    }
    embedding->set_topo_column(topo_column + 1);
    delete [] label_HW;
}

void expanding_total(Embedding* embedding, Pass* ps, Graph* P) {
    vector<EB_Point*>* selected = embedding->get_selected();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    vector<int>* old_node = ps->get_old_node();
    
    int rmax = -1, cmax = -1;
    vector<int> *r_select = new vector<int>();
    vector<int> *c_select = new vector<int>();
    r_select->push_back(-1);
    c_select->push_back(-1);
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        rmax = max(rmax, (*it)->x);
        cmax = max(cmax, (*it)->y);
    }
    int r = rmax + 1, c = cmax + 1;
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        if ((*it)->k > 3) {
            if ((*it)->x < r - 1 && find(r_select->begin(), r_select->end(), (*it)->x) == r_select->end()) {
                r_select->push_back((*it)->x);
            }
        } else {
            if ((*it)->y < c - 1 && find(c_select->begin(), c_select->end(), (*it)->y) == c_select->end()) {
                c_select->push_back((*it)->y);
            }
        }
        rmax = max(rmax, (*it)->x);
        cmax = max(cmax, (*it)->y);
    }
    sort(r_select->begin(), r_select->end(), greater<int>());
    sort(c_select->begin(), c_select->end(), greater<int>());

    for (auto it = r_select->begin(); it != r_select->end(); ++it) {
        expanding_vertical(*it, old_node, embedding, P);
    }
    
    for (auto it = c_select->begin(); it != c_select->end(); ++it) {
        expanding_horizontal(*it, old_node, embedding, P);
    }
    delete r_select;
    delete c_select;
}

void compressing_horizontal(int pivot, Embedding* embedding) {
    vector<EB_Point*>* selected = embedding->get_selected();
    vector<EB_Point*>* new_selected = new vector<EB_Point*>();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    bool checkHW[topo_row][topo_column][8];
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                checkHW[i][j][k] = false;
            }
        }
    }
    
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if (x < pivot) {
            new_selected->push_back(new EB_Point(x, y, k, color));
        }
        if (x > pivot) {
            new_selected->push_back(new EB_Point(x - 1, y, k, color));
            checkHW[x - 1][y][k] = true;
        }
    }
    
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if (x == pivot) {
            if (!checkHW[x][y][k]) {
                new_selected->push_back(new EB_Point(x, y, k, color));
            }
        }
    }
    embedding->assign_selected(new_selected);
    embedding->set_topo_row(topo_row - 1);
}

void compressing_vertical(int pivot, Embedding* embedding) {
    vector<EB_Point*>* selected = embedding->get_selected();
    vector<EB_Point*>* new_selected = new vector<EB_Point*>();
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    bool checkHW[topo_row][topo_column][8];
    for (int i = 0; i < topo_row; ++i) {
        for (int j = 0; j < topo_column; ++j) {
            for (int k = 0; k < 8; ++k) {
                checkHW[i][j][k] = false;
            }
        }
    }
    
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if (y < pivot) {
            new_selected->push_back(new EB_Point(x, y, k, color));
        }
        if (y > pivot) {
            new_selected->push_back(new EB_Point(x, y - 1, k, color));
            checkHW[x][y - 1][k] = true;
        }
    }
    
    for (auto it = selected->begin(); it != selected->end(); ++it) {
        int x = (*it)->x, y = (*it)->y, k = (*it)->k, color = (*it)->color;
        if (y == pivot) {
            if (!checkHW[x][y][k]) {
                new_selected->push_back(new EB_Point(x, y, k, color));
            }
        }
    }
    embedding->assign_selected(new_selected);
    embedding->set_topo_column(topo_column - 1);
}

void compressing_total(Embedding* embedding) {
    int i = 0;
    int topo_row = embedding->get_topo_row(), topo_column = embedding->get_topo_column();
    while (i < topo_row - 1) {
        int* label_HW = embedding->generate_labelHW();
        bool check = true;
        for (int j = 0; j < topo_column; ++j){
            if (! check) {
                break;
            }
            for (int k = 0; k < 8; ++k) {
                if (!check) {
                    break;
                }
                int co = label_HW[i*topo_column*8+j*8+k];
                int co_prev = -1, co_next = -1;
                if (i >= 1) {
                    co_prev = label_HW[(i-1)*topo_column*8+j*8+k];
                }
                if (i <= topo_row - 2) {
                    co_next = label_HW[(i+1)*topo_column*8+j*8+k];
                }
                if (k <= 3) {
                    if (co != -1 && co != co_next && co != co_prev) {
                        check = false;
                    }
                }  else {
                    if (co != -1 && (co != co_next || co != co_prev)) {
                        check = false;
                    }
                }
            }
        }
        if (check) {
            compressing_horizontal(i, embedding);
            topo_row = embedding->get_topo_row();
        } else {
            i++;
        }
        delete [] label_HW;
    }
    
    int j = 0;
    topo_row = embedding->get_topo_row();
    topo_column = embedding->get_topo_column();
    while (j < topo_column - 1) {
        int* label_HW = embedding->generate_labelHW();
        bool check = true;
        for (int i = 0; i < topo_row; ++i){
            if (! check) {
                break;
            }
            for (int k = 0; k < 8; ++k) {
                if (!check) {
                    break;
                }
                int co = label_HW[i*topo_column*8+j*8+k];
                int co_prev = -1, co_next = -1;
                if (j >= 1) {
                    co_prev = label_HW[i*topo_column*8+(j-1)*8+k];
                }
                if (j <= topo_column - 2) {
                    co_next = label_HW[i*topo_column*8+(j+1)*8+k];
                }
                if (k > 3) {
                    if (co != -1 && co != co_next && co != co_prev) {
                        check = false;
                    }
                }  else {
                    if (co != -1 && (co != co_next || co != co_prev)) {
                        check = false;
                    }
                }
            }
        }
        if (check) {
            compressing_vertical(j, embedding);
            topo_column = embedding->get_topo_column();
        } else {
            j++;
        }
        delete [] label_HW;
    }
}


int main(int argc, char *argv[]) {
    srand (time(NULL));
    clock_t tStart = clock();
    auto begin = chrono::high_resolution_clock::now();
    //std::cout<< *(argv + 1);
    Graph* P = read_program(argc, argv);
    //Hardware* H = read_hardware(argc, argv);
    int topo_row = 1;
    int topo_column = 1;
    int seed_limit = 5;
    int n = P->get_n();
    Embedding* embedding = new Embedding(topo_row, topo_column);
    vector<Pass*>* list_passes = new vector<Pass*>();
    int seed_len = 5;
    vector<int>* seed_set = P->get_seed_set(seed_len);
    //cout<<seed_set->at(0);
    initialize_embedding(embedding, seed_set, seed_limit);
    extract_order(P, list_passes, seed_set, seed_limit);
    //P->print();
    embedding->expanding_border();
    auto it = list_passes->begin();
    
    while (it != list_passes->end()) {
        //cout << "NEW ROUND: " << endl;
        //embedding->print();
        /*cout << "INFO"<<endl;
        cout << (*it)->get_root()<<endl;
        for (auto it_p = (*it)->get_old_node()->begin(); it_p != (*it)->get_old_node()->end(); ++it_p) {
            cout << (*it_p) << ' ';
        }
        cout<<endl;*/
        int last_selected_size = (int)(embedding->get_selected()->size());
        //int a = a + 1;
        //cout << 1;
        expanding_chain(embedding, P, *it);
        //cout<<"EXPANDING"<<endl;
        //embedding->print();
        //bool kt = F_2(embedding, P, *it);
        int start_sz = embedding->get_selected()->size();
        if (!find_embedding_own(embedding, P, *it)) {
            expanding_total(embedding, *it, P);
        } else {
            //cout << "GAP: " <<embedding->get_selected()->size() - start_sz << endl;
            //embedding->print();
            for (auto it_p = (*it)->get_old_node()->begin(); it_p != (*it)->get_old_node()->end(); ++it_p) {
                P->reduce_degree((*it_p));
                P->reduce_degree( (*it)->get_root());
            }
            remove_trash(n, embedding, last_selected_size, P);
            compressing_total(embedding);
            embedding->expanding_border();
            it++;
        }
//break;
    }
    vector<EB_Point*> *selected_ = embedding->get_selected();
    
    int max_col = -1, max_row = -1;
    for (auto it = selected_->begin(); it != selected_->end(); it++) {
        max_row = max((*it)->x, max_row);
        max_col = max((*it)->y, max_col); 
    }
    int test_num = atoi(argv[4]);
    fstream fo;
    fo.open("Results.txt", ofstream::app);
    fo<<"Node_size"<<endl<<n<<endl;
    fo<<"Test"<<endl<<test_num<<endl;
    fo<<"Nums_qubits"<<endl<<embedding->compute_qubits_used()<<endl;
    fo<<"Max_row"<<endl<<max_row<<endl;
    fo<<"Max_column"<<endl<<max_col<<endl;
    fo<<"CPU_time"<<endl<<(double)(clock() - tStart)/CLOCKS_PER_SEC<<endl;
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    fo<<"Wall_time"<<endl<<(double)elapsed.count() * 1e-9<<endl;
    fo.close();

    
    return 0;
    //embedding->print();
}
