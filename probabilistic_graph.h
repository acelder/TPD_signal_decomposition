//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H
#define TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H

#include "gaussian_operations.h"
#include "factor_operations.h"
#include <vector>
#include <algorithm>
#include <list>
#include <queue>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;

namespace gml {  // "graphical model library"

    class clique_tree;
    typedef vector<variable<int>> var_vec;
    typedef vector<variable<int>*> pvar_vec;

    class probabilistic_graph {
        vector<variable<int>> discrete_variables;

        vector<int> elimination_order;


        void add_variable(const variable<int> &v);


        // These functions irreversibly alter the graph. They should only be called as part of the
        // ve_query_probability function.
        void sum_product_eliminate_var(variable<int> *Z);
        void sum_product_ve();

    public:
        vector<factor> factors;
        unordered_map<string, variable<int> *> variable_map;
        unordered_map<string, factor*> factor_map;


        probabilistic_graph() = default;

        probabilistic_graph(const probabilistic_graph &to_copy);

        // destructors
        ~probabilistic_graph() = default;

        // manipulation functions
        void add_factor(const variable<int> &v);
        void add_factor(const variable<int> &v, const var_vec &parents);

        // these constructors make use of the above functions
        explicit probabilistic_graph(const vector<variable<int>> &variable_list);
        explicit probabilistic_graph(const vector<pair<variable<int>, var_vec>> &variable_pairs);

        factor ve_query_probability( vector<int> elimination_order, bool normalize = true);

        void observe_evidence(variable<int> *target_variable, int var_assign);

        void compute_elimination_order();
        void input_elimination_order(vector<int> order);

        // misc
        void print_graph();


    };

    struct clique_edge;
    struct clique {
        vector<factor> members;
        vector<clique_edge *> edges;
        bool ready_to_transmit = false;
        factor belief;
    };

    struct clique_edge {
        pair<clique *, clique *> cliques;
        factor belief;

        clique *adjacent(clique *current) {
            if (current == cliques.first) {
                return cliques.second;
            } else if (current == cliques.second) {
                return cliques.first;
            } else {
                return NULL;
            }
        }
    };

    class clique_tree {
    private:
        clique *root;
        vector<clique> cliques;
        vector<clique_edge> edges;
        vector<vector<clique *>> sepsets;

        vector<factor> clique_beliefs;
        vector<factor> sepset_beliefs;

        void initialize_message_schedule();

        void initialize_cliques();

        pair<clique *, clique *> next_ready_clique();

        factor compute_message(const clique &clique_i, const clique &clique_j);

        void pass_message();

    public:
        const double default_onvergence_criterion = 0.01;

        clique_tree() = default;

        clique_tree(const vector<pair<unordered_set<factor *>, vector<int>>> &clique_list);

        void calibrate();
    };

    /**
     * Compute the maximum spanning tree of the graph using Kruskal's algorithm
     * @param adjacency_list : edge list in for form of a vector of pairs, the first entry is the weight, the second entry is a pair of nodes
     * @return : adjacency list of maximum spanning tree
     */

    vector<pair<int, pair<int, int>>> maximum_spanning_tree(vector<pair<int, pair<int, int>>> edge_list, int V) {
        auto cmp = [](pair<int, pair<int, int>> a, pair<int, pair<int, int>> b) {
            return a.first > b.first;
        };
        typedef pair<int, pair<int, int>> adj_edge;
        priority_queue<adj_edge, vector<adj_edge>, decltype(cmp)> q(cmp);
        for (auto edge: edge_list) {
            q.push(edge);
        }

        adj_edge e = q.top();
        q.pop();
        vector<adj_edge> T;
        unordered_set<int> X;

        T.push_back(e);
        X.insert(e.second.first);
        X.insert(e.second.second);
        while (T.size() != V - 1) {
            e = q.top();
            q.pop();
            auto it_a = X.find(e.second.first);
            auto it_b = X.find(e.second.first);
            if (it_a == X.end() or it_b == X.end()) {
                T.push_back(e);
            }
            if (q.empty()) {
                cout << "Warning: graph is disconnected.\n";
            }
        }
        return T;
    }


    probabilistic_graph::probabilistic_graph(const probabilistic_graph &to_copy) {
        for (auto v: to_copy.discrete_variables) {
            add_factor(v);
        }
    }

    void probabilistic_graph::add_factor(const variable<int> &v) {
        factor phi(v);
        factors.push_back(phi);
        factor_map[phi.name] = &factors.back();
        phi.canonical_table.resize(v.cardinality);
        add_variable(v);
    }

    void probabilistic_graph::add_factor(const variable<int> &v, const var_vec &parents){
        // add child variable, must be added after each parent has been added
        factor phi(v);
        add_variable(v);
        factors.push_back(phi);
        factor_map[phi.name] = &factors.back();


        int cardinality = v.cardinality;

        for(auto parent: parents){
            // add pointers between parents and children
            factors.back().parents.push_back(factor_map[parent.name]);
            factors.back().discrete_scope.push_back(parent);

            // compute product of cardinalities
            cardinality *= parent.cardinality;
        }

        factors.back().canonical_table.resize(cardinality);
    }

    void probabilistic_graph::add_variable(const variable<int> &v) {
        discrete_variables.push_back(v);
        variable_map[v.name] = &discrete_variables.back();
    }


    probabilistic_graph::probabilistic_graph(const vector<variable<int>> &variable_list) {
        // make copies of the input factors
        for (auto v: variable_list) {
            add_factor(v);
        }
    }

    probabilistic_graph::probabilistic_graph(const vector<pair<variable<int>, var_vec>> &variable_pairs){
        // assumes variable_pairs are in order
        for(auto p: variable_pairs){
            add_factor(p.first, p.second);
        }
        // add pointers from children to parents
        for(auto phi: factors){
            for(auto c: phi.children){
                c->parents.push_back(&phi);
            }
        }
    }

    factor probabilistic_graph::ve_query_probability( vector<int> elim_order, bool normalize){
        elimination_order = elim_order;
        sum_product_ve();
        factor phi = factors[0];
        // multiply remaining factors together
        while(factors.size() > 1){
            phi = factor_product(phi, factors.back(), "phi");
            factors.pop_back();
        }

        if(normalize){
            double t = 0.0;
            for(int i=0; i<phi.canonical_table.size(); ++i){
                t += exp(phi.canonical_table[i].g);
            }
            for(int i=0; i<phi.canonical_table.size(); ++i){
                phi.canonical_table[i].g -= log(t);
            }
        }

        return phi;;
    }

    /**
     * Compute the elimination ordering using the maximum cardinality search algorithm.
     */
//
//    void probabilistic_graph::compute_elimination_order() {
//        int n_factors = factors.size();
//        elimination_order.resize(n_factors);
//        unordered_map<factor*, bool> mark;
//        for(auto &f: factors){
//            mark[&f] = false;
//        }
//
//        auto marked_neighbors = [](factor* phi, unordered_map<factor*, bool> m){
//            int n_marked_neighbors;
//            for(auto c: phi->children){ if(m[c]){++n_marked_neighbors;}}
//            for(auto p: phi->parents){ if(m[p]){++n_marked_neighbors;}}
//            return n_marked_neighbors;
//        };
//
//
//        auto cmp = [](factor* left, factor* right){
//            return marked_neighbors(left) > marked_neighbors(right);
//        };
//
//        vector<factor*> factor_heap;
//        for(auto phi: factors){
//            factor_heap.push_back(&phi);
//        }
//        make_heap(factor_heap.begin(), factor_heap.end(), cmp);
//
//        int k = n_factors - 1;
//        while(!factor_heap.empty()){
//            pop_heap(factor_heap.begin(), factor_heap.end(), cmp);
//            factor* largest = factor_heap.back();
//            factor_heap.pop_back();
//
//            int pos = distance(largest, &factors[0]);
//            elimination_order[pos] = k;
//            --k;
//            mark[largest] = true;
//        }
//    }

    void probabilistic_graph::input_elimination_order(vector<int> order) {
        elimination_order = order;
    }

//     Z should be a pointer to a variable stored in the *graph* (not a copy stored by a factor)
    void probabilistic_graph::sum_product_eliminate_var(variable<int> *Z) {
        // determnine which factors will remain and which will be multiplied
        vector<factor> factors_to_multiply;
        vector<factor> factors_to_remain;
        string name = Z->name;
        for (auto &phi: factors) {
            bool variable_in_scope = false;
            for (int i = 0; i < phi.discrete_scope.size(); ++i) {
                if (!name.compare(phi.discrete_scope[i].name)) { // not sure why CLion doesn't like that
                    variable_in_scope = true;
                }
            }

            if (variable_in_scope) {
                factors_to_multiply.push_back(phi);
            } else {
                factors_to_remain.push_back(phi);
            }
        }

        // compute the product of the factors with Z in their scope
        factor psi = unit_factor();
        for (auto &phi: factors_to_multiply) {
            psi = factor_product(psi, phi, "psi_" + Z->name);
        }
        // get a pointer to the variable to be marginalized (the sum step)
        variable<int> *to_marginalize;
        for (auto &v: psi.discrete_scope) {
            if (!Z->name.compare(v.name)) {
                to_marginalize = &v;
            }
        }

        factor tau = discrete_marginal(psi, to_marginalize, "tau_" + Z->name);

        // Replace the factor list with the remaining factors as well as the new one.
        factors_to_remain.push_back(tau);

        factors = factors_to_remain;
        string message = "";

        for(auto f: factors_to_multiply){
            message += f.name + ", ";
        }
        cout << "Multiplied factors: " << message << endl;
    }

    void probabilistic_graph::sum_product_ve() {

        for( auto o: elimination_order ){
            sum_product_eliminate_var(&discrete_variables[o]);
        }

    }

    void probabilistic_graph::print_graph(){
        for(auto phi: factors){
            phi.print();
        }
    }

//    void probabilistic_graph::construct_clique_tree(const string& algorithm) {
//        if( !algorithm.compare("Bron–Kerbosch") ){
//
//            // begin by making sure the variable elimination ordering has been computed
//            if(!elimination_order.size()){
//                compute_elimination_order();
//            }
//
//            list<pair<set<variable*>, variable*>> factor_scopes;
//            for(int i=0; i<factors.size(); ++i){
//                factor_scopes.push_back(make_pair(factors[i], factors[i].total_scope()));
//            }
//
//            vector<vector<variable*>> intermediate_factor_scopes;
//            // construct the neighbor map representation of the induced graph-begin by simulating variable elimination
//            for(int i=0; i<factors.size(); ++i){
//                // collect all factors such that elimination_order[i] is in their scope
//                vector<variable*> intermediate_factor_scope;
//                variable* to_elim = ordered_variables[elimination_order[i]];
//
//
//                for(auto scope: factor_scopes){
//                    if(scope.first.find(to_elim) != scope.first.end()){
//                        intermediate_factor_scope.push_back(scope.second);
//                    }
//                }
//                intermediate_factor_scopes.push_back(intermediate_factor_scope);
//            }
//
//            // Variables that end up in shared intermediate factors define edges. Bron-Kerbosch is not needed.
//
//            int n_cliques = intermediate_factor_scopes.size();
//            vector<unordered_set<vector*>> maximal_cliques (n_cliques);
//            for(int i=0; i<n_cliques; ++i){
//                for(int j=0; j<intermediate_factor_scopes[i].size(); ++j){
//                    maximal_cliques[i].insert(intermediate_factor_scopes[i][j]);
//                }
//            }
//
//
//            // To apply Kruskal's algorithm, we'll store the graph as an edge list where each edge weight is the size
//            // of the intersection of the two cliques. We'll switch over to index representation. Recall that the
//            // variables and their primitive (vocab?) factors are in the same order in their respective probabilistic_graph
//            // member containers.
//            typedef pair<int, pair<int, int>> edge;
//            vector<edge> edge_list;
//            set<pair<int, int>> edge_check;
//            int i = 0;
//            for(auto &psi_1: maximal_cliques){
//                int j = 0;
//                for(auto &psi_2: maximal_cliques){
//                    // if an edge already exists, move on
//                    pair<int, int> nodes = make_pair(i, j);
//                    pair<int, int> rev_nodes = make_pair(j, i);
//
//                    if(edge_check.find(nodes) == edge_check.end() and edge_check.find(rev_nodes) == edge_check.end()){
//
//                        // compute intersection size
//                        int weight = 0;
//                        for(auto member: psi_1){
//                            if(psi_2.find(member) != psi_2.end()){
//                                weight++; // add 1 to the weight for each shared factor
//                            }
//                        }
//
//                        edge e = make_pair(weight, nodes);
//                        edge_list.push_back(e);
//                    }
//                    j++;
//                }
//                i++;
//            }
//
//            // apply Kruskal's algorithm to find the maximum spanning tree
//            vector<edge> MST = maximum_spanning_tree(edge_list, n_cliques);
//
//            //*********************************************************************************************
//            //         Apply pruning routine here- remove leaves that are subsets of their parents
//            //*********************************************************************************************
//
//            // construct the adjacency list to pass as an argument to the clique tree constructor
//            typedef vector<pair<unordered_set<factor*>, vector<int>>> clique_adjacency;
//            clique_adjacency clique_list (n_cliques);
//            i = 0;
//            for(auto c: maximal_cliques){
//                clique_list[i].first = c;
//                ++i;
//            }
//            for(auto e: MST){
//                int node_1 = e.second.first;
//                int node_2 = e.second.second;
//                clique_list[node_1].second.push_back(node_2);
//                clique_list[node_2].second.push_back(node_1);
//            }
//
//            clique_tree* tree = new clique_tree(clique_list);
//            c_tree = *tree;// this should prevent the new clique tree from going out of scope
//        }
//    }



//    clique_tree::clique_tree(const vector<pair<unordered_set<factor*>, vector<int>>>& clique_list){
//        // Construct clique tree from factors, adjacencies are stored in the second element of the list of pairs.
//        // Note: each node in the clique tree is constructed as the product of the cliques in the parent graph.
//
//        size_t n_cliques = clique_list.size();
//        cliques.resize(n_cliques);
//        set<pair<int, int>> edge_check;
//        for(int i=0; i<n_cliques; ++i){
//            for(auto phi: clique_list[i].first){
//                cliques[i].members.push_back(*phi);
//            }
//            for(auto c: clique_list[i].second){
//                auto check_rev = edge_check.find(make_pair(c, i)); // check that the reverse edge has not already been entered
//                if(check_rev == edge_check.end()){
//                    clique_edge e;
//                    e.cliques.first = &cliques[i];
//                    e.cliques.second = &cliques[c];
//                    edge_check.insert(make_pair(i, c));// note that the edge has been stored
//                    edges.push_back(e);
//                }
//            }
//        }
//
//        // set the root as the highest cardinality clique
//        root = &cliques.front();
//        for(auto c: cliques){
//            if(c.edges.size() > root->edges.size()){
//                root = &c;
//            }
//        }
//
//    }
//
//    void clique_tree::initialize_message_schedule(){
//        // mark leaves as ready to transmit
//        for(auto c: cliques){
//            if(c.edges.size() == 1){
//                c.ready_to_transmit = true;
//            }
//            else{
//                c.ready_to_transmit = false;
//            }
//        }
//    }





}

#endif //TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H
