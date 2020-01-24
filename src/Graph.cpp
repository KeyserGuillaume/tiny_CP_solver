#include <form.h>
#include "Graph.h"


void Vertex::enable_value(const unsigned int &i) {
    if (is_possible_[i])
        throw std::logic_error("Value " + std::to_string(domain[i]) + " is already possible");
    is_possible_[i] = true;
    nb_possible_values++;
}

void Vertex::disable_value(const unsigned int &i) {
    if (!is_possible_[i])
        throw std::logic_error("Value " + std::to_string(domain[i]) + " is already impossible");
    is_possible_[i] = false;
    nb_possible_values--;
    if (nb_possible_values == 1){
        for (value_index = 0; !is_possible_[value_index]; value_index++){
            if (value_index == domain.size())
                throw std::logic_error("only one possible value left but they are all impossible");
        }
    }
}

void Vertex::set_value(const unsigned int &i, std::vector<uint_pair> &diff, const bool &FC, const bool &MAC) {
    if (!is_possible_[i])
        throw std::logic_error("Value " + std::to_string(domain[i]) + " is impossible");
    // reduce possible domain to the chosen value
    for (unsigned int j = 0; j < domain.size(); j++){
        if (is_possible_[j] && j != i){
            disable_value(j);
            diff.push_back(uint_pair(id, j));
        }
    }
    // forward_checking
    if (FC && !MAC)
        do_forward_checking(i, diff);
    // maintain_arc_consistency
    if (MAC){
        std::vector<ConstraintArc*> stack(0);
        for (unsigned int s = 0; s < delta.size(); s++)
            stack.push_back(delta[s]->get_symmetric());
        arc_consistency_routine(stack, diff);
    }
}

void Vertex::do_forward_checking(const unsigned int &i, std::vector<uint_pair> &diff) const {
    for (unsigned int m = 0; m < delta.size(); m++){
        ConstraintArc* a = delta[m];
        Vertex* v = a->get_v();
        for (unsigned int s = 0; s < v->domain.size(); s++){
            if (v->is_possible_[s] && !a->are_compatible(i, s)){
                v->disable_value(s);
                diff.push_back(uint_pair(v->id, s));
            }
        }
    }
}

unsigned int Vertex::get_value() const {
    if (nb_possible_values > 1)
        throw std::logic_error("there are more than 1 values");
    return domain[value_index];
}


Graph::Graph(const unsigned int &n) {
    // n-queens
    // we add one small constraint : the last queen must be below the first one.
    n_V = n;
    V = new Vertex [n_V];
    A = std::vector<ConstraintArc*> (0);
    unsigned int next_arc_index = 0;

    std::vector<unsigned int> domain(n, 0);
    for (unsigned int i = 0; i < n; i++)
        domain[i] = i;
    for (unsigned int i = 0; i < n; i++)
        V[i] = Vertex(i, domain);

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < i; j++) {
            if (i == n - 1 && j == 0 && avoid_symmetries){
                A.push_back(new QueenSpecialArc(next_arc_index, V + i, V + j));
                next_arc_index++;
                A.push_back(new QueenSpecialArcBis(next_arc_index, V + j, V + i));
                A[next_arc_index]->set_symmetric(A[next_arc_index - 1]);
                A[next_arc_index - 1]->set_symmetric(A[next_arc_index]);
                next_arc_index++;
            } else {
                A.push_back(new QueenArc(next_arc_index, V + i, V + j));
                next_arc_index++;
                A.push_back(new QueenArc(next_arc_index, V + j, V + i));
                A[next_arc_index]->set_symmetric(A[next_arc_index - 1]);
                A[next_arc_index - 1]->set_symmetric(A[next_arc_index]);
                next_arc_index++;
            }
        }
    }

    for (unsigned int i = 0; i < A.size(); i++){
        A[i]->fill_vertex_delta();
    }
}

Graph::Graph(const unsigned int &n, const std::vector<uint_pair> &edges, const unsigned int &K) {
    // graph coloring with at most K colors
    // we assume that each edge is gven twice, like in the file "anna.col".
    // we assume absence of duplicate edges...
    // we choose the first given edge and impose two different colors for the nodes.
    // we should at least look for triangles though.
    n_V = n;
    V = new Vertex [n_V];
    A = std::vector<ConstraintArc*> (0);
    unsigned int next_arc_index = 0;

    std::vector<unsigned int> domain(K, 0);
    for (unsigned int i = 0; i < K; i++)
        domain[i] = i;
    for (unsigned int i = 0; i < n; i++)
        V[i] = Vertex(i, domain);

    for (unsigned int i = 0; i < edges.size(); i++) {
        if (edges[i].first > edges[i].second) continue;
        A.push_back(new DiffArc(next_arc_index, V + edges[i].first, V + edges[i].second));
        next_arc_index++;
        A.push_back(new DiffArc(next_arc_index, V + edges[i].second, V + edges[i].first));
        A[next_arc_index]->set_symmetric(A[next_arc_index - 1]);
        A[next_arc_index - 1]->set_symmetric(A[next_arc_index]);
        next_arc_index++;
    }

    for (unsigned int i = 0; i < A.size(); i++){
        A[i]->fill_vertex_delta();
    }

    if (!avoid_symmetries) return;

    std::vector<uint_pair> diff(0);
    unsigned int epsilon;
    unsigned int u, v, w;
    if (detect_triangle(edges, u, v, w)){
        V[u].set_value(0, diff, forward_checking, maintain_arc_consistency);
        V[v].set_value(1, diff, forward_checking, maintain_arc_consistency);
        V[w].set_value(2, diff, forward_checking, maintain_arc_consistency);
        epsilon = 3;
    } else if (edges.size() > 0) {
        V[edges[0].first].set_value(0, diff, forward_checking, maintain_arc_consistency);
        V[edges[0].second].set_value(1, diff, forward_checking, maintain_arc_consistency);
        epsilon = 2;
    } else { // this means the pb is trivial...
        epsilon = 0;
    }

    // epsilon is the nb of preset colors (2 or 3)
    // Select K - epsilon variables x_i and impose x_i <= epsilon + i. If a solution exists, we can shuffle the colors to
    // verify this condition, because there are at most K colors, and the following algorithm works:
    // for i in 1 ... K - epsilon do:
    //     we denote l_i the color of x_i. if l_i > epsilon + i:
    //         exchange l_i with color epsilon + i for all the variables. some variables take color l_i but none of the x_j
    //         with j < i because those have a color lower than epsilon + j and epsilon + j < epsilon + i.
    // the x_i are chosen among variables with the highest number of possible values left in the domain.
    for (unsigned int i = 0; i < K - epsilon; i++){
        Vertex* u = V;
//        unsigned int max = 0;
//        while(u != V + n){
//            if (u->get_nb_possible_values() > max)
//                max = u->get_nb_possible_values();
//            u++;
//        }
        u = V + i;
//        while (u->get_nb_possible_values() != max)
//            u++;
        for (unsigned int j = epsilon + i; j < K; j++)
            if (u->is_possible(j))
                u->disable_value(j);
    }
}


status Graph::solve(const clock_t &time_limit, unsigned int& nb_nodes) {
    if (nb_nodes == 0 && maintain_arc_consistency)
        make_arc_consistent();
    nb_nodes++;
    if (clock() > time_limit)
        return ABORT;
    // choose the variable on which we will do the branching
    status s = UNDECIDED;
    Vertex *current_var;
    if (branching_strategy == 1)
        branch_strat_1(s, current_var);
    else if (branching_strategy == 2)
        branch_strat_2(s, current_var);
    else
        branch_strat_3(s, current_var);

    if (s == NEED_CHECKING){
        if (check_solution())
            return FOUND;
        else
            return NO_SOLUTION;
    } else if (s == NO_SOLUTION)
        return NO_SOLUTION;

    // branch
    for (unsigned int i = 0; i < current_var->get_domain_size(); i++){
        if (!current_var->is_possible(i)) continue;
        std::vector<uint_pair> diff(0);
        current_var->set_value(i, diff, forward_checking, maintain_arc_consistency);
        status s = solve(time_limit, nb_nodes);
        if (s != NO_SOLUTION)
            return s;
        // revert changes made
        for (unsigned int m = 0; m < diff.size(); m++)
            V[diff[m].first].enable_value(diff[m].second);
    }

    // tell the previous call that no solution has been found
    return NO_SOLUTION;
}

void Graph::branch_strat_1(status &s, Vertex *&current_var) const {
    // branch on first variable we can find with at least 2 possible values left
    current_var = V;
    while (current_var->get_nb_possible_values() <= 1 && current_var != V + n_V) {
        if (current_var->get_nb_possible_values() == 0) {
            s = NO_SOLUTION;
            return;
        }
        current_var++;
    }
    if (current_var == V + n_V)
            s = NEED_CHECKING;
    // check all variables have at least one possible value otherwise prune search tree
    for (Vertex *u = current_var; u != V + n_V; u++){
        if (u->get_nb_possible_values() == 0) {
            s = NO_SOLUTION;
            return;
        }
    }
}

void Graph::branch_strat_2(status &s, Vertex *&current_var) const {
    // branch on a variable with as few possible values left as possible (what we would instinctively do on sudoku)
    unsigned int min = V->get_domain_size() + 1;
    current_var = V;
    for (unsigned int i = 0; i < n_V; i++){
        unsigned int nb_val = V[i].get_nb_possible_values();
        if (nb_val == 0) {
            s = NO_SOLUTION;
            return;
        }
        else if (nb_val > 1 && nb_val < min){
            min = nb_val;
            current_var = V + i;
        }
    }
    if (min == V->get_domain_size() + 1)
            s = NEED_CHECKING;
}

void Graph::branch_strat_3(status &s, Vertex *&current_var) const {
    // branch on a variable with as few possible values left as possible (what we would instinctively do on sudoku)
    // and (secondary objective) as many neighbors as possible (what we would do in graph coloring)
    unsigned int min_nb_val = V->get_domain_size() + 1;
    unsigned int max_nb_nei = 0;
    current_var = V;
    for (unsigned int i = 0; i < n_V; i++){
            unsigned int nb_val = V[i].get_nb_possible_values();
            if (nb_val == 0) {
                s = NO_SOLUTION;
                return;
            }
            else if (nb_val == 1)
                continue;
            else if (nb_val < min_nb_val){
                min_nb_val = nb_val;
                current_var = V + i;
                max_nb_nei = V[i].get_delta_size();
            }
            else if (nb_val == min_nb_val && V[i].get_delta_size() > max_nb_nei){
                max_nb_nei = V[i].get_delta_size();
                current_var = V + i;
            }
        }
    if (min_nb_val == V->get_domain_size() + 1)
            s = NEED_CHECKING;
}

std::vector<unsigned int> Graph::get_solution() const {
    std::vector<unsigned int> solution(0);
    for (unsigned int i = 0; i < n_V; i++){
        solution.push_back(V[i].get_value());
    }
    return solution;
}

bool Graph::check_solution() const {
    // we check all the constraints but in fact they are identical by pairs so it would also
    // be possible to check one for each pair.
    for (unsigned int i = 0; i < A.size(); i++){
        if (!A[i]->are_compatible(A[i]->get_u()->get_value(), A[i]->get_v()->get_value()))
            return false;
    }
    return true;
}

bool Graph::detect_triangle(const std::vector<uint_pair> &edges, unsigned int &u, unsigned int &v, unsigned int &w) const {
    // intended to limit symmetry of graph coloring
    // very inefficient but no time right now
    for (unsigned int i = 0; i < edges.size(); i++){
        u = edges[i].first;
        v = edges[i].second;
        Vertex* var_u = V + u;
        Vertex* var_v = V + v;
        for (unsigned int s = 0; s < var_u->get_delta_size(); s++){
            for (unsigned int m = 0; m < var_v->get_delta_size(); m++) {
                if (var_u->get_ith_constraint(s)->get_v()->get_id() ==
                    var_v->get_ith_constraint(m)->get_v()->get_id()) {
                    w = var_u->get_ith_constraint(s)->get_v()->get_id();
                    return true;
                }
            }
        }
    }
    return false;
}

void arc_consistency_routine(std::vector<ConstraintArc *> &stack, std::vector<uint_pair> &diff) {
    while (stack.size() > 0){
        ConstraintArc* a = stack[stack.size() - 1];
        Vertex* u = a->get_u();
        Vertex* v = a->get_v();
        stack.pop_back();
        for (unsigned int m = 0; m < u->get_domain_size(); m++){
            if (!u->is_possible(m)) continue;
            bool found = false;
            for (unsigned int s = 0; s < v->get_domain_size() && !found; s++)
                if (v->is_possible(s) && a->are_compatible(m, s))
                    found = true;
            if (!found){
                u->disable_value(m);
                diff.push_back(uint_pair(u->get_id(), m));
                for (unsigned int t = 0; t < u->get_delta_size(); t++)
                    if (u->get_ith_constraint(t)->get_id() != a->get_id())
                        stack.push_back(u->get_ith_constraint(t)->get_symmetric());
            }
        }
    }
}

void Graph::make_arc_consistent(std::vector<uint_pair> &diff) {
    std::vector<ConstraintArc*> stack(0);
    for (unsigned int i = 0; i < A.size(); i++)
        stack.push_back(A[i]);
    arc_consistency_routine(stack, diff);
}

void Graph::make_arc_consistent() {
    std::vector<uint_pair> diff(0);
    make_arc_consistent(diff);
}

//bool Graph::detect_odd_cycle(unsigned int &u, unsigned int &v, unsigned int &w) const {
//    // intended to limit symmetry of graph coloring
//    it is wrong though, only cliques can impose colors.
//    std::vector<unsigned int> parity(n_V, 0);
//    std::vector<unsigned int> predecessor(n_V, 0);
//    std::vector<unsigned int> stack(1, 0);
//    parity[0] = 1;
//    while (stack.size() > 0){
//        unsigned int curr_id = stack[stack.size() - 1];
//        Vertex* current = V + curr_id;
//        stack.pop_back();
//        for (unsigned int m = 0; m < current->get_delta_size(); m++){
//            unsigned int neighbor_id = current->get_ith_constraint(m)->get_v()->get_id();
//            unsigned int s = parity[neighbor_id];
//            if (parity[curr_id] == s){
//                u = predecessor[curr_id];
//                v = curr_id;
//                w = neighbor_id;
//                return true;
//            } else if (s == 0) {
//                predecessor[neighbor_id] = curr_id;
//                parity[neighbor_id] = parity[curr_id] == 1 ? 2 : 1;
//                stack.push_back(neighbor_id);
//            }
//        }
//    }
//    return false;
//}
