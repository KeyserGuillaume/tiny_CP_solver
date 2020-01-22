#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>


enum status {FOUND, NO_SOLUTION, UNDECIDED, NEED_CHECKING, ABORT};
typedef std::pair<unsigned int, unsigned int> uint_pair;

// Vertex must know about ConstraintArc and vice-versa
class ConstraintArc;

class Vertex {
    // it will be the variables so they have a domain, and whether the values of the domain are deemed possible or not.
    // we assume when contructed, the domain has at least two different values.
    // we assume the domain contains unique values.
private:
    unsigned int id;
    std::vector<unsigned int> domain;
    std::vector<bool> is_possible_;
    unsigned int nb_possible_values;
    unsigned int value_index;
    std::vector<ConstraintArc*> delta;
public:
    Vertex(){}
    Vertex(const unsigned int &id, const std::vector<unsigned int> &domain): id(id), domain(domain){
        is_possible_ = std::vector<bool>(domain.size(), true);
        nb_possible_values = domain.size();
    }
    void add_arc(ConstraintArc* a){delta.push_back(a);}
    unsigned int get_id() const{return id;}
    unsigned int get_delta_size() const{return delta.size();}
    ConstraintArc* get_ith_constraint(const unsigned int& i) const{return delta[i];}
    unsigned int get_domain_size() const{return domain.size();}
    unsigned int get_ith_value(const unsigned int& i) const{return domain[i];}
    void disable_value(const unsigned int& i);
    void enable_value(const unsigned int& i);
    bool is_possible(const unsigned int &i) const{return is_possible_[i];}
    unsigned int get_nb_possible_values()const{return nb_possible_values;}
    void set_value(const unsigned int &i, std::vector<uint_pair> &diff, const bool &FC, const bool& MAC);
    unsigned int get_value() const;

    void do_forward_checking(const unsigned int &i, std::vector<uint_pair> &diff) const;
};

class ConstraintArc{
    // it will be the constraints, which are binary, so they come in different types and this is the class from which
    // they will inherit
protected:
    unsigned int id;
    Vertex *u, *v;
    ConstraintArc* symmetric;
public:
    ConstraintArc(){}
    ConstraintArc(unsigned int id, Vertex* u, Vertex* v): id(id), u(u), v(v){}
    virtual bool are_compatible(const unsigned int &i, const unsigned int &j) const = 0;
    void fill_vertex_delta(){
        u->add_arc(this);
    }
    Vertex* get_u() const {return u;}
    Vertex* get_v() const {return v;}
    unsigned int get_id()const{return id;}
    void set_symmetric(ConstraintArc* a){symmetric = a;}
    ConstraintArc* get_symmetric(){return symmetric;}
};

class DiffArc: public ConstraintArc{
public:
    DiffArc(unsigned int id, Vertex* u, Vertex* v): ConstraintArc(id, u, v){}
    bool are_compatible(const unsigned int &i, const unsigned int &j) const{
        return (u->get_ith_value(i) != v->get_ith_value(j));}
};

class QueenArc: public ConstraintArc{
public:
    QueenArc(unsigned int id, Vertex* u, Vertex* v):
            ConstraintArc(id, u, v){}
    bool are_compatible(const unsigned int &i, const unsigned int &j) const{
        unsigned int value_v = v->get_ith_value(j);
        unsigned int value_u = u->get_ith_value(i);
        unsigned int id_u = u->get_id();
        unsigned int id_v = v->get_id();
        return (value_u != value_v && value_u - value_v != id_u - id_v && value_u - value_v != id_v - id_u);}
};

class QueenSpecialArc: public QueenArc{
    // must verify the usual queen constraint, and one variable is lower than the other
public:
    QueenSpecialArc(unsigned int id, Vertex* u, Vertex* v):
            QueenArc(id, u, v){}
    bool are_compatible(const unsigned int &i, const unsigned int &j) const{
        return (QueenArc::are_compatible(i, j) && v->get_ith_value(j) < u->get_ith_value(i));}
};

class QueenSpecialArcBis: public QueenArc{
    // must verify the usual queen constraint, and one variable is lower than the other
public:
    QueenSpecialArcBis(unsigned int id, Vertex* u, Vertex* v):
            QueenArc(id, u, v){}
    bool are_compatible(const unsigned int &i, const unsigned int &j) const{
        return (QueenArc::are_compatible(i, j) && v->get_ith_value(j) > u->get_ith_value(i));}
};


void arc_consistency_routine(std::vector<ConstraintArc *> &stack, std::vector<uint_pair> &diff);


class Graph{
private:
    Vertex* V;
    std::vector<ConstraintArc*> A;
    unsigned int n_V;

    //bool detect_odd_cycle(unsigned int &u, unsigned int& v, unsigned int& w) const;
    bool detect_triangle(const std::vector<uint_pair> &edges, unsigned int &u, unsigned int& v, unsigned int& w) const;

public:
    bool debug = false;
    bool forward_checking = true;
    bool maintain_arc_consistency = true;
    unsigned int branching_strategy = 3;
    Graph(const unsigned int &n);
    Graph(const unsigned int &n, const std::vector<uint_pair> &edges, const unsigned int &K);
    ~Graph(){
        delete [] V;
        for (unsigned int i = 0; i < A.size(); i++)
            delete [] A[i];
    }
    bool check_solution() const;
    status solve(const clock_t &time_limit, unsigned int& nb_nodes);
    std::vector<unsigned int> get_solution() const;
    void make_arc_consistent(std::vector<uint_pair>& diff);
    void make_arc_consistent();
    void branch_strat_1(status &s, Vertex *&current_var) const;
    void branch_strat_2(status &s, Vertex *&current_var) const;
    void branch_strat_3(status &s, Vertex *&current_var) const;
};
