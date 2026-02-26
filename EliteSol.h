/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef ELITESOL_H
#define ELITESOL_H

#include <random>

class EliteSol {
public:
    EliteSol() = default;

    EliteSol(const int &capa, const int &nbNode) : capacity(capa), nbNode(nbNode) {}

    auto get_begin() {
        return this->dict.begin();
    }

    auto get_end() {
        return this->dict.end();
    }

    int get_size() {
        return static_cast<int>(this->dict.size());
    }

    auto get_sol_worst() {
        auto itr = this->get_end(); --itr; return itr;
    }

    void erase_sol_by_val(double solObj) {
        this->dict.erase(solObj);
    }

    auto get_random_member(std::mt19937 &myRnd) {
        auto randPtr = std::uniform_int_distribution<int>(0, this->get_size() - 1);
        auto it = this->get_begin();
        std::advance(it, randPtr(myRnd));
        return it;
    }

    void reset() {
        this->dict.clear();
    }

    bool set_sols(Solution &solution) {
        if (this->get_size() < this->capacity) {
            auto res = this->dict.insert({solution.SolObjcost, solution});
            return (res.second);  // return true if successfully insert
        } else {
            auto itrWorst = this->get_sol_worst();
            if (solution.SolObjcost < itrWorst->second.SolObjcost) {
                this->dict.insert({solution.SolObjcost, solution});
                if (this->get_size() > this->capacity){
                    this->dict.erase(itrWorst);
                }

                return true;
            } else {

                return false;
            }

        }

    }


private:
    int capacity;
    int nbNode;
    std::map<double, Solution> dict;
};

#endif //ELITESOL_H