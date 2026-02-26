/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef EVALUATIONMEMORY_H
#define EVALUATIONMEMORY_H

class EvaluationMemory {
public:
    int node;
    int neighbor;
    int way;
    int longPathMove;
    double addInsert;
    double addRemove;
    double result;
    double operaTime;
    double charge;
    int category;

    EvaluationMemory() {
        this->node = this->neighbor = this->way = this->longPathMove = this->category = -1u;
        this->addRemove = this->charge = this->operaTime = this->addInsert = this->result = -1.0f;
    }

    ~EvaluationMemory() = default;

    [[nodiscard]] bool is_feasible(const double &Q, const int &L) const {
        return (this->charge <= Q && this->operaTime <= L);
    }

    friend std::ostream &operator<<(std::ostream &os, EvaluationMemory &A) {
        os << "[node: " << A.node << " neighbor " << A.neighbor << " result: " << A.result << " addInsert: " << A.addInsert << " addRemove: " << A.addRemove << " Way: " << A.way << " longPathMove: " << A.longPathMove << " charge: " << A.charge << " operaTime: " << A.operaTime << "] ";
        return os;
    }
};

#endif //EVALUATIONMEMORY_H