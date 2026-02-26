/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef PAIRBITMATRIX_H
#define PAIRBITMATRIX_H

class PairBitMatrix  {
public:
    PairBitMatrix() = default;

    PairBitMatrix(const int &capacity_, const int &node_num_) {
        this->tabuCapacity = capacity_;
        this->data.resize(node_num_, node_num_);
    }

    void put(const int &node, const int &neighbor) {
        const auto already_here = this->is_any(node, neighbor);
        if (!already_here) {
            const int nodes_size = static_cast<int>(this->nodes.size());
            if (nodes_size >= this->tabuCapacity) {
                this->nodes.pop_back();
            }

            this->insert_to_dictionary(node, neighbor);
        }

    }

    void insert_to_dictionary(const int &node, const int &neighbor) {
        this->data.at(node, neighbor, true);
        this->nodes.emplace_front(node, neighbor);
    }

    bool is_any(const int &node, const int &neighbor) {
        return this->data.at(node, neighbor); // tabu check
    }

    void reset() {
        for (auto& pair : this->nodes) {
            this->data.at(pair.first, pair.second, false);
        }

        this->nodes.clear();
    }

private:
    Matrix2DFlat<bool> data;    // node x node
    std::list<std::pair<int, int>> nodes;
    int tabuCapacity;
};

#endif //PAIRBITMATRIX_H