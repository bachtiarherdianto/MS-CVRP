/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef BITMATRIX_H
#define BITMATRIX_H

class BitMatrix {
public:
    BitMatrix() = default;

    BitMatrix(const int &node_num_) : tabuCapacity(node_num_) {
        this->data.resize(node_num_);  // hash map
    }

    void put(const int &node) {
        const auto already_here = this->is_any(node);
        if (!already_here) {
            const int nodes_size = static_cast<int>(this->nodes.size());
            if (nodes_size >= this->tabuCapacity) {
                this->nodes.pop_back();
            }

            this->insert_to_dictionary(node);
        }

    }

    void insert_to_dictionary(const int &node)  {
        this->data[node] = true;
        this->nodes.push_front(node);
    }

    bool is_any(const int &node) {
        return this->data[node]; // tabu check
    }

    void clear_data_on(const int &node) {
        this->data[node] = false;
    }

    void clear_all_nodes() {
        this->nodes.clear();
    }

    void reset()  {
        for (int& node : this->nodes) {
            this->clear_data_on(node);
        }

        this->clear_all_nodes();
    }

private:
    std::vector<bool> data;
    std::list<int> nodes;
    int tabuCapacity;
};

#endif //BITMATRIX_H