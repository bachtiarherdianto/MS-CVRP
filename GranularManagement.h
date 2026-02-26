/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef GRANULARMANAGEMENT_H
#define GRANULARMANAGEMENT_H

class GranularManagement {
public:
    GranularManagement() = default;

    GranularManagement(const int &node_num, const int &size) {
        this->data.resize(node_num);
        this->nNeighbor.resize(node_num);
        this->marks.resize(node_num);
        this->maxSize = size;
    }

    void put(const int &node) {
        const auto already_here = this->is_any(node);
        if (!already_here) {
            this->insert_to_dictionary(node);
        }
    }

    void insert_to_dictionary(const int &node) {
        this->data[node] = true;
        this->nodes.push_back(node);
    }

    bool is_any(const int &node) { return this->data[node]; }

    int get_neighbor(const int &node) { return this->nNeighbor[node]; }

    void set_growth(const int &node, const int growth) { this->nNeighbor[node] += growth; }

    void clear_data_on(const int &node) {
        this->data[node] = false;
    }

    void clear_all_nodes() {
        this->nodes.clear();
    }

    void initialization(const int &size_init) {
        for (int& neighbor : this->nNeighbor) {
            neighbor = size_init;
        }
    }

    void granular_reset(const int &size_init) {
        for (int& node : changed) {
            this->nNeighbor[node] = size_init;
            this->marks[node] = false;
        }

        this->changed.clear();
    }

    void growth_neighbors(const int &growth) {
        for (int& node : this->nodes) {
            if (this->get_neighbor(node) + growth <= this->maxSize) {
                this->set_growth(node, growth);

                if (this->marks[node] != true) {
                    this->changed.push_back(node);
                    this->marks[node] = true;
                }

                this->clear_data_on(node);
            }
        }

        this->clear_all_nodes();
    }

private:
    std::vector<bool> data;
    std::vector<int> nodes;
    std::vector<int> nNeighbor;
    std::vector<int> marks;
    std::vector<int> changed;
    int maxSize;
};

#endif //GRANULARMANAGEMENT_H