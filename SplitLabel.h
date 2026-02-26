/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef SPLITLABEL_H
#define SPLITLABEL_H

struct SplitLabel {
    int index;                  // label position
    double cost = my::dummy;    // cost to get to the label
    int C;                      // category of the last tour
    double lengthTrips;         // cost of the last trips
    std::vector<int> vhlUsed;   // no of vehicle used
    SplitLabel *previousLabel;  // pointer for predecessor label

    SplitLabel() = default;

    SplitLabel(const SplitLabel &sl) {
        this->index = sl.index;
        this->cost = sl.cost;
        this->C = sl.C;
        this->lengthTrips = sl.lengthTrips;
        this->vhlUsed = sl.vhlUsed;
        this->previousLabel = sl.previousLabel;
    }

    ~SplitLabel() {
        this->previousLabel = NULL;
        delete previousLabel;
    }

};

#endif //SPLITLABEL_H