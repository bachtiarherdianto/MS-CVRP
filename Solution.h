/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef SOLUTION_H
#define SOLUTION_H

#include <iomanip>

#include "Welford.h"

class Solution {
public:
  std::vector<int> Category;
  std::vector<std::vector<int>> Soltours;
  std::vector<std::vector<double>> Costours;
  std::vector<std::vector<double>> Charge;
  std::vector<double> Operational;
  std::vector<int> routesNode;
  std::vector<int> positionNode;

  std::vector<CircleSector> routeSector;

  int nbNodes;

  double SolObjcost;

  Solution() {
    this->SolObjcost = std::numeric_limits<float>::max();
  }

  Solution(const Solution &B) {
    copy(B);
  }

  ~Solution() = default;

  Solution &operator=(const Solution &source) {
    if (this == &source) {
      return *this;
    }

    copy(source);
    return *this;
  }

  bool operator==(const Solution &other) const {
    if (std::fabs(this->SolObjcost - other.SolObjcost) >= EPSILON) {
      return false;
    }

    for (int i = 1; i < this->nbNodes; ++i) {   // 0 is depot
      if (get_node_prev(i) != other.get_node_prev(i) || get_node_next(i) != other.get_node_next(i)) {
        return false;
      }

    }

    return true;
  }

  bool operator!=(const Solution &other) const {
    return !(*this == other);
  }

  friend std::ostream &operator<<(std::ostream &os, Solution &A);

  inline int get_size() const { return static_cast<int>(this->Soltours.size()); }

  inline int get_operational_size() const { return this->Operational.size(); }

  inline int get_route_size(const int route) { return static_cast<int>(this->Soltours[route].size()); }

  inline int get_route_index(const int node) { return this->routesNode[node]; }

  inline int get_position_node(const int node) { return this->positionNode[node]; }

  inline double get_route_cost_at(const int &route, const int node) { return this->Costours[route][node]; }

  inline double get_route_cost(const int &route) { return this->Costours[route][0]; }

  inline double get_route_load(const int &route) { return this->Charge[route][0]; }

  inline double get_route_time(const int &route) { return this->Operational[route]; }

  inline bool feasible_tour(const int &tour, const double &Q, const int &L) { return (this->get_route_load(tour) <= Q && this->get_route_time(tour) <= L); }  // indicates whether the tour respects the (potentially relaxed) capacity and duration constraints

  inline double capacity_utilization(const int &tour, const double &Q) { return (this->get_route_load(tour) / Q); }  // calculate capacity utilization

  inline std::vector<int> get_route(const int &i) { return this->Soltours[i]; }

  inline std::vector<double> get_route_cost_list(const int &i) { return this->Costours[i]; }

  inline std::vector<double> get_route_load_list(const int &i) { return this->Charge[i]; }

  [[nodiscard]] inline int get_node_next(const int &node) const {
    assert(node != my::dummy);
    return this->Soltours[this->routesNode[node]][this->positionNode[node]+1];
  }

  [[nodiscard]] inline int get_node_prev(const int &node) const {
    const auto routeIndex = this->routesNode[node];
    const auto pos = this->positionNode[node];
    return this->Soltours[routeIndex][pos-1];
  }

  inline void set_route(const int i, std::vector<int> route) { this->Soltours[i] = route; }

  inline int node_at(const int indexRoute, const int indexPos) {
    return this->Soltours[indexRoute][indexPos];
  }

  double get_s19(const double& Q) {
    Welford sol;
    for (int i = 0; i < this->get_size(); ++i) {
      sol.update(this->capacity_utilization(i, Q));
    }

    return sol.get_mean();
  }

  double get_s20(const double& Q) {
    Welford sol;
    for (int i = 0; i < this->get_size(); ++i) {
      sol.update(this->capacity_utilization(i, Q));
    }

    return sol.get_std_dev();
  }

  int n_moves(Solution &guiding, std::vector<int> &dict) {
    int count = 0u;
    int ptr1 = (this->get_route_size(my::indexGiantTour) % 2 != 0) ? static_cast<int>(std::ceil(this->get_route_size(my::indexGiantTour) / 2.0)) : static_cast<int>(this->get_route_size(my::indexGiantTour) / 2.0) - 1;
    int ptr2 = (this->get_route_size(my::indexGiantTour) % 2 != 0) ? ptr1 : ptr1 + 1;

    while (ptr1 >= 0 && ptr2 < this->get_route_size(my::indexGiantTour)) {
      if (this->get_route_size(my::indexGiantTour) % 2 != 0 && ptr1 == ptr2) {
        if (this->node_at(my::indexGiantTour, ptr1) != guiding.node_at(my::indexGiantTour, ptr1)) {
          ++count;
          dict.push_back(this->node_at(my::indexGiantTour, ptr1));
        }

        --ptr1;
        ++ptr2;
      }

      if (this->node_at(my::indexGiantTour, ptr1) != guiding.node_at(my::indexGiantTour, ptr1)) {
        ++count;
        dict.push_back(this->node_at(my::indexGiantTour, ptr1));
      }

      if (this->node_at(my::indexGiantTour, ptr2) != guiding.node_at(my::indexGiantTour, ptr2)) {
        ++count;
        dict.push_back(this->node_at(my::indexGiantTour, ptr2));
      }

      --ptr1;
      ++ptr2;
    }

    return count;
  }

  int proximity(Solution &benchmarkTour, const int &nbNode) {
    int count = 0u;
    int ptr1 = 1;
    int ptr2 = nbNode - 1;
    while (ptr1 < nbNode && ptr2 > 0 && ptr1 < ptr2 && ptr1 != ptr2) {  // two pointer technique
      const int iPosNode1 = this->positionNode[ptr1];
      const int iTourNode1 = this->routesNode[ptr1];
      const int iNextNode1 = this->node_at(iTourNode1, iPosNode1 + 1);
      const int iPrevNode1 = this->node_at(iTourNode1, iPosNode1 - 1);

      const int jPosNode1 = benchmarkTour.positionNode[ptr1];
      const int jTourNode1 = benchmarkTour.routesNode[ptr1];
      const int jNextNode1 = benchmarkTour.node_at(jTourNode1, jPosNode1 + 1);
      const int jPrevNode1 = benchmarkTour.node_at(jTourNode1, jPosNode1 - 1);

      if (iPrevNode1 != jPrevNode1 || iNextNode1 != jNextNode1) {
        ++count;
      }

      ++ptr1;

      const int iPosNode2 = this->positionNode[ptr2];
      const int iTourNode2 = this->routesNode[ptr2];
      const int iNextNode2 = this->node_at(iTourNode2, iPosNode2 + 1);
      const int iPrevNode2 = this->node_at(iTourNode2, iPosNode2 - 1);

      const int jPosNode2 = benchmarkTour.positionNode[ptr2];
      const int jTourNode2 = benchmarkTour.routesNode[ptr2];
      const int jNextNode2 = benchmarkTour.node_at(jTourNode2, jPosNode2 + 1);
      const int jPrevNode2 = benchmarkTour.node_at(jTourNode2, jPosNode2 - 1);

      if (iPrevNode2 != jPrevNode2 || iNextNode2 != jNextNode2) {
        ++count;
      }

      --ptr2;
    }

    if ((nbNode - 1) % 2 != 0) {
      const int ptr = static_cast<int>(std::ceil(nbNode / 2.0));
      const int iPosNode = this->positionNode[ptr];
      const int iTourNode = this->routesNode[ptr];
      const int iNextNode = this->node_at(iTourNode, iPosNode + 1);
      const int iPrevNode = this->node_at(iTourNode, iPosNode - 1);

      const int jPosNode = benchmarkTour.positionNode[ptr];
      const int jTourNode = benchmarkTour.routesNode[ptr];
      const int jNextNode = benchmarkTour.node_at(jTourNode, jPosNode + 1);
      const int jPrevNode = benchmarkTour.node_at(jTourNode, jPosNode - 1);

      if (iNextNode != jNextNode && iNextNode != jPrevNode) {
        ++count;
      }

      if (iPrevNode == 0 && jPrevNode != 0 && jNextNode != 0) {
        ++count;
      }
    }

    return count;
  }

private:
  void copy(const Solution &source) {
    this->Category = source.Category;
    this->Soltours = source.Soltours;
    this->Costours = source.Costours;
    this->Charge = source.Charge;
    this->SolObjcost = source.SolObjcost;
    this->Operational = source.Operational;
    this->routesNode = source.routesNode;
    this->positionNode = source.positionNode;
    this->routeSector = source.routeSector;
    this->nbNodes = source.nbNodes;
  }
};

inline std::ostream &operator<<(std::ostream &os, Solution &A) {  // print ostream
  os << "\n Tours-------" << std::endl;
  my::display_matrix(A.Soltours);

  os << "Lenght of tours (s)" << std::endl;
  my::display_matrix(A.Costours);
  std::cout << std::endl;

  os << "Charge of tours " << std::endl;
  my::display_matrix(A.Charge);
  std::cout << std::endl;

  os << "Operational time " << std::endl;
  my::display_vector(A.Operational);
  std::cout << std::endl;

  std::cout.setf(std::ios::left);  // left frame
  os << "Total length=" << A.SolObjcost;
  os.setf(std::ios::right);  //  right frame
  os << std::setw(20) << " Number of tours = " << A.Soltours.size() << std::endl;

  return os;
}

#endif //SOLUTION_H