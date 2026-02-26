/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef INSTANCE_H
#define INSTANCE_H

#include "LRUCache.h"

class Instance {
public:
  Instance() = default;

  ~Instance() = default;

  static int skip_line(FILE *fp) {
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
    return c;
  }

  void read_data(const char* fileIn) {
    this->sumQ = 0;
    this->C = 1;
    this->Q = std::vector<double>(this->C);

    FILE* file;

    int num_nodes_file;

    double vehicle_capa_file;
    double coord_x_file;
    double coord_y_file;
    double demand_node_file;
    double best_solution;

    int id_node;

    std::string instanceName = fileIn;
    std::regex target(".vrp");
    this->name = std::regex_replace(instanceName, target, "");

    file = fopen(fileIn, "r");

    if(file == NULL){
      printf("ERROR! Instance file not found \n");
      exit(1);
    }

    this->skip_line(file);
    this->skip_line(file);
    this->skip_line(file);

    fscanf(file, "DIMENSION : %d\n", &num_nodes_file);
    this->skip_line(file);
    this->skip_line(file);

    fscanf(file, "CAPACITY : %lf\n", &vehicle_capa_file);
    // fscanf(file, "BKS : %lf\n", &best_solution); // TODO debugging
    fscanf(file, "BKS : %lf\n", &best_solution);
    this->skip_line(file);

    this->N = num_nodes_file;
    this->Q[CVRP] = vehicle_capa_file;
    // this->BKS = best_solution; // TODO debugging

    this->demand.resize(this->N);
    for (auto i = 0; i < this->N; ++i){
      this->demand[i] = 0;
    }

    this->CoordinateNodes.resize(this->N, 2);
    this->customerList.reserve(num_nodes_file);

    for(auto i = 0; i < num_nodes_file; ++i){

      if (i > 0) {
        this->customerList.emplace_back(i);
      }

      fscanf(file, "%d %lf %lf\n", &id_node, &coord_x_file, &coord_y_file);
      this->CoordinateNodes.at(i, 0, coord_x_file);
      this->CoordinateNodes.at(i, 1, coord_y_file);
    }

    this->skip_line(file);

    for (auto i = 0; i < num_nodes_file; ++i){
      fscanf(file, "%d %lf\n", &id_node, &demand_node_file);
      this->demand[i] = demand_node_file;
      this->sumQ += this->demand[i];
    }

    std::vector<double> weightCost = {1.0};
    int precision = 0;

    this->distances.resize(this->N, this->N);

    for (int i = 0; i < this->N; ++i) {

      for (int j = i + 1; j < this->N; ++j) {
        const double tmp = my::round_val(std::sqrt(std::pow(this->CoordinateNodes.at(i, 0) - this->CoordinateNodes.at(j, 0), 2) + std::pow(this->CoordinateNodes.at(i, 1) - this->CoordinateNodes.at(j, 1), 2)) * weightCost[0],precision);
        this->distances.at(i, j, tmp);
        this->distances.at(j, i, tmp);
      }

    }

    this->polarAngle.resize(this->N);
    for (int i = 0; i < this->N; ++i) {
      this->polarAngle[i] = CircleSector::positive_mod(
                              static_cast<int>(std::round(32768.0 * atan2(this->CoordinateNodes.at(i, 1) - this->CoordinateNodes.at(my::depot, 1), this->CoordinateNodes.at(i, 0) - this->CoordinateNodes.at(my::depot, 0)) / 3.14159265359))
                            );
    }

    for (int c = 0; c < this->C; ++c) {
      if (this->Kmax[0] > this->N) {
        this->Kmax[0] = this->N;
      }

    }

  }  // save

  double distance(const int &category, const int &a, const int &b) { return this->distances.at(a, b); }

  void update_ksum(const int &newKsum) { this->Ksum = newKsum; }

  // void verification_all(const std::string &functionName) {
  //   this->verification_constraint_relaxed(functionName);
  //   std::vector<int> Ktrue;
  //   Ktrue.assign(this->C, 0);
  //   for (int c = 0; c < this->Ksum; ++c) {
  //     ++Ktrue[this->currSol.Category[c]];
  //   }
  //
  //   for (int c = 0; c < this->C; ++c) {
  //     if (Ktrue[c] > this->Kmax[c]) {
  //       std::cout << "VERIFICATION PROBLEMS: " << functionName << " number of vehicles exceed to Ksum " << Ktrue[c] << " >  " << this->Kmax[c] << std::endl;
  //       exit(-1);
  //     }
  //   }
  //
  //   Ktrue.clear();
  //
  //   for (int k = 0; k < this->Ksum; ++k) {
  //     int category = this->currSol.Category[k];
  //     if (this->currSol.get_route_load(k) > this->Q[category] + EPSILON) {
  //       std::cout << this->currSol << std::endl;
  //       std::cout << this->currSol.get_route_load(k) << " > " << this->Q[category] << std::endl;
  //       my::error_function(functionName, " Capacity violated by tour " + my::itos(k) + " Charge " + my::itos(this->currSol.get_route_load(k)) + " > " + my::itos(this->Q[category]));
  //     }
  //
  //     if (this->currSol.get_route_time(k) > this->L[category] + EPSILON) {
  //       std::cout << " currSol " << this->currSol << std::endl;
  //       my::error_function(functionName, " Max Length violated by tour  " + my::itos(k) + " Length Cost " +  my::itos(this->currSol.get_route_time(k)) + " << " + my::itos(this->L[category]));
  //     }
  //   }
  // }

  double length_path(const int &tour, const int &category, const int &posBegin, const int &posEnd) {
    double length = 0.0f;
    if (posBegin < posEnd) {
      if (category == this->currSol.Category[tour]) {
        (posBegin != 0) ? length = this->currSol.Costours[tour][posEnd] - this->currSol.Costours[tour][posBegin] : length = this->currSol.Costours[tour][posEnd];
      } else {
        for (int i = posBegin; i < posEnd; ++i) {
          length += this->distance(category, this->currSol.Soltours[tour][i], this->currSol.Soltours[tour][i + 1]);
        }

      }

    } else if (posBegin > posEnd) {
      my::error_function("length_path", " positions " + my::itos(posBegin) + " > " + my::itos(posEnd));
    }

    return length;
  }

  double length_reverse_path(const int &tour, const int &category, const int &posEnd, const int &posBegin) {
    std::string functionName = "length_reverse_path";
    double cost = 0.0f;
    if (posEnd > posBegin) {
      for (int i = posEnd; i > posBegin; --i) {
        cost += this->distance(category, this->currSol.Soltours[tour][i], this->currSol.Soltours[tour][i - 1]);
      }

    } else if (posEnd < posBegin) {
      my::error_function(functionName, " positions " + my::itos(posEnd) + " < " + my::itos(posBegin));
    }

    return cost;
  }

  private:
  void compute_cost_tours_charge() {  // calculate the sum of the length of each tour
    this->Ksum = this->currSol.get_size();
    double TotalLength = 0.0l;
    for (int k = 0; k < this->Ksum; ++k) {
      TotalLength += this->compute_cost_tours_charge_tour(k) + this->F[this->currSol.Category[k]];
    }

    this->currSol.SolObjcost = TotalLength;
  }

  void verification_constraint_relaxed(const std::string &functionName) {         // TODO debugging
    auto once = std::vector<int>(this->N, 0);
    for (auto& Soltour : this->currSol.Soltours) {
      for (int& i : Soltour) {
        once[i]++;
      }

    }

    for (int i = 1; i < this->N; ++i) {
      if (once[i] != 1) {
        (once[i] == 0) ? std::cout << "VERIFICATION PROBLEMS: " << functionName << " Node " << i << " is NOT visited " << std::endl : std::cout << "VERIFICATION PROBLEMS: " << functionName << " Node " << i << " is visited MORE than ONCE: " << once[i] << " times " << std::endl;
        exit(-1);
      }

    }

    Solution Saux(this->currSol);
    this->compute_cost_tours_charge();
    if (std::abs(this->currSol.SolObjcost - Saux.SolObjcost) > EPSILON) {
      std::cout << "ABS: " << std::abs(this->currSol.SolObjcost - Saux.SolObjcost) << std::endl;
      std::cout << "VERIFICATION PROBLEMS: " << functionName << " Solution Objective " << this->currSol.SolObjcost << " !=  " << Saux.SolObjcost << std::endl;
      std::cout << " Correct " << this->currSol << std::endl;
      exit(-1);
    }

    auto Nbtour = std::vector<int>(this->C);
    for (int k = 0; k < this->Ksum; ++k) {
      ++Nbtour[this->currSol.Category[k]];
    }

    for (int c = 0; c < this->C; ++c) {
      if (this->K[c] != Nbtour[c]) {
        my::error_function(functionName, " Problem categories: " + my::itos(c) + " size: " + my::itos(this->K[c]) + " != NbTour: " + my::itos(Nbtour[c]));
      }

    }

    for (int k = 0; k < this->Ksum; ++k) {
      if (this->currSol.Costours[k] != Saux.Costours[k]) {
        for (int j = 0; j < Saux.get_route_size(k); ++j) {
          if (std::abs(currSol.Costours[k][j] - Saux.Costours[k][j]) > EPSILON) {
            my::error_function(functionName, "Problems Costours: tour " + my::itos(k) + " node: " + my::itos(j) + " " + my::itos(currSol.Costours[k][j]) + "!=" + my::itos(Saux.Costours[k][j]));
          }

        }

      }

      if (this->currSol.Charge[k] != Saux.Charge[k]) {
        for (int j = 0; j < Saux.get_route_load(k); ++j) {
          if (std::abs(this->currSol.Charge[k][j] - Saux.Charge[k][j]) > EPSILON) {
            my::error_function(functionName, "Problems in Charge: tour " + my::itos(k) + " node " + my::itos(j) + " " + my::itos(this->currSol.Charge[k][j]) + " != " + my::itos(Saux.Costours[k][j]));
          }

        }

      }

    }

    if (this->currSol.get_operational_size() != this->Ksum) {
      std::cout << this->currSol << std::endl;
      my::error_function(functionName, " size of Operational time " + my::itos(this->currSol.get_operational_size()) + " Ksum " + my::itos(this->Ksum));
    }

    double operationalTime;
    for (int k = 0; k < this->Ksum; ++k) {
      operationalTime = this->D[this->currSol.Category[k]] * (this->currSol.get_route_size(k) - 2) + this->currSol.Costours[k][0];
      if (std::abs(operationalTime - this->currSol.get_route_time(k)) > EPSILON) {
        my::error_function(functionName,"Tour  " + my::itos(k) + " computed Operational time: " + my::itos(operationalTime) + " current: " + my::itos(this->currSol.get_route_time(k)));
      }

      if (this->currSol.get_route_time(k) > this->L[this->currSol.Category[k]]) {
        my::error_function(functionName, "Operational time: tour  " + my::itos(k) + " : " + my::itos(this->currSol.get_route_time(k)) + " > " + my::itos(this->L[this->currSol.Category[k]]));
      }

    }

    Saux = this->currSol;
    this->solution_from_vertices(functionName);
  }

  void solution_from_vertices(const std::string &functionName) {
    auto aux = std::vector<std::vector<int>>();
    if (this->currSol.routesNode[0] != this->N) {
      my::error_function(functionName, "Vertices 0 != N");
    }

    for (int k = 0; k < this->Ksum; ++k) {
      for (int i = 1; i < this->N; ++i) {
        if (this->currSol.routesNode[i] == k) {
          aux.push_back({i, this->currSol.positionNode[i]});
        }

      }

      std::sort(aux.begin(), aux.end(), [this](std::vector<int> i, std::vector<int> j) {
        if (i[1] == j[1]) {
          return (i[0] < j[0]);
        } else {
          return (i[1] < j[1]);
        }

      });

      const int aux_size = static_cast<int>(aux.size());
      for (int i = 0; i < aux_size; ++i) {
        if (static_cast<int>(aux[i][0]) != this->currSol.node_at(k, i + 1)) {
          my::error_function("solution_from_vertices in " + functionName,"tour " + my::itos(k) + " vertices " + my::itos(aux[i][0]) + "!=" + my::itos(this->currSol.node_at(k, i + 1)));
        }

        if (static_cast<int>(aux[i][0]) != 0 && this->currSol.positionNode[aux[i][0]] != i + 1) {
          my::error_function("Position Vertices solution_from_vertices in " + functionName, " tour " + my::itos(k) + " position vertex " + my::itos(aux[i][0]) + " " + my::itos(this->currSol.positionNode[aux[i][0]]) + "!=" + my::itos(i + 1));
        }

      }

      aux.clear();
    }

  }

  double compute_cost_tours_charge_tour(const int &tour) {
    int category = this->currSol.Category[tour];
    this->currSol.Costours[tour] = {0.0f};
    this->currSol.Charge[tour] = {0.0f};

    double charge = 0.0f;
    double auxLenght = 0.0f;
    for (auto itline = this->currSol.Soltours[tour].begin() + 1; itline != this->currSol.Soltours[tour].end(); ++itline) {
      auxLenght += this->distance(category, *(itline - 1), *itline);
      charge += this->demand[*itline];
      this->currSol.Costours[tour].push_back(auxLenght);
      this->currSol.Charge[tour].push_back(charge);
    }

    this->currSol.Costours[tour][0] = auxLenght;
    this->currSol.Charge[tour][0] = charge;
    return auxLenght;
  }

protected:
  int N;
  int C; // number of vehicle categories: For the CVRP, only one vehicle category is considered.
  int Ksum;

  std::vector<double> F;
  std::vector<double> Q;
  std::vector<double> L;
  std::vector<double> D;
  std::vector<double> demand;

  std::vector<int> Kmax;
  std::vector<int> K;
  std::vector<int> permNeighborhood;
  std::vector<int> polarAngle; // to measure cross route (inspired from HGS Vidal et al.)

  Matrix2DFlat<double> CoordinateNodes;
  Matrix2DFlat<double> distances;

public:
  clock_t globalTime;
  int vehEstimated;
  double stopsPerRoute;   // Kool, W., et al. (2022, April). In 12th DIMACS Implementation Challenge Workshop
  // double BKS; // TODO debugging
  double runTimeMax;
  double sumQ;
  unsigned int nbRun;
  std::vector<std::vector<double>> timeBestSols;

  bool epona;
  bool generated;
  std::vector<int> ellipseRed;
  std::vector<int> ellipseOrange;
  std::vector<int> nodeArea;

  std::vector<std::vector<std::vector<std::vector<int>>>> pathX;
  std::vector<std::vector<std::vector<std::vector<int>>>> pathY;

  Solution currSol;
  Solution localBest;
  Solution globalBest;
  Solution initialSol;
  Solution guideSol;

#ifdef GUIDANCE
  Solution reffSol;
#endif

  Matrix2DFlat<int> neighborsDistance;
  LRUCache customersCache;

  std::vector<int> customerList;

  std::string name;
};

#endif //INSTANCE_H