/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>

#include "base.h"
#include "Solution.h"
#include "GranularManagement.h"
#include "BitMatrix.h"
#include "PairBitMatrix.h"
#include "Instance.h"
#include "EvaluationMemory.h"
#include "EliteSol.h"
#include "SplitLabel.h"

class Solver : public Instance {
public:

  struct Saving {
    int node;
    int neighbor;
    double value;
  };

  Solver(Parameters &parameters) :  params(parameters) {
    this->generated = false;
    this->Ksum = 0;
    this->K = std::vector<int>(1);  // number of vehicles used for each category
    this->Kmax = std::vector<int>(1, INT_MAX);  // number of vehicles available for each category
    this->L = std::vector<double>(1, INT_MAX);  // maximum trip length available for each vehicle
    this->F = std::vector<double>(1);  // vehicles fixed costs for each category (unused yet)
    this->D = std::vector<double>(1);  // vehicles delivery time for each category (unused yet)
    // this->BKS = my::dummy;  // default value of the best known solution

    this->read_data(this->params.get_file_input());

    if (this->Q.size() > 1) { // TODO design for future research
      this->vehEstimated = my::dummy;
      this->stopsPerRoute = my::dummy;
    } else { // for CVRP
      this->vehEstimated = this->greedy_first_fit();
      this->stopsPerRoute = this->Q[CVRP] / (this->sumQ / (this->N - 1));
    }

    this->regretSelect = std::uniform_real_distribution<float>(0, 1);
    this->booleanDist = std::uniform_int_distribution<int>(0, 1);
    this->repairOrder = std::uniform_int_distribution<int>(0, 3);
    this->selectingNode = std::uniform_int_distribution<int>(Solver::get_cust_begin(), this->get_cust_end() - 1);

    for (int i = 0; i < NB_NEIGHBORHOOD; ++i) {
      this->permNeighborhood.push_back(i);
    }

    this->customersCache = LRUCache(this->params.get_cache_capacity(), this->N);
    this->runTimeMax = this->params.get_running_time_max();
    this->randEngine = SplitMix64Engine(this->params.get_seed());
    // this->primalIntegral = 0;        // TODO debugging
    // this->lastTimeSpend = 0;
    // this->lastCostFound = my::dummy;
    this->maxEliteSolSize = this->params.get_elite_size(this->stopsPerRoute);
    this->population = EliteSol(this->maxEliteSolSize, this->N);
    this->globalTime = clock();
  }

  ~Solver() = default;

  int greedy_first_fit() {
    std::sort(this->customerList.begin(), this->customerList.end(), [this](int i, int j) {
      return this->demand[i] > this->demand[j]; // sort based on high demand
    });

    struct tmp {
      double load{};
      std::vector<int> cust{};
    };

    auto tmps = std::vector<tmp>(this->N);
    auto res = 0;

    for (auto i : this->customerList) {

      int tmpSize = static_cast<int>(tmps.size());

      for (auto p = 0; p < tmpSize; p++) {

        auto& tmp = tmps[p];

        if (tmp.load + this->demand[i] <= this->Q[CVRP]) {

          tmp.load += this->demand[i];
          tmp.cust.push_back(i);

          if (p + 1 > res) {
            res = p + 1;
          }

          break;
        }

      }

    }

    return res;
  }

  // void finalize_primal_integral()  {       // TODO debugging
  //   this->primalIntegral += ((this->lastCostFound * (this->runTimeMax  - this->lastTimeSpend)) / (this->runTimeMax * this->BKS));
  //   this->primalIntegral -= 1;
  //   this->primalIntegral *= 100;
  // }

  void initialize_vertices()  {
    KDTree kd_tree(this->CoordinateNodes, this->N);

    const int neighborhoodTableSize = this->neighborsDistance.get_cols();
    int NbMaxNeighbors = this->neighborsDistance.get_cols();

    this->currSol.nbNodes = this->N;
    this->currSol.routesNode.resize(this->N, my::dummy);
    this->currSol.positionNode.resize(this->N, my::dummy);

    std::vector<std::vector<double>> distToDepot = std::vector<std::vector<double>>(this->get_cust_end(), std::vector<double>(2));  // distance from node to the depot    // TODO matrix2Dflat std::cout << "-> NbMaxNeighbors: " << NbMaxNeighbors << std::endl;  // debugging
    distToDepot[0] = {0, 0};

    for (int node = Solver::get_cust_begin(); node < this->get_cust_end(); ++node) {
      distToDepot[node] = {static_cast<double>(node), this->distance(CVRP, my::depot, node)};
    }

    std::sort(distToDepot.begin() + 1, distToDepot.begin() + this->N, [this](std::vector<double> i, std::vector<double> j) {
      if (i[1] == j[1]) {
        return (i[0] < j[0]);
      } else {
        return (i[1] < j[1]);
      }

    });

    for (int node = Solver::get_cust_begin(); node < this->get_cust_end(); ++node) {
      this->neighborsDistance.at(0, node, static_cast<int>(distToDepot[node][0]));
    }

    distToDepot.clear();
    int c = 0;

    for (int node = Solver::get_cust_begin(); node < this->get_cust_end(); ++node) {

      kd_tree.FIllNeighborhoodTable(
        this->CoordinateNodes.at(node, 0),
        this->CoordinateNodes.at(node, 1),
        neighborhoodTableSize, node,
        this->neighborsDistance
      );

      this->currSol.routesNode[node] = my::dummy;  // declare an un-routed node
      this->currSol.positionNode[node] = my::dummy;  // declare an un-routed node

    }

    this->currSol.routesNode[0] = this->get_cust_end();
    this->currSol.positionNode[0] = 0;
  }

  void update_vertices()  {
    int k = 0u;
    this->currSol.routeSector.reserve(this->currSol.Soltours.size());
    this->K = std::vector<int>(this->C);

    for (auto &Soltour : this->currSol.Soltours) {

      int i = 1;

      for (auto itt = Soltour.begin() + 1; itt != Soltour.end() - 1; ++itt) {

        this->currSol.routesNode[*itt] = k;
        this->currSol.positionNode[*itt] = i;

        if (*itt != my::depot) {
          if (i == 1) {
            this->currSol.routeSector[k].initialize(this->polarAngle[*itt]);
          } else {
            this->currSol.routeSector[k].extend(this->polarAngle[*itt]);
          }

        }

        ++i;
      }

      ++this->K[this->currSol.Category[k]];
      ++k;
    }

    this->Ksum = k;
  }

  void savings()  {
    this->savingTable.resize(this->C);

    int savingSize = this->params.get_saving_matrix_size(this->stopsPerRoute);
    int maxAdded = std::min(this->N - 1, savingSize);

    for (int c = 0; c < this->C; ++c) {

      auto savingTableCategorize = std::vector<Saving>();

      for (int node = Solver::get_cust_begin(); node < this->get_cust_end(); ++node) {

        for (int j = 0, added = 0; added < maxAdded && j < this->neighborsDistance.get_cols(); ++j) {

          int neighbor = this->neighborsDistance.at(node, j);

          if (node == neighbor || neighbor == my::depot) { continue; }

          const double value = this->distance(c, node, my::depot) + this->distance(c, my::depot, neighbor) - this->distance(c, node, neighbor);
          savingTableCategorize.push_back({node, neighbor, value});
          ++added;
        }

      }

      std::sort(savingTableCategorize.begin(), savingTableCategorize.end(),
                [](const Saving &i, const Saving &j) {
                  return i.value > j.value;
                });

      this->savingTable[c] = savingTableCategorize;
    }
  }

  void ts_2opt_intra()  {
    double oldCost = this->currSol.SolObjcost;
    int node = my::depot;
    int neighbor = my::depot;
    double minResult = 0.0f;
    while (node != my::dummy) {

      this->search_2opt(node, neighbor, minResult);

      if (node != my::dummy && neighbor != my::dummy) {

        int tour = this->currSol.routesNode[node];
        int possA = this->currSol.positionNode[node];
        int possB = this->currSol.positionNode[neighbor];

        this->taboo.put(node, neighbor);
        this->execute_2opt(tour, possA, possB, minResult);
        this->customersCache.put(node);
        this->customersCache.put(neighbor);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_realloc()  {
    double oldCost = this->currSol.SolObjcost;
    double result = 0.0f;
    int node = my::depot;
    int neighbor = my::depot;
    int tourNeighbor = my::dummy;

    while (node != my::dummy) {

      this->search_realloc(node, neighbor, tourNeighbor, result);

      if (node != my::dummy) {

        this->taboo.put(node, neighbor);
        this->execute_reallocation(node, neighbor, tourNeighbor);
        this->customersCache.put(node);
        this->customersCache.put(neighbor);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_swap()  {
    double oldCost = this->currSol.SolObjcost;
    double result = 0.0f;
    int node = my::depot;
    int neighbor = my::depot;
    while (node != my::dummy) {

      this->search_swap(node, neighbor, result);

      if (node != my::dummy) {

        this->taboo.put(node, neighbor);
        this->execute_swap(node, neighbor);
        this->customersCache.put(node);
        this->customersCache.put(neighbor);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_cross_exchange()  {
    EvaluationMemory moveTmp;
    std::vector<EvaluationMemory> bestMoves = {moveTmp};

    double result = 0.0f;
    double oldCost = this->currSol.SolObjcost;
    bestMoves[0].node = my::depot;
    bestMoves[0].neighbor = my::depot;

    while (bestMoves[0].node != my::dummy) {

      this->search_cross_exchange(bestMoves, result);

      if (bestMoves[0].node != my::dummy) {

        this->execute_moves(bestMoves);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          bestMoves[0].node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_2optstar()  {
    EvaluationMemory moveTmp;
    std::vector<EvaluationMemory> bestMoves = {moveTmp};

    double result = 0.0f;
    double oldCost = this->currSol.SolObjcost;
    bestMoves[0].node = my::depot;
    bestMoves[0].neighbor = my::depot;

    while (bestMoves[0].node != my::dummy) {

      this->search_2optstar(bestMoves, result);

      if (bestMoves[0].node != my::dummy) {

        this->execute_moves(bestMoves);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          bestMoves[0].node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_path_move(const int &sizeMin, const int &sizePath)  {
    EvaluationMemory bestMove;

    double oldCost = this->currSol.SolObjcost;
    bestMove.node = my::depot;

    while (bestMove.node != my::dummy) {

      this->search_path_move(bestMove, sizeMin, sizePath);

      if (bestMove.node != my::dummy) {

        this->taboo.put(bestMove.node, bestMove.neighbor);
        this->execute_move(bestMove);
        this->customersCache.put(bestMove.node);
        this->customersCache.put(bestMove.neighbor);

        if (this->currSol.SolObjcost < oldCost - EPSILON) {
          oldCost = this->currSol.SolObjcost;
          this->update_best();
        } else {
          bestMove.node = my::dummy;
          break;
        }

      }

    }

  }

  void ts_ec_realloc()  {
    int bestNode = my::depot;
    double resultTmp = 0.0f;

    Solution evalSol = this->currSol;

    const double oldCost = this->currSol.SolObjcost;

    this->search_ec_realloc(resultTmp, bestNode, evalSol);

    if (bestNode != my::dummy) {

      this->currSol = evalSol;
      this->update_ksum(this->currSol.get_size());

      if (this->currSol.SolObjcost < oldCost - EPSILON) {
        this->update_best();
      }

    }

  }

  void evaluation_movement(EvaluationMemory &EvalMove) {

    if (this->currSol.positionNode[EvalMove.node] == 1 && EvalMove.neighbor == my::depot) {
      EvalMove.result = INT_MAX;
      EvalMove.way = my::dummy;
    } else {
      this->evaluation_make_move(EvalMove);
    }

    this->evaluation_charge_opera_time(EvalMove);

  }

  void evaluation_make_move(EvaluationMemory &EvalMove) {
    std::string name = "Solver::evaluation_make_move";

    if (EvalMove.longPathMove <= 0) {
      my::error_function(name, " EvalMove.longPathMove " + my::itos(EvalMove.longPathMove));
    }

    if (EvalMove.node < 1 || EvalMove.neighbor < 0) {
      my::error_function(name, " node: " + my::itos(EvalMove.node) + "  neighbor " + my::itos(EvalMove.neighbor));
    } else {

      int node_next, poss_insert;

      if (EvalMove.neighbor == my::depot) {  // a new tour must be created
        node_next = my::depot;
        poss_insert = 0;
      } else {  // no new routes
        const int tourInsert = currSol.routesNode[EvalMove.neighbor];
        poss_insert = currSol.positionNode[EvalMove.neighbor];
        node_next = currSol.Soltours[tourInsert][poss_insert + 1];
      }

      const int node_head = EvalMove.node;
      const int tourRemove = this->currSol.routesNode[node_head];
      const int categoryNode = this->currSol.Category[tourRemove];
      int categoryNeighbor;

      (EvalMove.neighbor == my::depot) ? categoryNeighbor = EvalMove.category : categoryNeighbor = this->currSol.Category[this->currSol.routesNode[EvalMove.neighbor]];

      const int poss_pathMove = this->currSol.positionNode[node_head];
      const int node_tail = this->currSol.Soltours[tourRemove][poss_pathMove + EvalMove.longPathMove - 1];
      const int node_a = this->currSol.Soltours[tourRemove][poss_pathMove - 1];
      const int node_b = this->currSol.Soltours[tourRemove][poss_pathMove + EvalMove.longPathMove];

      double AddInsert = - this->distance(categoryNeighbor, EvalMove.neighbor, node_next) + this->distance(categoryNeighbor, EvalMove.neighbor, node_head) + this->distance(categoryNeighbor, node_tail, node_next) + this->length_path(tourRemove, categoryNeighbor, poss_pathMove, poss_pathMove + EvalMove.longPathMove - 1);

      double addRemove = this->distance(categoryNode, node_a, node_b) - this->distance(categoryNode, node_a, EvalMove.node) - this->distance(categoryNode, node_tail, node_b) - (this->currSol.Costours[tourRemove][poss_pathMove + EvalMove.longPathMove - 1] - this->currSol.Costours[tourRemove][poss_pathMove]);

      double AddInsertPrevious = INT_MAX;
      double addInsertPreviousRetro = INT_MAX;
      double AddInsertRetro = INT_MAX;

      if (poss_insert == 1) {
        AddInsertPrevious = - this->distance(categoryNeighbor, 0, EvalMove.neighbor) + this->distance(categoryNeighbor, 0, node_head) + this->distance(categoryNeighbor, node_tail, EvalMove.neighbor) + this->length_path(tourRemove, categoryNeighbor, poss_pathMove, poss_pathMove + EvalMove.longPathMove - 1);  // route size
        addInsertPreviousRetro = - this->distance(categoryNeighbor, 0, EvalMove.neighbor) + this->distance(categoryNeighbor, 0, node_tail) + this->distance(categoryNeighbor, node_head, EvalMove.neighbor) + this->length_reverse_path(tourRemove, categoryNeighbor, poss_pathMove + EvalMove.longPathMove - 1, poss_pathMove);  // gain total
      }

      if (AddInsertPrevious < AddInsert - EPSILON) {
        EvalMove.way = 1;  // move before to neighbor with Retro=false
        EvalMove.addInsert = AddInsertPrevious;
        EvalMove.addRemove = addRemove;
        EvalMove.result = EvalMove.addInsert + EvalMove.addRemove;
      } else {
        EvalMove.way = 0;
        EvalMove.addInsert = AddInsert;
        EvalMove.addRemove = addRemove;
        EvalMove.result = EvalMove.addInsert + EvalMove.addRemove;
      }

      if (EvalMove.longPathMove > 1) {
        AddInsertRetro = - this->distance(categoryNeighbor, EvalMove.neighbor, node_next) + this->distance(categoryNeighbor, EvalMove.neighbor, node_tail) + this->distance(categoryNeighbor, node_head, node_next) + this->length_reverse_path(tourRemove, categoryNeighbor, poss_pathMove + EvalMove.longPathMove - 1, poss_pathMove);

        if (AddInsertRetro + addRemove < EvalMove.result - EPSILON) {
          EvalMove.way = 2;
          EvalMove.addInsert = AddInsertRetro;
          EvalMove.addRemove = addRemove;
          EvalMove.result = EvalMove.addInsert + EvalMove.addRemove;
        }

        if (addInsertPreviousRetro + addRemove < EvalMove.result - EPSILON) {
          EvalMove.way = 3;
          EvalMove.addInsert = addInsertPreviousRetro;
          EvalMove.addRemove = addRemove;
          EvalMove.result = EvalMove.addInsert + EvalMove.addRemove;
        }

      }

    }

  }

  void execute_moves(std::vector<EvaluationMemory> &evalMoves) {

    for (auto &evalMove : evalMoves) {

      this->taboo.put(evalMove.node, evalMove.neighbor);
      this->execute_move(evalMove);
      this->customersCache.put(evalMove.node);
      this->customersCache.put(evalMove.neighbor);

    }

  }

  void execute_move(EvaluationMemory &EvalMove) {
    std::string name = "execute_move";

    if (EvalMove.longPathMove < 1) {
      my::error_function(name, " EvalMove.longPathMove " + my::itos(EvalMove.longPathMove));
    }

    if (EvalMove.node < 1 || EvalMove.neighbor < 0) {
      my::error_function(name, " path to Move or Receptor is " + my::itos(EvalMove.node));
    } else {

      int tourInsert, poss_insert;

      if (EvalMove.neighbor == 0 && (EvalMove.way == 0 || EvalMove.way == 2)) {
        this->make_new_tour(EvalMove.category);
        tourInsert = this->Ksum - 1;
        poss_insert = 0;
      } else {
        tourInsert = this->currSol.routesNode[EvalMove.neighbor];
        poss_insert = this->currSol.positionNode[EvalMove.neighbor];
      }

      int node_head = EvalMove.node;
      int tourRemove = this->currSol.routesNode[EvalMove.node];
      int poss_pathMove = this->currSol.positionNode[node_head];

      this->currSol.SolObjcost += EvalMove.result;

      auto it = this->currSol.Soltours.begin() + tourInsert;
      double charge, length;
      int node;

      if (EvalMove.way == 0 || EvalMove.way == 2) {
        ++poss_insert;
      }

      std::vector<int> line(this->currSol.Soltours[tourRemove].begin() + poss_pathMove,this->currSol.Soltours[tourRemove].begin() + poss_pathMove + EvalMove.longPathMove);

      if (EvalMove.way > 1) {
        std::reverse(line.begin(), line.end());
      }

      (*it).insert((*it).begin() + poss_insert, line.begin(), line.end());

      if (poss_insert == 1 && (EvalMove.way == 1 || EvalMove.way == 3)) {
        charge = 0;
        length = 0;
      } else {
        charge = this->currSol.Charge[tourInsert][poss_insert - 1];
        length = this->currSol.Costours[tourInsert][poss_insert - 1];
      }

      int Auxposs_insert = poss_insert;

      for (int i = 0; i < EvalMove.longPathMove; ++i) {

        node = this->currSol.Soltours[tourInsert][Auxposs_insert];

        this->currSol.routesNode[node] = tourInsert;
        this->currSol.positionNode[node] = Auxposs_insert;

        charge += this->demand[node];
        length += this->distance(this->currSol.Category[tourInsert], this->currSol.Soltours[tourInsert][Auxposs_insert - 1], node);

        this->currSol.Charge[tourInsert].insert(this->currSol.Charge[tourInsert].begin() + Auxposs_insert, charge);
        this->currSol.Costours[tourInsert].insert(this->currSol.Costours[tourInsert].begin() + Auxposs_insert, length);

        ++Auxposs_insert;

      }

      if (poss_insert != 1) {
        charge -= this->currSol.Charge[tourInsert][poss_insert - 1];
      }

      for (int i = Auxposs_insert; i < static_cast<int>((*it).size()); ++i) {
        this->currSol.positionNode[this->currSol.Soltours[tourInsert][i]] = i;
        this->currSol.Charge[tourInsert][i] += charge;
        this->currSol.Costours[tourInsert][i] += EvalMove.addInsert;
      }

      this->currSol.Costours[tourInsert][0] += EvalMove.addInsert;
      this->currSol.Charge[tourInsert][0] += charge;
      this->currSol.Operational[tourInsert] += EvalMove.addInsert + this->D[this->currSol.Category[tourInsert]] * (EvalMove.longPathMove);

      if (tourInsert == tourRemove && poss_insert < poss_pathMove) {
        poss_pathMove += EvalMove.longPathMove;
      }

      it = this->currSol.Soltours.begin() + tourRemove;
      (*it).erase((*it).begin() + poss_pathMove, (*it).begin() + poss_pathMove + EvalMove.longPathMove);

      this->currSol.Charge[tourRemove].erase(this->currSol.Charge[tourRemove].begin() + poss_pathMove, this->currSol.Charge[tourRemove].begin() + poss_pathMove + EvalMove.longPathMove);
      this->currSol.Costours[tourRemove].erase(this->currSol.Costours[tourRemove].begin() + poss_pathMove, this->currSol.Costours[tourRemove].begin() + poss_pathMove + EvalMove.longPathMove);

      for (int i = poss_pathMove; i < static_cast<int>((*it).size()); ++i) {
        this->currSol.Charge[tourRemove][i] -= charge;
        this->currSol.Costours[tourRemove][i] += EvalMove.addRemove;
        this->currSol.positionNode[this->currSol.Soltours[tourRemove][i]] = i;
      }

      this->currSol.Costours[tourRemove][0] += EvalMove.addRemove;
      this->currSol.Charge[tourRemove][0] -= charge;
      this->currSol.Operational[tourRemove] += EvalMove.addRemove - this->D[this->currSol.Category[tourRemove]] * (EvalMove.longPathMove);

      this->currSol.positionNode[0] = 0;
      this->currSol.routesNode[0] = this->N;

      if (this->currSol.Soltours[tourRemove].size() == 2) {
        this->erase_empty_tour(tourRemove);
      }

    }

  }

  std::vector<EvaluationMemory> cross_over_to_eval_move(const int &node, const int &neighbor, const bool &split, const bool &swapNode, const bool &swapNextNode) {
    EvaluationMemory Eval;
    auto Moves = std::vector<EvaluationMemory>();

    const int tourNode = this->currSol.routesNode[node];
    const int tourNeighbor = this->currSol.routesNode[neighbor];
    const int categoryNode = this->currSol.Category[tourNode];
    const int categoryNeighbor = this->currSol.Category[tourNeighbor];
    const int posNode = this->currSol.positionNode[node];
    const int posNeighbor = this->currSol.positionNode[neighbor];
    const int nextNode = this->currSol.Soltours[tourNode][posNode + 1];
    const int previous = this->currSol.Soltours[tourNeighbor][posNeighbor - 1];
    const int tailNeighbor = this->currSol.Soltours[tourNeighbor][this->currSol.Soltours[tourNeighbor].size() - 2];
    const int headNeighbor = (previous != 0) ? this->currSol.Soltours[tourNeighbor][1] : my::depot;

    if (node != my::dummy) {

      if (nextNode != my::depot) {

        Eval.node = nextNode;
        Eval.neighbor = my::depot;
        Eval.way = 0;
        Eval.longPathMove = this->ask_size_path(nextNode, this->N);
        Eval.addRemove = - this->length_path(tourNode, categoryNode, posNode, this->currSol.get_route_size(tourNode) - 1) + this->distance(categoryNode, node, my::depot);
        Eval.addInsert = this->length_path(tourNode, categoryNeighbor, posNode + 1, this->currSol.get_route_size(tourNode) - 1) + this->distance(categoryNeighbor, 0, nextNode);
        Eval.result = Eval.addInsert + Eval.addRemove;
        Eval.category = categoryNeighbor;

        Moves.push_back(Eval);
      }

      if (previous != my::depot) {

        Eval.node = neighbor;
        Eval.neighbor = my::depot;
        Eval.way = 0;
        Eval.longPathMove = this->ask_size_path(neighbor, this->N);
        Eval.addRemove = - this->length_path(tourNeighbor, categoryNeighbor, posNeighbor - 1, this->currSol.get_route_size(tourNeighbor) - 1) + this->distance(categoryNeighbor, previous, my::depot);
        Eval.category = this->currSol.Category[tourNeighbor];
        Eval.addInsert = this->length_path(tourNeighbor, categoryNeighbor, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 1) + this->distance(categoryNeighbor, my::depot, neighbor);
        Eval.result = Eval.addInsert + Eval.addRemove;

        Moves.push_back(Eval);
      }

      if (!split) {

        Eval.node = neighbor;
        Eval.neighbor = node;
        Eval.longPathMove = this->ask_size_path(neighbor, N);
        Eval.addRemove = -(this->currSol.Costours[tourNeighbor][0] - this->currSol.Costours[tourNeighbor][posNeighbor] + this->distance(categoryNeighbor, 0, neighbor));

        if (!swapNode) {
          Eval.way = 0;
          Eval.addInsert = this->distance(categoryNode, node, neighbor) + this->length_path(tourNeighbor, categoryNode, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 2) - this->distance(categoryNode, node, 0) + this->distance(categoryNode, tailNeighbor, 0);
        } else {
          Eval.way = 2;
          Eval.addInsert = this->distance(categoryNode, node, tailNeighbor) + this->length_reverse_path(tourNeighbor, categoryNode,this->currSol.get_route_size(tourNeighbor) - 2,posNeighbor) - this->distance(categoryNode, node, 0) + this->distance(categoryNode, neighbor, 0);
        }

        Eval.result = Eval.addInsert + Eval.addRemove;

        Moves.push_back(Eval);

        if (nextNode != my::depot && previous != my::depot) {

          Eval.node = this->currSol.Soltours[tourNeighbor][1];
          Eval.neighbor = nextNode;
          Eval.longPathMove = posNeighbor - 1;
          Eval.addRemove = -(this->currSol.Costours[tourNeighbor][posNeighbor - 1] + distance(categoryNeighbor, previous, my::depot));

          if (!swapNextNode) {
            Eval.way = 1;
            Eval.addInsert = this->currSol.Costours[tourNeighbor][posNeighbor - 1] + this->distance(categoryNeighbor, previous, nextNode) - this->distance(categoryNeighbor, my::depot, nextNode);
          } else {
            Eval.way = 3;
            Eval.addInsert = this->length_reverse_path(tourNeighbor, categoryNeighbor,posNeighbor - 1, 1) + this->distance(categoryNeighbor, my::depot, previous) + this->distance(categoryNeighbor, headNeighbor, nextNode) - this->distance(categoryNeighbor, my::depot, nextNode);
          }

          Eval.result = Eval.addInsert + Eval.addRemove;

          Moves.push_back(Eval);
        }
      } else {
        if (previous != my::depot) {

          Eval.node = headNeighbor;
          Eval.neighbor = node;
          Eval.longPathMove = posNeighbor - 1;
          Eval.addRemove = -this->length_path(tourNeighbor, categoryNeighbor, 0, posNeighbor - 1) - this->distance(categoryNeighbor, previous, my::depot);

          if (!swapNode) {
            Eval.way = 2;
            Eval.addInsert = this->length_reverse_path(tourNeighbor, categoryNode, posNeighbor - 1, 0) + this->distance(categoryNode, node, previous) - this->distance(categoryNode, node, my::depot);
          } else {
            Eval.way = 0;
            Eval.addInsert = this->length_path(tourNeighbor, categoryNode, 1,posNeighbor - 1) + this->distance(categoryNode, node, headNeighbor) + this->distance(categoryNode, previous, my::depot) - this->distance(categoryNode, node, my::depot);
          }

          Eval.result = Eval.addInsert + Eval.addRemove;

          Moves.push_back(Eval);
        }

        if (nextNode != my::depot) {

          Eval.node = neighbor;
          Eval.neighbor = nextNode;
          Eval.longPathMove = this->ask_size_path(neighbor, this->get_cust_end());
          Eval.addRemove = -(this->currSol.Costours[tourNeighbor][0] - this->currSol.Costours[tourNeighbor][posNeighbor] + this->distance(categoryNeighbor, 0, neighbor));

          if (!swapNextNode) {
            Eval.way = 3;
            Eval.addInsert = this->length_reverse_path(tourNeighbor, categoryNeighbor, this->currSol.get_route_size(tourNeighbor) - 1, posNeighbor) + this->distance(categoryNeighbor, neighbor, nextNode) - this->distance(categoryNeighbor, my::depot, nextNode);
          } else {
            Eval.way = 1;
            Eval.addInsert = this->length_path(tourNeighbor, categoryNeighbor, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 2) + this->distance(categoryNeighbor, tailNeighbor, nextNode) - this->distance(categoryNeighbor, my::depot, nextNode) + this->distance(categoryNeighbor, my::depot, neighbor);
          }

          Eval.result = Eval.addInsert + Eval.addRemove;
          Moves.push_back(Eval);
        }
      }
    } else {
      Eval.result = INT_MAX;
      Moves.push_back(Eval);
      Moves.push_back(Eval);
    }

    return Moves;
  }

  std::vector<EvaluationMemory> cross_exchange_to_eval_move(const int &node, const int &neighbor, const bool &swapNode, const bool &swapNeighbor) {
    EvaluationMemory move;
    auto moves = std::vector<EvaluationMemory>();

    if (node != my::dummy) {
      const int tour = this->currSol.routesNode[node];
      const int tourNeighbor = this->currSol.routesNode[neighbor];
      const int categoryNode = this->currSol.Category[tour];
      const int categoryNeighbor = this->currSol.Category[tourNeighbor];
      const int possNode = this->currSol.positionNode[node];
      const int possNeighbor = this->currSol.positionNode[neighbor];
      const int preNode = this->currSol.Soltours[tour][possNode - 1];
      const int tailNode = this->currSol.Soltours[tour][possNode + 1];
      const int postNode = this->currSol.Soltours[tour][possNode + 2];
      const int preNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor - 1];
      const int tailNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor + 1];
      const int postNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor + 2];

      move.node = neighbor;
      move.neighbor = tailNode;
      move.longPathMove = 2;

      if (!swapNeighbor) {
        move.way = 0;
        move.addInsert = this->distance(categoryNode, tailNode, neighbor) + this->distance(categoryNode, neighbor, tailNeighbor) + this->distance(categoryNode, tailNeighbor, postNode) - this->distance(categoryNode, tailNode, postNode);
      } else {
        move.way = 2;
        move.addInsert = this->distance(categoryNode, tailNode, tailNeighbor) + this->distance(categoryNode, tailNeighbor, neighbor) + this->distance(categoryNode, neighbor, postNode) - this->distance(categoryNode, tailNode, postNode);
      }

      move.addRemove = - this->distance(categoryNeighbor, preNeighbor, neighbor) - this->distance(categoryNeighbor, neighbor, tailNeighbor) - this->distance(categoryNeighbor, tailNeighbor, postNeighbor) + this->distance(categoryNeighbor, preNeighbor, postNeighbor);

      move.result = move.addInsert + move.addRemove;
      moves.push_back(move);

      move.node = node;
      move.longPathMove = 2;

      if (preNeighbor != my::depot || postNeighbor == my::depot) {
        move.neighbor = preNeighbor;
        move.way = 0;
        move.category = this->currSol.Category[this->currSol.routesNode[neighbor]];
      } else {
        move.neighbor = postNeighbor;
        move.way = 1;
      }

      if (!swapNode) {
        move.addInsert = this->distance(categoryNeighbor, preNeighbor, node) + this->distance(categoryNeighbor, node, tailNode) + this->distance(categoryNeighbor, tailNode, postNeighbor) - this->distance(categoryNeighbor, preNeighbor, postNeighbor);
      } else {
        move.addInsert = this->distance(categoryNeighbor, preNeighbor, tailNode) + this->distance(categoryNeighbor, tailNode, node) + this->distance(categoryNeighbor, node, postNeighbor) - this->distance(categoryNeighbor, preNeighbor, postNeighbor);
        move.way += 2;
      }

      if (!swapNeighbor) {
        move.addRemove = - this->distance(categoryNode, preNode, node) - this->distance(categoryNode, node, tailNode) - this->distance(categoryNode, tailNode, neighbor) + this->distance(categoryNode, preNode, neighbor);
      } else {
        move.addRemove = - this->distance(categoryNode, preNode, node) - this->distance(categoryNode, node, tailNode) - this->distance(categoryNode, tailNode, tailNeighbor) + this->distance(categoryNode, preNode, tailNeighbor);
      }

      move.result = move.addInsert + move.addRemove;
      moves.push_back(move);
    } else {
      move.result = INT_MAX;
      moves.push_back(move);
      moves.push_back(move);
    }

    return moves;
  }

  void evaluation_charge_opera_time(EvaluationMemory &EvalMove) {
    if (EvalMove.result != INT_MAX) {

      if (EvalMove.neighbor == my::depot) {

        EvalMove.operaTime = EvalMove.addInsert + this->D[0] * (EvalMove.longPathMove);

        if (this->currSol.positionNode[EvalMove.node] > 1) {
          EvalMove.charge = this->currSol.Charge[this->currSol.routesNode[EvalMove.node]] [this->currSol.positionNode[EvalMove.node] + EvalMove.longPathMove - 1] - this->currSol.Charge[this->currSol.routesNode[EvalMove.node]][this->currSol.positionNode[EvalMove.node] - 1];
        } else {
          EvalMove.charge = this->currSol.Charge[this->currSol.routesNode[EvalMove.node]][this->currSol.positionNode[EvalMove.node] + EvalMove.longPathMove - 1];
        }

      } else if (this->currSol.positionNode[EvalMove.node] > 1) {
        EvalMove.charge = this->currSol.Charge[this->currSol.routesNode[EvalMove.neighbor]][0] + this->currSol.Charge[this->currSol.routesNode[EvalMove.node]][this->currSol.positionNode[EvalMove.node] + EvalMove.longPathMove - 1] - this->currSol.Charge[this->currSol.routesNode[EvalMove.node]][this->currSol.positionNode[EvalMove.node] - 1];
        EvalMove.operaTime = this->currSol.Operational[this->currSol.routesNode[EvalMove.neighbor]] + EvalMove.addInsert + this->D[0] * (EvalMove.longPathMove);
      } else {
        EvalMove.charge = this->currSol.Charge[this->currSol.routesNode[EvalMove.neighbor]][0] + this->currSol.Charge[this->currSol.routesNode[EvalMove.node]][this->currSol.positionNode[EvalMove.node] + EvalMove.longPathMove - 1];
        EvalMove.operaTime = this->currSol.Operational[this->currSol.routesNode[EvalMove.neighbor]] + EvalMove.addInsert + this->D[0] * (EvalMove.longPathMove);
      }

    } else {
      EvalMove.charge = INT_MAX;
      EvalMove.operaTime = INT_MAX;
    }

  }

  void randomizeConcat() {
    std::vector<int> nodeHead = std::vector<int>(this->Ksum);
    for (int k = 0; k < this->Ksum; k++) {
      nodeHead[k] = this->currSol.Soltours[k][1];
    }

    std::shuffle(nodeHead.begin(), nodeHead.end(), this->randEngine);

    unsigned int oldKsum = this->Ksum;

    for (unsigned int k = 0; k < oldKsum; k++) {

      const int routenode = this->currSol.routesNode[nodeHead[k]];
      const int node = this->currSol.Soltours[routenode][this->currSol.Soltours[routenode].size() - 2];
      const int categoryNode = this->currSol.Category[routenode];

      for (int routeNeighbor = 0; routeNeighbor < this->Ksum; routeNeighbor++) {

        int neighbor = this->currSol.Soltours[routeNeighbor][1];
        routeNeighbor = this->currSol.routesNode[neighbor];
        const int categoryNeighbor = this->currSol.Category[routeNeighbor];

        if (routenode != routeNeighbor) {

          double moveDistance = this->distance(categoryNode, node, my::depot) + this->distance(categoryNeighbor, my::depot, neighbor) - this->distance(categoryNode, node, neighbor) + this->length_path(routeNeighbor, categoryNeighbor, 1,this->currSol.Soltours[routeNeighbor].size() - 1) - this->length_path(routeNeighbor, categoryNode, 1,this->currSol.Soltours[routeNeighbor].size() - 1);

          if (this->currSol.Soltours[routeNeighbor].size() > 3) {
            neighbor = this->currSol.Soltours[routeNeighbor][this->currSol.Soltours[routeNeighbor].size() - 2];
            moveDistance = this->distance(categoryNode, node, my::depot) + this->distance(categoryNeighbor, neighbor, my::depot) - this->distance(categoryNode, node, neighbor) - this->length_reverse_path(routeNeighbor, categoryNode, this->currSol.Soltours[routeNeighbor].size() - 2, 0) + this->length_path(routeNeighbor, categoryNeighbor, 0, this->currSol.Soltours[routeNeighbor].size() - 2);
          }

          EvaluationMemory evalMove;

          evalMove.node = this->currSol.Soltours[this->currSol.routesNode[neighbor]][1];
          evalMove.neighbor = node;
          evalMove.way = (this->currSol.positionNode[neighbor] == 1) ? 0 : 2;

          evalMove.longPathMove = this->ask_size_path(evalMove.node, this->N);
          evalMove.result = - moveDistance;
          evalMove.addInsert = evalMove.result + this->currSol.Costours[this->currSol.routesNode[neighbor]][0];
          evalMove.addRemove = - this->currSol.Costours[this->currSol.routesNode[neighbor]][0];

          this->execute_move(evalMove);

        }
      }

      if (this->currSol.Soltours.size() == 1) {
        break;  // early stop
      }
    }
  }

  bool is_dominated(SplitLabel &label1, SplitLabel &label2) {

    if (label1.cost < label2.cost - EPSILON) {
      return false;
    } else {

      for (int c = 0; c < this->C; c++) {

        if (label1.vhlUsed[c] < label2.vhlUsed[c] && this->Kmax[c] < this->get_cust_end() - 1) {
          return false;
        }

      }

    }

    return true;
  }

  bool is_dominated(SplitLabel &label, std::vector<SplitLabel> &labels) {
    unsigned int i = 0u;
    bool tmp = false;

    while (i < labels.size() && !tmp) {
      tmp = this->is_dominated(label, labels[i]);
      i++;
    }

    return tmp;
  }

  void add_split_label(std::vector<SplitLabel> &labels, SplitLabel &newLabel) {
    unsigned const int sizeMax = 10u;
    labels.emplace_back(newLabel);

    if (labels.size() > sizeMax) {
      std::sort(labels.begin(), labels.end(), [this](SplitLabel a, SplitLabel b) {
        return a.cost < b.cost;
      });

      labels.erase(labels.begin() + sizeMax, labels.end());
    }

  }

  void split_label_to_current_solution(SplitLabel *bestSplit, const int &tour, const int &oldCategory) {
    SplitLabel *labelTmp = bestSplit;
    SplitLabel *prevLabel = bestSplit->previousLabel;

    while (labelTmp->index > 0) {

      EvaluationMemory evalMove;

      evalMove.node = this->currSol.Soltours[tour][prevLabel->index + 1];
      evalMove.neighbor = my::depot;
      evalMove.way = 0;
      evalMove.longPathMove = this->ask_size_path(evalMove.node, this->N);
      evalMove.addInsert = labelTmp->lengthTrips;
      evalMove.addRemove = -this->length_path(tour, oldCategory, prevLabel->index,labelTmp->index + 1) + this->distance(oldCategory,this->currSol.Soltours[tour][prevLabel->index],my::depot);
      evalMove.result = evalMove.addInsert + evalMove.addRemove;
      evalMove.category = labelTmp->C;

      this->execute_move(evalMove);

      labelTmp = prevLabel;
      prevLabel = prevLabel->previousLabel;

    }

    labelTmp = NULL;
    delete labelTmp;
    prevLabel = NULL;
    delete prevLabel;
  }

  bool split_giant_tour(const int &tour) {
    const int oldCategory = this->currSol.Category[tour];
    const int size = this->currSol.Soltours[tour].size() - 2;
    auto labels = std::vector<std::vector<SplitLabel>>(size + 1);

    SplitLabel labelTmp = SplitLabel();

    labelTmp.index = 0;
    labelTmp.cost = 0;
    labelTmp.C = my::dummy;
    labelTmp.lengthTrips = 0;
    labelTmp.previousLabel = &labelTmp;

    for (int c = 0; c < this->C; c++) {
      labelTmp.vhlUsed.push_back(0);
    }

    labels[0].push_back(labelTmp);

    double dist;

    for (int i = 1; i <= size; i++) {

      for (unsigned int l = 0; l < labels[i - 1].size(); l++) {

        for (int c = 0; c < this->C; c++) {

          if (labels[i - 1][l].vhlUsed[c] < this->Kmax[c]) {

            int j = i;
            int node = this->currSol.Soltours[tour][j];
            int prev = this->currSol.Soltours[tour][j - 1];
            dist = labels[i - 1][l].cost;
            double load = this->demand[node];

            while (load <= this->Q[c] && j <= size ) {

              if (i == j) {
                dist += this->F[c] + this->distance(c, my::depot, node) + this->distance(c, node, my::depot);
              } else {
                dist += this->distance(c, prev, node) + this->distance(c, node, my::depot) - this->distance(c, prev, my::depot);
              }

              SplitLabel child = SplitLabel(labels[i - 1][l]);
              child.cost = dist;

              for (int c1 = 0; c1 < this->C; c1++) {

                if (c1 == c) {
                  child.vhlUsed[c1] = labels[i - 1][l].vhlUsed[c1] + 1;
                } else {
                  child.vhlUsed[c1] = labels[i - 1][l].vhlUsed[c1] + 0;
                }
              }

              if (!this->is_dominated(child, labels[j]) ){

                int s = 0;

                while (s < static_cast<int>(labels[j].size())) {
                  if (this->is_dominated(labels[j][s], child)) {
                    labels[j].erase(labels[j].begin() + s);
                  } else {
                    s++;
                  }

                }

                child.index = j;
                child.C = c;
                child.lengthTrips = dist - this->F[c] - labels[i - 1][l].cost;
                child.cost = dist;
                child.previousLabel = &labels[i - 1][l];
                this->add_split_label(labels[j], child);
              }

              j++;
              prev = node;
              node = this->currSol.Soltours[tour][j];
              load += this->demand[node];
            }

          }

        }

      }

    }

    int indBest = my::dummy;
    double minCost = std::numeric_limits<float>::max();
    const int labelSize = static_cast<int>(labels[size].size());

    for (int i = 0; i < labelSize; i++) {

      if (labels[size][i].cost < minCost - EPSILON) {
        indBest = i;
        minCost = labels[size][i].cost;
      }

    }

    if (indBest != my::dummy) {
      this->split_label_to_current_solution(&labels[size][indBest], tour, oldCategory);
      return true;
    } else {
      return false;
    }

  }

  void neighborhood_improvement(const double &index) {
    const int maxIter = this->population.get_size();
    const double regret = this->regretIndex * (REGRET_CONST - index);
    const int tourIndex = INDEX_PERTURBATION(index);

    auto perturbSols = std::vector<Solution>();

    for (int iter = 0; iter < maxIter; ++iter) {

      this->customersCache.reset();
      this->currSol = this->population.get_begin()->second;

      const double lastSolCost = this->currSol.SolObjcost;

      this->population.erase_sol_by_val(lastSolCost);
      this->update_ksum(this->currSol.get_size());

      this->destroy_tour(tourIndex, REPAIR_CONST, regret);

    __repeat__:
      for (int level = 1; level <= 2; ++level) {

        const bool improv = this->run_level(level);

        if (improv) {
          goto __repeat__;
        }

      }

      perturbSols.push_back(this->currSol);
      if (this->currSol.SolObjcost < lastSolCost - EPSILON) {
        break;
      }

    }

    for (auto& sol : perturbSols) {
      const bool isInserted = this->population.set_sols(sol);
    }

    const bool isInserted = this->population.set_sols(this->globalBest);

  }

  bool destroy_paired_tour() {
    this->unaffectedNodes.clear();

    const int first = this->selectingNode(this->randEngine);
    int second = my::dummy;

    for (int n = 0; n < this->maxGranularSize; ++n) {

      int eval = this->neighborsDistance.at(first, n);

      if (first == eval || eval == my::depot){ continue; }

      if (this->currSol.routesNode[first] != this->currSol.routesNode[eval]) {
        second = eval;
        break;
      }

    }

    // destroy
    this->erase_tour(this->currSol.routesNode[first]);

    if (second != my::dummy && second != my::depot) {
      this->erase_tour(this->currSol.routesNode[second]);
    }

    if (!this->remainingFreeNodes.empty()) {

      for (auto &remaining : this->remainingFreeNodes) {

        if (this->currSol.routesNode[remaining] != my::dummy) {
          this->remove_node(remaining);
          this->unaffectedNodes.push_back(remaining);
        }

      }

      this->remainingFreeNodes.clear();
    }

    // rebuild
    switch (this->booleanDist(this->randEngine)) {
      case 0:
        std::shuffle(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), this->randEngine);
        break;
      case 1:
        std::sort(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), [this](int i, int j) {
          return this->demand[i] > this->demand[j];
        });
        break;
    }

    for (int &node : this->unaffectedNodes) {
      double minCost = std::numeric_limits<float>::max();
      double curCost = my::dummy;
      int bestPos = my::dummy;
      int bestRoute = my::dummy;

      for (int k = 0; k < this->Ksum; ++k) {

        int category = this->currSol.Category[k];

        if (this->demand[node] + this->currSol.Charge[k][0] <= this->Q[category]) {

          for (int i = Solver::get_cust_begin(); i < this->currSol.get_route_size(k); ++i) {

            curCost = this->cost_insertion(node, k, i);

            if (curCost < minCost) {
              bestPos = i;
              bestRoute = k;
              minCost = curCost;
            }

          }

        }

      }

      for (int c = 0; c < this->C; ++c) {

        if (this->K[c] < this->Kmax[c] && this->demand[node] <= this->Q[c]) {

          curCost = this->distance(c, 0, node) + this->distance(c, node, 0) + this->F[c];

          if (curCost < minCost) {
            bestPos = 0;
            minCost = curCost;
            bestRoute = c;
          }

        }

      }

      if (bestPos > 0) {
        this->currSol.routesNode[node] = bestRoute;
        this->insert_node(node, this->currSol.node_at(bestRoute, bestPos - Solver::get_cust_begin()));
        this->customersCache.put(node);
      } else {
        if (this->booleanDist(this->randEngine)) {
          this->make_new_tour(node, 0); // TODO improved it for HVRP
          this->customersCache.put(node);
        } else {
          this->remainingFreeNodes.push_back(node);
        }

      }

    }

    for (int &tmp : this->unaffectedNodes) {

      if (this->currSol.routesNode[tmp] == my::dummy) {
        return false;
      }

    }

    this->unaffectedNodes.clear();
    return true;

  }

  bool destroy_paired_tour2() {
    this->unaffectedNodes.clear();
    this->update_vertices();

    int first = this->currSol.routesNode[this->selectingNode(this->randEngine)];
    int second = my::dummy;

    std::vector<int> route1(this->currSol.get_size());
    std::iota(route1.begin(), route1.end(), 0);

    std::vector<int> route2(this->currSol.get_size());
    std::iota(route2.begin(), route2.end(), 0);

    std::shuffle(route1.begin(), route1.end(), this->randEngine);
    std::shuffle(route2.begin(), route2.end(), this->randEngine);

    for (int i1 : route1) {

      for (int i2 : route2) {

        if (i1 != i2 && i1 < i2) {
          if (CircleSector::overlap(this->currSol.routeSector[i1], this->currSol.routeSector[i2])) {
            first = i1;
            second = i2;
            goto __destroy_paired2__;;
          }

        }

      }

    }

  __destroy_paired2__:

    this->erase_tour(first);

    if (second != my::dummy) {
      this->erase_tour(second);
    }

    if (!this->remainingFreeNodes.empty()) {

      for (auto &remaining : this->remainingFreeNodes) {

        if (this->currSol.routesNode[remaining] != my::dummy) {
          this->remove_node(remaining);
          this->unaffectedNodes.push_back(remaining);
        }

      }

      this->remainingFreeNodes.clear();
    }

    // rebuild
    switch (this->booleanDist(this->randEngine)) {
      case 0:
        std::shuffle(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), this->randEngine);
        break;
      case 1:
        std::sort(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), [this](int i, int j) {
          return this->demand[i] > this->demand[j];
        });
        break;
    }

    for (int &node : this->unaffectedNodes) {

      double minCost = std::numeric_limits<float>::max();
      double curCost = my::dummy;
      int bestPos = my::dummy;
      int bestRoute = my::dummy;

      for (int k = 0; k < this->Ksum; ++k) {

        int category = this->currSol.Category[k];

        if (this->demand[node] + this->currSol.Charge[k][0] <= this->Q[category]) {

          for (int i = Solver::get_cust_begin(); i < this->currSol.get_route_size(k); ++i) {

            curCost = this->cost_insertion(node, k, i);

            if (curCost < minCost) {
              bestPos = i;
              bestRoute = k;
              minCost = curCost;
            }

          }

        }

      }

      if (bestPos == my::dummy) {

        for (int c = 0; c < this->C; ++c) {

          if (this->K[c] < this->Kmax[c] && this->demand[node] <= this->Q[c]) {

            curCost = this->distance(c, 0, node) + this->distance(c, node, 0) + this->F[c];

            if (curCost < minCost) {
              bestPos = 0;
              minCost = curCost;
              bestRoute = c;
            }

          }

        }

      }

      if (bestPos > 0) {
        this->currSol.routesNode[node] = bestRoute;
        this->insert_node(node, this->currSol.node_at(bestRoute, bestPos - Solver::get_cust_begin()));
        this->customersCache.put(node);
      } else {
        if (this->booleanDist(this->randEngine)) {
          this->make_new_tour(node, 0); // TODO improved it for HVRP (randomize category)
          this->customersCache.put(node);
        } else {
          this->remainingFreeNodes.push_back(node);
        }
      }

    }

    for (int &tmp : this->unaffectedNodes) {
      if (this->currSol.routesNode[tmp] == my::dummy) {
        return false;
      }
    }

    this->unaffectedNodes.clear();
    return true;
  }

  void destroy_tour(const int &tourIndex, const int &repairIndex, const double &regretIndex){
    this->unaffectedNodes.clear();
    const Solution oldSol = this->currSol;

    int seed1 = this->selectingNode(this->randEngine);
    int seed2 = my::dummy;

    if (tourIndex > 1 && seed1 != my::dummy) {

      for (int n = 0; n < this->granularSizeNodes.get_neighbor(seed1); ++n) {
        int eval2 = this->neighborsDistance.at(seed1, n);

        if (eval2 == seed1 || eval2 == my::depot) { continue; }

        if (this->currSol.routesNode[seed1] != this->currSol.routesNode[eval2]) {
          seed2 = eval2;
          goto __find_seed2__;
        }

      }

    }

    __find_seed2__:

    int seed3 = my::dummy;

    if (tourIndex > 2 && seed2 != my::dummy) {

      for (int n = 0; n < this->granularSizeNodes.get_neighbor(seed2); ++n) {

        int eval3 = this->neighborsDistance.at(seed2, n);

        if (eval3 == seed1 || eval3 == seed2 || eval3 == my::depot) { continue; }

        if (this->currSol.routesNode[seed1] != this->currSol.routesNode[eval3] && this->currSol.routesNode[seed2] != this->currSol.routesNode[eval3]) {
          seed3 = eval3;
          goto __find_seed3__;
        }

      }

    }

  __find_seed3__:

    // destroy routes
    this->erase_tour(this->currSol.routesNode[seed1]);
    if (seed2 != my::dummy && seed2 != my::depot) {
      this->erase_tour(this->currSol.routesNode[seed2]);
    }

    if (seed3 != my::dummy && seed3 != my::depot) {
      this->erase_tour(this->currSol.routesNode[seed3]);
    }

    // rebuild
    const bool repaired = this->repair_all_unaffected_nodes(repairIndex, regretIndex);

    if (!repaired) {
      this->currSol = oldSol;
      this->update_ksum(this->currSol.get_size());
    }

    if (this->currSol.SolObjcost < this->globalBest.SolObjcost - EPSILON) {
      this->update_best();
    }

  }

  void destroy_tour2(const int &tourIndex, const int &repairIndex, const double &regretIndex) {
    this->unaffectedNodes.clear();
    const Solution oldSol = this->currSol;

    this->update_vertices();

    int i = this->currSol.routesNode[this->selectingNode(this->randEngine)];
    int j = my::dummy;

    if (tourIndex > 1) {

      std::vector<int> route1(this->currSol.get_size());

      std::iota(route1.begin(), route1.end(), 0);
      std::shuffle(route1.begin(), route1.end(), this->randEngine);

      std::vector<int> route2(this->currSol.get_size());

      std::iota(route2.begin(), route2.end(), 0);
      std::shuffle(route2.begin(), route2.end(), this->randEngine);

      for (int i1 : route1) {
        for (int i2 : route2) {
          if (i1 != i2 && i1 < i2) {
            if (CircleSector::overlap(this->currSol.routeSector[i1], this->currSol.routeSector[i2])) {
              i = i1;
              j = i2;
              goto __destroy2__;
            }
          }
        }
      }

    }

  __destroy2__:
    // destroy
    this->erase_tour(i);
    if (j != my::dummy) {
      this->erase_tour(j);
    }

    // rebuild
    const bool repaired = this->repair_all_unaffected_nodes(repairIndex, regretIndex);
    if (!repaired) {
      this->currSol = oldSol;
      this->update_ksum(this->currSol.get_size());
    }

    if (this->currSol.SolObjcost < this->globalBest.SolObjcost - EPSILON) {
      this->update_best();
    }

  }

  void clarke_and_wrigth_solver() {
#ifdef GUIDANCE
    this->reffSol.SolObjcost = 0.0;
#endif

    this->update_vertices();

    this->population.reset();
    this->Ksum = 0;

    auto isInBeginningTour = std::vector<bool>(this->N);
    auto isInEndingTour = std::vector<bool>(this->N);

    EvaluationMemory evalMove;

    this->erase_current_solution();

    this->currSol.SolObjcost = 0.0f;

    for (int c = 0; c < this->C; ++c) {

      this->build_one_customer_routes(isInBeginningTour, isInEndingTour, c);
      int idxSavings = 0;

      while (idxSavings < static_cast<int>(this->savingTable[c].size()) - 1) {

        const int ii = this->savingTable[c][idxSavings].node;
        const int jj = this->savingTable[c][idxSavings].neighbor;
        const int ir = this->currSol.routesNode[ii];
        const int jr = this->currSol.routesNode[jj];
        const double savingVal = this->savingTable[c][idxSavings].value;

        if ( ir != jr &&
            (isInEndingTour[ii] && isInBeginningTour[jj]) &&
            this->currSol.get_route_load(jr) + this->currSol.get_route_load(ir) <= this->Q[c] &&
            this->currSol.get_route_time(ir) + this->currSol.get_route_time(jr) - this->savingTable[c][idxSavings].value <= this->L[c]) {

          int routeMove =  my::dummy;
          if (isInEndingTour[ii] && isInBeginningTour[jj]) {
            evalMove.node = jj;
            evalMove.neighbor = ii;
            routeMove = jr;
          } else if (isInEndingTour[jj] && isInBeginningTour[ii]) {
            evalMove.node = ii;
            evalMove.neighbor = jj;
            routeMove = ii;
          } else {
            continue;
          }

          evalMove.way = 0;
          evalMove.longPathMove = this->ask_size_path(evalMove.node, this->get_cust_end());
          evalMove.result = -savingVal;
          evalMove.addInsert = evalMove.result + this->currSol.Costours[routeMove][0];
          evalMove.addRemove = -this->currSol.Costours[routeMove][0];

          this->execute_move(evalMove);

          if (isInEndingTour[ii] && isInBeginningTour[jj]) {
            isInEndingTour[ii] = false;
            isInBeginningTour[jj] = false;
          } else if (isInEndingTour[jj] && isInBeginningTour[ii]) {
            isInEndingTour[jj] = false;
            isInBeginningTour[ii] = false;
          }

        }

        ++idxSavings;
      }

      if (this->K[c] > this->Kmax[c]) {
        this->erase_exceeded_routes(c);
      }

    }

    bool availableSol = this->repair_all_unaffected_nodes(1, 0);

    if (availableSol) {
      this->update_vertices();
      this->update_best();

      unsigned int counter = 0;
      Solution lastSol = this->currSol;
      double lastSolCost = std::numeric_limits<float>::max();

      while (this->population.get_size() < 2) {  // adjustable

        this->currSol = lastSol;
        this->update_ksum(this->currSol.get_size());

        if (this->population.get_size() == 0 || this->currSol.get_size() <= this->vehEstimated) {

        __end__:

          bool is_inserted = this->population.set_sols(this->currSol);

          if (is_inserted) {

#ifdef GUIDANCE
            if (this->currSol.SolObjcost > this->reffSol.SolObjcost) {
              this->reffSol = this->currSol;
            }

#endif
            if (this->population.get_size() == 2) {
              break;
            }

          }

        }

        if (counter >= this->params.get_init_recombination()) {
          break;
        }

        const bool feasibleSol = this->destroy_paired_tour2();

        if (feasibleSol) {

          if (this->currSol.SolObjcost < lastSolCost - EPSILON || this->currSol.proximity(lastSol, this->N) >= 4) {

            lastSol = this->currSol;
            lastSolCost = lastSol.SolObjcost;
            this->remainingFreeNodes.clear();

            if (this->currSol.get_size() <= this->vehEstimated) {
              goto __end__;
            }

          }

        }

        ++counter;
      }
    }
  }

  bool stopping_criteria() {

    const clock_t currentTime = clock();

    const double timeSpend = static_cast<double>(currentTime - this->globalTime) / static_cast<double>(CLOCKS_PER_SEC);

    if ((this->nbRun != 0 && this->nbRun % 500 == 0) && this->runTimeMax - EPSILON >= timeSpend) {
#ifdef VERBOSE
      std::cout << "- current solution: " << this->globalBest.SolObjcost << "\t| time: " << timeSpend << "s\t| iteration: " << nbRun << std::endl;
#endif
      this->timeBestSols.push_back({timeSpend, this->globalBest.SolObjcost});
    }

    return (timeSpend < this->runTimeMax);
  }

  bool run_level(const int &level);

  bool neighborhood_operator(const int &index);

  void iteratively_evaluate_neighbour(Solution &initial, Solution &guiding, const int &maxIter, std::vector<int> &dict, const double &regret);

  void gt_path_relinking(const double &index);

  void generate_solutions();

  void run();

  void write_solution_file() {
    std::string file = "./results/" + this->name + ".sol";

    if (file != "none") {

      unsigned int numberTours = this->globalBest.Soltours.size();

      std::ofstream fa;
      fa.open(file);
      fa.precision(10);
      fa << this->globalBest.SolObjcost << std::endl;
      fa << numberTours << std::endl;

      // operational time (moving time + depository time) of each route
      fa << "operational time of each route -------" << std::endl;
      for (unsigned int i = 0; i < numberTours; ++i) { fa << this->globalBest.Operational[i] << " "; }

      fa << std::endl;

      // load description of each route
      fa << "load description of each route -------" << std::endl;
      for (unsigned int i = 0; i < numberTours; ++i) { fa << this->globalBest.Charge[i][0] << " "; }

      fa << std::endl;

      // customers delivered in each route
      fa << "customers delivered in each route -------" << std::endl;
      for (unsigned int i = 0; i < numberTours; ++i) {
        for (int j : this->globalBest.Soltours[i]) { fa << j << " "; }
        fa << std::endl;
      }

      // write detail of tour length
      fa << "detail of tour length -------" << std::endl;
      for (unsigned int i = 0; i < numberTours; ++i) {
        fa << 0 << " ";
        for (unsigned int j = 1; j < this->globalBest.Soltours[i].size(); ++j) { fa << this->globalBest.Costours[i][j] << " "; }
        fa << std::endl;
      }

      // write detail of tour load
      fa << "detail of tour load -------" << std::endl;
      for (unsigned int i = 0; i < numberTours; ++i) {
        fa << 0 << " ";
        for (unsigned int j = 1; j < this->globalBest.Soltours[i].size(); ++j) { fa << this->globalBest.Charge[i][j] << " "; }
        fa << std::endl;
      }

      fa.close();
    }

  }

private:

  void update_best() {

    if (this->currSol.SolObjcost < this->localBest.SolObjcost - EPSILON) {
      this->localBest = Solution(this->currSol);
    }

    if (this->currSol.SolObjcost < this->globalBest.SolObjcost - EPSILON) {
      this->globalBest = Solution(this->currSol);
      const clock_t currentTime = clock();
      const double timeSpend = static_cast<double>((currentTime - this->globalTime)) / static_cast<double>(CLOCKS_PER_SEC);
#ifdef LOG_OUTPUT
      std::cout << "* Cost = " << this->globalBest.SolObjcost << " ; Gap: " << (this->globalBest.SolObjcost - this->BKS) * 100 / this->BKS << " ; Time = " << timeSpend << std::endl;  // debugging // TODO log
#endif
      // if (this->globalBest.SolObjcost < this->lastCostFound || this->lastCostFound == -1) {        // TODO debugging
      //   auto t_iMinus1 = this->lastTimeSpend;
      //   if (this->lastCostFound == -1) {
      //     this->lastCostFound = this->globalBest.SolObjcost;
      //     return; // empty return
      //   }
      //   this->primalIntegral += (this->lastCostFound * (timeSpend - t_iMinus1) / (this->runTimeMax * this->BKS));
      //   this->lastCostFound = this->globalBest.SolObjcost;
      //   this->lastTimeSpend = timeSpend;
      // }

    }

  }

  void current_to_best() {

    if (this->localBest.SolObjcost < this->currSol.SolObjcost - EPSILON) {
      this->currSol = Solution(this->localBest);
      this->update_vertices();
    }

  }

  void erase_current_solution() {
    this->currSol = Solution();

    this->Ksum = 0u;
    this->K = std::vector<int>(this->C);

    this->currSol.routeSector.clear();

    this->currSol.nbNodes = this->N;
    this->currSol.routesNode.resize(this->get_cust_end());
    this->currSol.positionNode.resize(this->get_cust_end());

    for (int i = 1; i < this->get_cust_end(); ++i) {
      this->currSol.routesNode[i] = my::dummy;
      this->currSol.positionNode[i] = my::dummy;
    }

  }

  int search_unfeasible_tour() {

    for (int k = 0; k < this->Ksum; ++k) {

      if (!this->currSol.feasible_tour(k, this->Q[this->currSol.Category[k]], this->L[this->currSol.Category[k]])) {
        return k;
      }

    }

    return my::dummy;
  }

  void make_new_tour(const int &node, const int &category) {
    double distanceTour = this->distance(category, my::depot, node) + this->distance(category, node, my::depot);

    this->currSol.routesNode[node] = this->Ksum;
    this->currSol.positionNode[node] = 1;

    this->currSol.Category.push_back(category);
    this->currSol.Soltours.push_back({my::depot, node, my::depot});
    this->currSol.Costours.push_back({distanceTour, this->distance(category, my::depot, node), distanceTour});

    this->currSol.Charge.emplace_back(3, this->demand[node]);

    this->currSol.Operational.push_back(distanceTour + this->D[category]);

    this->currSol.SolObjcost += distanceTour + this->F[category];
    this->Ksum++;
    this->K[category]++;
  }

  void make_new_tour(const int &category) {
    ++this->Ksum;
    ++this->K[category];

    this->currSol.Category.push_back(category);

    this->currSol.Soltours.emplace_back(2, my::depot);
    this->currSol.Costours.emplace_back(2, 0);
    this->currSol.Charge.emplace_back(2, 0);

    this->currSol.Operational.push_back(0);

    this->currSol.SolObjcost += this->F[category];
  }

  void erase_empty_tour(const int &tourRemove) {

    if (this->Ksum != this->currSol.get_size()) {
      my::error_function("erase_empty_tour", "AQUI different number of vehicles " + my::itos(this->Ksum) + " " + my::itos(this->currSol.Soltours.size()));
    }

    this->currSol.SolObjcost -= this->F[this->currSol.Category[tourRemove]];
    this->K[this->currSol.Category[tourRemove]]--;

    if (tourRemove != this->Ksum - 1) {
      this->currSol.set_route(tourRemove, this->currSol.get_route(this->Ksum - 1));
      this->currSol.Costours[tourRemove] = this->currSol.get_route_cost_list(this->Ksum - 1);
      this->currSol.Charge[tourRemove] = this->currSol.get_route_load_list(this->Ksum - 1);
      this->currSol.Operational[tourRemove] = this->currSol.Operational[this->Ksum - 1];
      this->currSol.Category[tourRemove] = this->currSol.Category[this->Ksum - 1];

      for (int i = Solver::get_cust_begin(); i < this->currSol.get_route_size(tourRemove) - 1; ++i) {
        this->currSol.routesNode[this->currSol.Soltours[tourRemove][i]] = tourRemove;
      }

    }

    this->currSol.Soltours.pop_back();
    this->currSol.Costours.pop_back();
    this->currSol.Charge.pop_back();
    this->currSol.Operational.pop_back();
    this->currSol.Category.pop_back();

    --this->Ksum;
    this->currSol.routesNode[0] = this->get_cust_end();
  }

  void erase_tour(const int tourRemove) {

    if (tourRemove >= this->currSol.get_size()) {
      return;
    }

    int sizeTour = this->currSol.get_route_size(tourRemove) - Solver::get_cust_begin();

    this->currSol.SolObjcost -= this->currSol.get_route_cost(tourRemove);

    for (int i = 1; i < sizeTour; ++i) {
      const int node = this->currSol.node_at(tourRemove, i);
      this->currSol.positionNode[node] = my::dummy;
      this->currSol.routesNode[node] = my::dummy;
      this->unaffectedNodes.push_back(node);
    }

    this->erase_empty_tour(tourRemove);
  }

  void erase_exceeded_routes(const int &category) {
    auto routes = std::vector<std::vector<double>>(this->K[category], std::vector<double>(2));
    int index = 0u;

    for (int k = 0; k < this->Ksum; ++k) {

      if (this->currSol.Category[k] == category) {

        double maxWeightTour = 0;
        for (auto i = Solver::get_cust_begin(); i < this->currSol.get_route_size(k) - 1; ++i) {

          const int node = this->currSol.node_at(k, i);

          if (maxWeightTour < this->demand[node]) {
            maxWeightTour = this->demand[node];
          }

        }

        routes[index][0] = static_cast<double>(this->currSol.node_at(k, 1));
        routes[index][1] = maxWeightTour;
        ++index;
      }

    }

    std::sort(routes.begin(), routes.begin() + this->K[category], [this](std::vector<double> i, std::vector<double> j) {
      return (i[1] > j[1]);
    });

    int size = this->K[category] - this->Kmax[category];
    for (int k = 1; k <= size; ++k) {
      const int node_index = static_cast<int>(routes[index - k][0]);
      this->erase_tour(this->currSol.get_route_index(node_index));
    }

  }

  int ask_size_path(const int &node, const int &sizePath) {
    int route_index = this->currSol.get_route_index(node);

    if (this->currSol.get_position_node(node) + sizePath - Solver::get_cust_begin() < this->currSol.get_route_size(route_index) - Solver::get_cust_begin()) {
      return sizePath;
    } else {
      int actualSize = this->currSol.get_route_size(route_index) - Solver::get_cust_begin() - this->currSol.get_position_node(node);
      return actualSize;
    }

  }

  double cost_insertion(const int &node, const int &route, const int &position) {
    if (position > 0 && position < this->currSol.get_route_size(route)) {

      int category = this->currSol.Category[route];
      int nodeBefore = this->currSol.node_at(route, position - Solver::get_cust_begin());
      int nodeAfter = this->currSol.node_at(route, position);

      if (route == this->Ksum) {
        nodeBefore = my::depot;
        nodeAfter = my::depot;
      }

      return (this->distance(category, nodeBefore, node) + this->distance(category, node, nodeAfter) - this->distance(category, nodeBefore, nodeAfter));
    } else {
      std::cout << "ERROR cost_insertion; pos: " << position << " sizeRoute: " << this->currSol.Soltours[route].size() << std::endl;
      exit(-1);
    }

  }

  void neighborhood_insertion_free_node(const int &node, const int &index, const double &regret) {
    double minCost = std::numeric_limits<float>::max();
    int bestNeighbor = my::dummy;

    const int neighborSize = (index == 1) ? this->granularSizeNodes.get_neighbor(node) : this->maxGranularSize;
    for (int i = 0; i < neighborSize; ++i) {

      const int neighbor = this->neighborsDistance.at(node, i);

      if (neighbor == node || neighbor == my::depot) { continue; }

      const int posNeigh = this->currSol.positionNode[neighbor];
      const int routeNeigh = this->currSol.routesNode[neighbor];

      if (routeNeigh == my::dummy) {
        continue;
      }

      const int categoryNeighbor = this->currSol.Category[routeNeigh];

      if (this->demand[node] + this->currSol.Charge[routeNeigh][0] <= this->Q[categoryNeighbor]) {
        int preNeigh = this->currSol.Soltours[routeNeigh][posNeigh - 1];
        int postNeigh = this->currSol.Soltours[routeNeigh][posNeigh];

        if (routeNeigh == this->Ksum) {
          preNeigh = my::depot;
          postNeigh = my::depot;
        }

        double curCost = this->distance(categoryNeighbor, preNeigh, node) + this->distance(categoryNeighbor, node, postNeigh) - this->distance(categoryNeighbor, preNeigh, postNeigh);

        if (curCost < minCost || this->regretSelect(this->randEngine) < regret) {
          bestNeighbor = neighbor;
          minCost = curCost;
        }

      }

    }

    if (bestNeighbor != my::dummy) {
      this->insert_node(node, bestNeighbor);
      this->customersCache.put(node);
    } else {
      this->make_new_tour(node, CVRP);
      this->customersCache.put(node);
    }

  }

  void best_insertion_free_node(const int &node, const double &regret) {
    double minCost = std::numeric_limits<float>::max();
    double curCost = my::dummy;
    int bestPos = my::dummy;
    int bestRoute = my::dummy;

    for (int k = 0; k < this->Ksum; ++k) {

      const int category = this->currSol.Category[k];

      if (this->demand[node] + this->currSol.Charge[k][0] <= this->Q[category]) {

        for (int i = Solver::get_cust_begin(); i < this->currSol.get_route_size(k); ++i) {

          curCost = this->cost_insertion(node, k, i);

          if (curCost < minCost || this->regretSelect(this->randEngine) < regret) {
            bestPos = i;
            bestRoute = k;
            minCost = curCost;
          }

        }

      }

    }

    for (int c = 0; c < this->C; ++c) {

      if (this->K[c] < this->Kmax[c] && this->demand[node] <= this->Q[c]) {

        curCost = this->distance(c, 0, node) + this->distance(c, node, 0) + this->F[c];

        if (curCost < minCost) {
          bestPos = 0;
          minCost = curCost;
          bestRoute = c;
        }

      }

    }

    if (bestPos > 0) {
      this->currSol.routesNode[node] = bestRoute;
      this->insert_node(node, this->currSol.node_at(bestRoute, bestPos - Solver::get_cust_begin()));
      this->customersCache.put(node);
    } else {  // a new route is created
      this->make_new_tour(node, CVRP);
      this->customersCache.put(node);
    }

  }

  bool repair_all_unaffected_nodes(const int &repairIndex, const double &regretIndex) {

    if (this->unaffectedNodes.empty()) {
      return true;
    }

    switch (this->repairOrder(this->randEngine)) {
      case 0:
        std::shuffle(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), this->randEngine);
        break;
      case 1:
        std::sort(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), [this](int i, int j) {
          return this->demand[i] > this->demand[j];
        });
        break;
      case 2:
        std::sort(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), [this](int i, int j) {
          return this->distance(this->currSol.Category[this->currSol.routesNode[i]], my::depot, i) > this->distance(this->currSol.Category[this->currSol.routesNode[j]], my::depot, j);
        });
        break;
      case 3:
        std::sort(this->unaffectedNodes.begin(), this->unaffectedNodes.end(), [this](int i, int j) {
          return this->distance(this->currSol.Category[this->currSol.routesNode[i]], my::depot, i) < this->distance(this->currSol.Category[this->currSol.routesNode[j]], my::depot, j);
        });
        break;

    }

    for (int &node : this->unaffectedNodes) {

      switch (repairIndex) {
        case 1:
          this->best_insertion_free_node(node, regretIndex);
          break;
        case 2:
          this->neighborhood_insertion_free_node(node, 2, regretIndex);
          break;
        case 3:
          this->neighborhood_insertion_free_node(node, 1, regretIndex);
          break;
      }

    }

    for (int &node : this->unaffectedNodes) {
      if (this->currSol.routesNode[node] == my::dummy) {
        return false;
      }
    }

    this->unaffectedNodes.clear();
    return true;
  }

  void build_one_customer_routes(std::vector<bool> &beginTour, std::vector<bool> &endTour, const int &category) {

    for (int node = Solver::get_cust_begin(); node < this->get_cust_end(); ++node) {

      if (this->currSol.routesNode[node] == my::dummy &&
          this->demand[node] < this->Q[category] &&
          this->distance(category, my::depot, node) < std::numeric_limits<float>::max()) {
        this->make_new_tour(node, category);
        beginTour[node] = true;
        endTour[node] = true;
      } else {
        beginTour[node] = false;
        endTour[node] = false;
      }

    }

  }

  bool ij_eval_2opt(const int &evalNode, const int &evalNeighbor, const int &tourNeighbor, double &costEval) {
    const int posNode = this->currSol.positionNode[evalNode];
    const int tourNode = this->currSol.routesNode[evalNode];
    const int postNode = this->currSol.node_at(tourNode, posNode + 1);
    const int tourSize = this->currSol.get_route_size(tourNode);

    if (this->currSol.Category[tourNode] != this->currSol.Category[tourNeighbor]) {
      return false;
    }

    const int posNeighbor = this->currSol.positionNode[evalNeighbor];

    const int preNeigh = this->currSol.node_at(tourNeighbor, posNeighbor - 1);

    const int category = this->currSol.Category[tourNode];

    const double current = this->length_path(tourNode, category, posNode, posNeighbor);
    const double possible = this->distance(category, evalNode, preNeigh) + this->length_reverse_path(tourNeighbor, category, posNeighbor - 1, posNode + 1) + this->distance(category, postNode, evalNeighbor);

    costEval = possible - current;

    return (tourNode == tourNeighbor && posNeighbor > posNode + 3 && posNeighbor < tourSize);
  }

  void search_2opt(int &bestNode, int &bestNeighbor, double &minResult) {
    bestNode = my::dummy;
    bestNeighbor = my::dummy;
    minResult = 0.0f;

    int previousBestNode = my::dummy;

    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      const int posNode = this->currSol.positionNode[node];
      const int tour = this->currSol.routesNode[node];
      const int tourSize = this->currSol.get_route_size(tour);

      if (tourSize > 4 && posNode + 3 < tourSize) {

        for (int j = posNode + 3; j < tourSize; ++j) {

          const int neighbor = this->currSol.node_at(tour, j);

          if (!this->taboo.is_any(node, neighbor)) {
            const int tourNeigh = this->currSol.routesNode[neighbor];
            double result = 0.0f;

            if (tour == tourNeigh && this->ij_eval_2opt(node, neighbor, tourNeigh, result)) {
              if (result < minResult - EPSILON) {

                minResult = result;
                bestNode = node;
                bestNeighbor = neighbor;
                if (previousBestNode != my::dummy && previousBestNode != node) {
                  this->granularSizeNodes.put(previousBestNode);
                }
              } else {
                previousBestNode = node;
              }
            }

          }

        }

      }

      if (previousBestNode != my::dummy && previousBestNode != bestNode) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

  }

  double cost_eval_realloc(const int &evalNode, const int &evalNeighbor, const int &tourNeighbor) {
    const int posNode = this->currSol.positionNode[evalNode];
    const int tourNode = this->currSol.routesNode[evalNode];

    const int preNode = this->currSol.node_at(tourNode, posNode - 1);
    const int postNode = this->currSol.node_at(tourNode, posNode + 1);

    const int categoryNode = this->currSol.Category[tourNode];

    const int categoryNeighbor = this->currSol.Category[tourNeighbor];
    const int posNeighbor = this->currSol.positionNode[evalNeighbor];

    const int postNeighbor = (evalNeighbor != my::depot) ? this->currSol.node_at(tourNeighbor, posNeighbor + 1) : this->currSol.node_at(tourNeighbor, 1);

    double addInsert = this->distance(categoryNeighbor, evalNeighbor, evalNode) + this->distance(categoryNeighbor, evalNode, postNeighbor) - this->distance(categoryNeighbor, evalNeighbor, postNeighbor);
    double addRemove = this->distance(categoryNode, preNode, postNode) - this->distance(categoryNode, preNode, evalNode) - this->distance(categoryNode, evalNode, postNode);

    return (addRemove + addInsert);
  }

  bool ij_eval_realloc(const int &evalNode, const int &evalNeighbor, const int &tourNeighbor, double &costEval) {
    const int posNode = this->currSol.positionNode[evalNode];
    const int tourNode = this->currSol.routesNode[evalNode];
    const int categoryNeighbor = this->currSol.Category[tourNeighbor];
    const int posNeighbor = this->currSol.positionNode[evalNeighbor];

    const int postNeighbor = (evalNeighbor != my::depot) ? this->currSol.node_at(tourNeighbor, posNeighbor + 1) : this->currSol.node_at(tourNeighbor, 1);

    double addInsert = 0.0f;
    if (evalNeighbor != my::depot && tourNode != tourNeighbor) {
      addInsert = this->distance(categoryNeighbor, evalNeighbor, evalNode) + this->distance(categoryNeighbor, evalNode, postNeighbor) - this->distance(categoryNeighbor, evalNeighbor, postNeighbor);
    }

    costEval = this->cost_eval_realloc(evalNode, evalNeighbor, tourNeighbor);

    if (evalNeighbor != my::depot) {
      return ((tourNode != tourNeighbor && this->currSol.get_route_load(tourNeighbor) + this->demand[evalNode] <= this->Q[categoryNeighbor] && this->currSol.get_route_cost(tourNeighbor) + addInsert + this->D[categoryNeighbor] * (this->currSol.get_route_size(tourNeighbor) - 1.0) >= this->L[categoryNeighbor]) || (tourNeighbor == tourNode && (posNeighbor > posNode || posNeighbor < posNode - 1)));
    } else {
      return ((tourNode != tourNeighbor && this->currSol.get_route_load(tourNeighbor) + this->demand[evalNode] <= this->Q[categoryNeighbor] && this->currSol.get_route_cost(tourNeighbor) + costEval + this->D[this->currSol.Category[tourNeighbor]] * (this->currSol.get_route_size(tourNeighbor) - 1.0) <= this->L[this->currSol.Category[tourNeighbor]]) || (tourNode == tourNeighbor));
    }

  }

  void search_realloc(int &bestNode, int &bestNeighbor, int &tourNeighbor, double &minResult) {
    bestNode = my::dummy;
    minResult = 0.0f;

    int previousBestNode = my::dummy;

    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node); ++j) {

        const int neighbor = this->neighborsDistance.at(node, j);

        if (node == neighbor || neighbor == my::depot ) { continue; }

        if (this->taboo.is_any(node, neighbor)) { continue; }

        const int tourNeighborTmp = this->currSol.routesNode[neighbor];
        double resultTmp = 0.0f;

        if (this->ij_eval_realloc(node, neighbor, tourNeighborTmp, resultTmp)) {
           if (resultTmp < minResult) {
             minResult = resultTmp;
             bestNode = node;
             bestNeighbor = neighbor;
             tourNeighbor = tourNeighborTmp;

             if (previousBestNode != my::dummy && previousBestNode != node) {
               this->granularSizeNodes.put(previousBestNode);
             }

           } else {
             previousBestNode = node;
           }
        }

      }

      for (int k = 0; k < this->Ksum; ++k) {

        double resultTmp = 0.0f;

        if ((this->currSol.routesNode[node] != k || (this->currSol.routesNode[node] == k && this->currSol.positionNode[node] != 1)) &&
            this->ij_eval_realloc(node, my::depot, k, resultTmp)) {

          if (resultTmp < minResult) {
            minResult = resultTmp;
            bestNode = node;
            bestNeighbor = my::depot;
            tourNeighbor = k;
            if (previousBestNode != my::dummy && previousBestNode != node) {
              this->granularSizeNodes.put(previousBestNode);
            }
          } else {
            previousBestNode = node;
          }

        }

        ++k;
      }

      if (previousBestNode != my::dummy && previousBestNode != bestNode) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

  }

  bool ij_eval_swap(const int &evalNode, const int &evalNeighbor, const int &tourNeighbor, double &costEval) {
    const int tourNode = this->currSol.routesNode[evalNode];
    const int categoryNode = this->currSol.Category[tourNode];
    const int posNode = this->currSol.positionNode[evalNode];

    const int preNode = this->currSol.node_at(tourNode, posNode - 1);
    const int nextNode = this->currSol.node_at(tourNode, posNode + 1);

    const int categoryNeighbor = this->currSol.Category[tourNeighbor];
    const int posNeigh = this->currSol.positionNode[evalNeighbor];

    const int preNeigh = this->currSol.node_at(tourNeighbor, posNeigh - 1);
    const int nextNeigh = this->currSol.node_at(tourNeighbor, posNeigh + 1);

    double addInsert = - this->distance(categoryNeighbor, preNeigh, evalNeighbor) - this->distance(categoryNeighbor, evalNeighbor, nextNeigh) + this->distance(categoryNeighbor, preNeigh, evalNode) + this->distance(categoryNeighbor, evalNode, nextNeigh);
    double addRemove = - this->distance(categoryNode, preNode, evalNode) - this->distance(categoryNode, evalNode, nextNode) + this->distance(categoryNode, preNode, evalNeighbor) + this->distance(categoryNode, evalNeighbor, nextNode);

    costEval = addInsert + addRemove;

    return ((tourNeighbor != tourNode && this->currSol.Charge[tourNode][0] - this->demand[evalNode] + this->demand[evalNeighbor] <= this->Q[this->currSol.Category[tourNode]] && this->currSol.Charge[tourNeighbor][0] + this->demand[evalNode] - this->demand[evalNeighbor] <= this->Q[this->currSol.Category[tourNeighbor]]) || (tourNeighbor == tourNode && evalNode != evalNeighbor && preNode != evalNeighbor && nextNode != evalNeighbor));
  }

  void search_swap(int &bestNode, int &bestNeighbor, double &minResult) {
    bestNode = my::dummy;
    bestNeighbor = my::dummy;
    minResult = 0.0f;

    int previousBestNode = my::dummy;
    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node); ++j) {

        const int neighbor = this->neighborsDistance.at(node, j);

        if (neighbor == node || neighbor == my::depot) { continue; }

        if (this->taboo.is_any(node, neighbor)) { continue; }

        const int tourNeighbor = this->currSol.routesNode[neighbor];
        double eval = 0.0f;

        if (this->ij_eval_swap(node, neighbor, tourNeighbor, eval)) {
          if (eval < minResult) {
            minResult = eval;
            bestNode = node;
            bestNeighbor = neighbor;

            if (previousBestNode != my::dummy && previousBestNode != node) {
              this->granularSizeNodes.put(previousBestNode);
            }

          }  else {
            previousBestNode = node;
          }
        }

      }

      if (previousBestNode != my::dummy && previousBestNode != bestNode) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

  }

  void search_cross_exchange(std::vector<EvaluationMemory> &bestMoves, double &minResult) {
    int bestNode = my::dummy;
    int bestNeighbor = my::dummy;
    minResult = 0.0f;

    int previousBestNode = my::dummy;
    bool swapNode, swapNeighbor;

    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      const int tourNode = this->currSol.routesNode[node];
      const int posNode = this->currSol.positionNode[node];

      if (posNode >= this->currSol.get_route_size(tourNode) - 2) {
        continue;
      }

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node); ++j) {

        const int neighbor = this->neighborsDistance.at(node, j);

        if (neighbor == node || neighbor == my::depot) { continue; }

        if (this->taboo.is_any(node, neighbor)) {
          continue;
        }

        const int tourNeighbor = this->currSol.routesNode[neighbor];
        const int possNeighbor = this->currSol.positionNode[neighbor];

        if (tourNode == tourNeighbor ||
            possNeighbor >= this->currSol.get_route_size(tourNeighbor) - 2) {
          continue;
        }

        const int categoryNode = this->currSol.Category[tourNode];

        const int preNode = this->currSol.node_at(tourNode, posNode - 1);
        const int tailNode = this->currSol.node_at(tourNode, posNode + 1);
        const int postNode = this->currSol.node_at(tourNode, posNode + 2);

        double chargeNodes = this->demand[node] + this->demand[tailNode];

        const int categoryNeighbor = this->currSol.Category[tourNeighbor];
        const int preNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor - 1];
        const int tailNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor + 1];
        const int postNeighbor = this->currSol.Soltours[tourNeighbor][possNeighbor + 2];

        double chargeNeighbors = this->demand[neighbor] + this->demand[tailNeighbor];
        double newChargeTour = this->currSol.Charge[tourNode][0] - chargeNodes + chargeNeighbors;
        double newChargeTourNeighbor = this->currSol.Charge[tourNeighbor][0] + chargeNodes - chargeNeighbors;

        if (newChargeTour > this->Q[this->currSol.Category[tourNode]] || newChargeTourNeighbor > this->Q[this->currSol.Category[tourNeighbor]]) {
          continue;
        }

        bool swapNodeTmp, swapNeighborTmp;

        double tmpResult = - this->distance(categoryNode, preNode, node) - this->distance(categoryNode, tailNode, postNode) - this->distance(categoryNeighbor, preNeighbor, neighbor) - this->distance(categoryNeighbor, tailNeighbor, postNeighbor);
        double resLeft = this->distance(categoryNeighbor, preNeighbor, node) + this->distance(categoryNeighbor, tailNode, postNeighbor) - this->distance(categoryNode, node, tailNode) + this->distance(categoryNeighbor, node, tailNode);
        double resRight = this->distance(categoryNeighbor, preNeighbor, tailNode) + this->distance(categoryNeighbor, tailNode, node) - this->distance(categoryNode, node, tailNode) + this->distance(categoryNeighbor, node, postNeighbor);

        if (resLeft < resRight) {
          swapNodeTmp = false;
          tmpResult += resLeft;
        } else {
          swapNodeTmp = true;
          tmpResult += resRight;
        }

        resLeft = this->distance(categoryNode, preNode, neighbor) + this->distance(categoryNode, neighbor, tailNeighbor) - this->distance(categoryNeighbor, neighbor, tailNeighbor) + this->distance(categoryNode, tailNeighbor, postNode);  // no swap
        resRight = this->distance(categoryNode, preNode, tailNeighbor) + this->distance(categoryNode, tailNeighbor, neighbor) - this->distance(categoryNeighbor, neighbor, tailNeighbor) + this->distance(categoryNode, neighbor, postNode);  // swap

        if (resLeft < resRight) {
          swapNeighborTmp = false;
          tmpResult += resLeft;
        } else {
          swapNeighborTmp = true;
          tmpResult += resRight;
        }

        double newLenghthTourNode = this->currSol.Costours[tourNode][0] + this->distance(categoryNode, preNode, neighbor) + this->distance(categoryNode, neighbor, tailNeighbor) + this->distance(categoryNode, tailNeighbor, postNode) - this->distance(categoryNode, preNode, node) - this->distance(categoryNode, node, tailNode) - this->distance(categoryNode, tailNode, postNode);
        double newLenghthTourNeighbor = this->currSol.Costours[tourNeighbor][0] + this->distance(categoryNeighbor, preNeighbor, node) + this->distance(categoryNeighbor, node, tailNode) + this->distance(categoryNeighbor, tailNode, postNeighbor) - this->distance(categoryNeighbor, preNeighbor, neighbor) - this->distance(categoryNeighbor, neighbor, tailNeighbor) - this->distance(categoryNeighbor, tailNeighbor, postNeighbor);

        if (newLenghthTourNode <= this->L[this->currSol.Category[tourNode]] - EPSILON &&
            newLenghthTourNeighbor <= this->L[this->currSol.Category[tourNeighbor]] - EPSILON) {

          if (tmpResult < minResult - EPSILON) {
            bestNode = node;
            bestNeighbor = neighbor;
            minResult = tmpResult;
            swapNode = swapNodeTmp;
            swapNeighbor = swapNeighborTmp;
            if (previousBestNode != my::dummy && previousBestNode != node) {
              this->granularSizeNodes.put(previousBestNode);
            }
          } else {
            previousBestNode = node;
          }

        }

      }

      if (previousBestNode != my::dummy && previousBestNode != bestNode) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

    EvaluationMemory move;

    move.node = my::dummy;
    if (bestNode != my::dummy) {
      bestMoves = this->cross_exchange_to_eval_move(bestNode, bestNeighbor, swapNode, swapNeighbor);
    } else {
      bestMoves = {move};
    }

  }

  void search_2optstar(std::vector<EvaluationMemory> &bestMoves, double &minResult) {
    int bestNode = my::dummy;
    int bestNeighbor = my::dummy;
    minResult = 0.0f;

    int previousBestNode = my::dummy;
    bool split = false;
    bool swapNode = false;
    bool swapNextNode = false;

    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node); ++j) {

        const int tourNode = this->currSol.routesNode[node];
        const int neighbor = this->neighborsDistance.at(node, j);

        if (neighbor == node || neighbor == my::depot) { continue; }

        const int tourNeighbor = this->currSol.routesNode[neighbor];

        if (tourNeighbor == tourNode || this->taboo.is_any(node, neighbor)) {
          continue;
        }

        const int posNode = this->currSol.positionNode[node];
        const int nextNode = this->currSol.node_at(tourNode, posNode + 1);
        const int categoryNode = this->currSol.Category[tourNode];

        double chargePrevious = 0;

        const double chargeNodes = this->currSol.Charge[tourNode][posNode];
        const double chargeNextNodes = this->currSol.Charge[tourNode][0] - chargeNodes;

        const int tailNeighbor = this->currSol.node_at(tourNeighbor, this->currSol.get_route_size(tourNeighbor) - 2);
        const int categoryNeighbor = this->currSol.Category[tourNeighbor];
        const int posNeighbor = this->currSol.positionNode[neighbor];
        const int prevNeigh = this->currSol.node_at(tourNeighbor, posNeighbor - 1);

        const int headNeighbor = (prevNeigh != 0) ? this->currSol.Soltours[tourNeighbor][1] : my::depot;

        double costReversePreviouscNeighbor = my::dummy;
        double costReversePreviouscNode = (prevNeigh != 0) ? this->length_reverse_path(tourNeighbor, categoryNode, posNeighbor - 1, 0) - this->length_path(tourNeighbor, categoryNeighbor, 0, posNeighbor - 1) : 0;

        if (categoryNode == categoryNeighbor) {
          costReversePreviouscNeighbor = costReversePreviouscNode;
        } else {
          if (prevNeigh != 0) {
            costReversePreviouscNeighbor = this->length_reverse_path(tourNeighbor, categoryNeighbor, posNeighbor - 1, 0) - this->length_path(tourNeighbor, categoryNeighbor, 0, posNeighbor - 1);
          } else {
            costReversePreviouscNode = 0;
          }
        }

        const double costReverseNeighborcNode = this->length_reverse_path(tourNeighbor, categoryNode, this->currSol.get_route_size(tourNeighbor) - 1, posNeighbor) - this->length_path(tourNeighbor, categoryNeighbor, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 1);
        const double costReverseNeighborcNeighbor = (categoryNode == categoryNeighbor) ? costReverseNeighborcNode : this->length_reverse_path(tourNeighbor, categoryNeighbor, this->currSol.get_route_size(tourNeighbor) - 1, posNeighbor) - this->length_path(tourNeighbor, categoryNeighbor, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 1);

        chargePrevious = (posNeighbor != 1) ? this->currSol.Charge[tourNeighbor][posNeighbor - 1] : 0;

        double chargeNeighbors = this->currSol.Charge[tourNeighbor][0] - chargePrevious;

        bool splitTmp = false;
        bool swapNextNode1 = false;
        bool swapNextNode2 = false;
        bool swapNode1 = false;
        bool swapNode2 = false;

        double resSplit1, resSplit2, resultTmp;

        if (prevNeigh != 0 || nextNode != 0) {

          if (chargeNodes + chargeNeighbors <= this->Q[categoryNode] && chargeNextNodes + chargePrevious <= this->Q[categoryNeighbor]) {

            double resLeft = this->distance(categoryNode, node, neighbor);
            double resRight = this->distance(categoryNode, node, tailNeighbor) - this->distance(categoryNeighbor, tailNeighbor, my::depot) + this->distance(categoryNode, neighbor, my::depot) + costReverseNeighborcNode;

            if (resLeft < resRight) {
              swapNode1 = false;
              resSplit1 = resLeft;
            } else {
              swapNode1 = true;
              resSplit1 = resRight;
            }

            if (this->currSol.get_route_cost_at(tourNode, posNode) + resSplit1 + this->length_path(tourNeighbor, categoryNode, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 1) + this->D[categoryNode] * (posNode + this->ask_size_path(neighbor, this->get_cust_end())) > this->L[categoryNode]) {
              resSplit1 = INT_MAX;
            } else {
              resLeft = this->distance(categoryNeighbor, prevNeigh, nextNode);
              resRight = this->distance(categoryNeighbor, headNeighbor, nextNode) + this->distance(categoryNeighbor, my::depot, prevNeigh) - this->distance(categoryNeighbor, my::depot, headNeighbor) + costReversePreviouscNeighbor;
              if (resLeft < resRight) {
                swapNextNode1 = false;
                resSplit1 += resLeft;
              } else {
                swapNextNode1 = true;
                resSplit1 += resRight;
              }

              if (this->currSol.Costours[tourNeighbor][posNeighbor - 1] + std::min(resLeft, resRight) +
                  this->length_path(tourNode, categoryNode, posNode + 1, this->currSol.Soltours[tourNode].size() - 1) +
                  this->D[categoryNeighbor] * (posNeighbor - 1 + this->ask_size_path(nextNode, this->get_cust_end())) > this->L[categoryNeighbor]) {
                resSplit1 = INT_MAX;
              }

            }

          } else {

            resSplit1 = INT_MAX;

          }

          if (chargeNodes + chargePrevious <= this->Q[categoryNode] && chargeNextNodes + chargeNeighbors <= this->Q[categoryNeighbor]) {

            double resLeft = this->distance(categoryNode, node, prevNeigh) + costReversePreviouscNode;
            double resRight = this->distance(categoryNode, node, headNeighbor) + this->distance(categoryNode, prevNeigh, my::depot) - this->distance(categoryNode, my::depot, headNeighbor);

            if (resLeft < resRight) {
              swapNode2 = false;
              resSplit2 = resLeft;
            } else {
              swapNode2 = true;
              resSplit2 = resRight;
            }

            if (this->currSol.Costours[tourNode][posNode] + resSplit2 + this->currSol.Costours[tourNeighbor][posNeighbor - 1] + this->D[categoryNode] * (posNeighbor - 1 + posNode) > this->L[categoryNode]) {
              resSplit2 = INT_MAX;
            } else {
              resLeft = this->distance(categoryNeighbor, nextNode, neighbor) + costReverseNeighborcNeighbor;
              resRight = this->distance(categoryNeighbor, nextNode, tailNeighbor) + this->distance(categoryNeighbor, my::depot, neighbor) - this->distance(categoryNeighbor, tailNeighbor, my::depot);

              if (resLeft < resRight) {
                swapNextNode2 = false;
                resSplit2 += resLeft;
              } else {
                swapNextNode2 = true;
                resSplit2 += resRight;
              }

              if (this->length_path(tourNode, categoryNode, posNode + 1, this->currSol.get_route_size(tourNode) - 1) + std::min(resLeft, resRight) + this->length_path(tourNeighbor, categoryNeighbor, posNeighbor, this->currSol.get_route_size(tourNeighbor) - 1) + this->D[categoryNeighbor] * (this->ask_size_path(nextNode, this->get_cust_end()) + this->ask_size_path(neighbor, this->get_cust_end())) > this->L[categoryNeighbor]) {
                resSplit2 = INT_MAX;
              }

            }

          } else {

            resSplit2 = INT_MAX;

          }

          if (resSplit1 < resSplit2) {
            splitTmp = false;
            resultTmp = resSplit1;
          } else {
            splitTmp = true;
            resultTmp = resSplit2;
          }

        } else {

          resultTmp = INT_MAX;

        }

        if (resultTmp != INT_MAX) {
          resultTmp -= (this->distance(categoryNode, node, nextNode) + this->distance(categoryNeighbor, prevNeigh, neighbor));
        }

        if (resultTmp < minResult - EPSILON) {
          minResult = resultTmp;
          split = splitTmp;

          if (!split) {
            swapNode = swapNode1;
            swapNextNode = swapNextNode1;
          } else {
            swapNode = swapNode2;
            swapNextNode = swapNextNode2;
          }

          bestNode = node;
          bestNeighbor = neighbor;

          if (previousBestNode != my::dummy && previousBestNode != node) {
            this->granularSizeNodes.put(previousBestNode);
          }

        } else {
          previousBestNode = node;
        }
      }

      if (previousBestNode != my::dummy && previousBestNode != bestNode) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

    if (bestNode != my::dummy) {
      bestMoves = this->cross_over_to_eval_move(bestNode, bestNeighbor, split, swapNode, swapNextNode);
    } else {
      EvaluationMemory bestMove;
      bestMoves = {bestMove};
    }

  }

  void search_path_move(EvaluationMemory &bestMove, const int &sizeMin, const int &sizePath) {
    EvaluationMemory EvalMove;

    bestMove.node = -1;
    bestMove.result = INT_MAX;

    double minResult = 0.0f;

    int previousBestNode = my::dummy;

    for (int node = this->customersCache.begin(); node != my::dummy; node = this->customersCache.get_next(node)) {

      const int sizePathTmp = this->ask_size_path(node, sizePath);

      if (sizePathTmp < sizeMin) {
        continue;
      }

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node); ++j) {

        const int neighbor = this->neighborsDistance.at(node, j);

        if (neighbor == node || neighbor == my::depot) { continue; }

        if (this->taboo.is_any(node, neighbor)) {
          continue;
        }

        const int posNode = this->currSol.positionNode[node];
        const int tourNode = this->currSol.routesNode[node];
        const int categoryNode = this->currSol.Category[tourNode];
        const int posTail = posNode + sizePathTmp - 1;
        const int tailNode = this->currSol.Soltours[tourNode][posTail];

        double chargeNodes = 0;

        for (int i = 0; i < sizePathTmp; i++) {
          chargeNodes += this->demand[this->currSol.Soltours[tourNode][posNode + i]];
        }

        const int posNeighbor = this->currSol.positionNode[neighbor];
        const int tourNeighbor = this->currSol.routesNode[neighbor];
        const int categoryNeighbor = this->currSol.Category[tourNeighbor];

        if (tourNeighbor == tourNode || this->currSol.Charge[tourNeighbor][0] + chargeNodes > this->Q[categoryNeighbor]) {
          continue;
        }

        const double costReversePath = this->length_reverse_path(tourNode, categoryNeighbor, posTail, posNode);
        const double costPath = this->length_path(tourNode, categoryNeighbor, posNode, posTail);

        const int postNeighbor = this->currSol.Soltours[tourNeighbor][posNeighbor + 1];

        const double removeTmp = -this->length_path(tourNode, categoryNode, posNode - 1,posNode + sizePathTmp) + this->distance(categoryNode, this->currSol.Soltours[tourNode][posNode - 1], this->currSol.Soltours[tourNode][posNode + sizePathTmp]);

        double d0 = this->distance(categoryNeighbor, neighbor, node) + costPath + this->distance(categoryNeighbor, tailNode, postNeighbor) - this->distance(categoryNeighbor, neighbor, postNeighbor);
        double d2 = this->distance(categoryNeighbor, neighbor, tailNode) + costReversePath + this->distance(categoryNeighbor, node, postNeighbor) - this->distance(categoryNeighbor, neighbor, postNeighbor);

        double insertTmp = (d0 < d2) ? d0 : d2;

        int wayTmp = (d0 < d2) ? 0 : 2;

        if (posNeighbor == 1) {
          double d1 = this->distance(categoryNeighbor, my::depot, node) + costPath + this->distance(categoryNeighbor, tailNode, neighbor) - this->distance(categoryNeighbor, my::depot, neighbor);
          double d3 = this->distance(categoryNeighbor, my::depot, tailNode) + costReversePath + this->distance(categoryNeighbor, node, neighbor) - this->distance(categoryNeighbor, my::depot, neighbor);

          if (d1 < insertTmp) {
            insertTmp = d1;
            wayTmp = 1;
          }

          if (d3 < insertTmp) {
            insertTmp = d3;
            wayTmp = 3;
          }

        }

        double resultsTmp = insertTmp + removeTmp;

        if (resultsTmp < minResult) {
          bestMove.node = node;
          bestMove.neighbor = neighbor;
          bestMove.longPathMove = sizePathTmp;
          bestMove.addRemove = removeTmp;
          bestMove.addInsert = insertTmp;
          bestMove.way = wayTmp;
          bestMove.result = resultsTmp;
          bestMove.charge = this->currSol.Charge[tourNeighbor][0] + chargeNodes;
          minResult = resultsTmp;

          if (previousBestNode != my::dummy && previousBestNode != node) {
            this->granularSizeNodes.put(previousBestNode);
          }

        } else {
          previousBestNode = node;
        }

      }

      if (previousBestNode != my::dummy && previousBestNode != bestMove.node) {
        this->granularSizeNodes.put(previousBestNode);
      }

    }

  }

  void search_ec_realloc(double &minResult, int &bestNode, Solution &bestSolution) {  // TODO need special treatment
    const Solution oldSolution = this->currSol;
    const double oldCost = this->currSol.SolObjcost;

    bestSolution = this->currSol;
    minResult = 0.0f;
    bestNode = my::dummy;

    bool stop = false;

    for (int node = this->customersCache.begin(); node != my::dummy && !stop; node = this->customersCache.get_next(node)) {

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node) && !stop; ++j) {

        int infeasible_tour = my::dummy;
        int neighbor = this->neighborsDistance.at(node, j);

        if (neighbor == node || neighbor == my::depot) { continue; }

        if (this->taboo.is_any(node, neighbor)) {
          continue;
        }

        int tourNeighbor = this->currSol.routesNode[neighbor];

        double resultTmp = 0.0f;
        const bool isPossibleRealloc = this->ij_eval_realloc(node, neighbor, tourNeighbor, resultTmp);

        if (isPossibleRealloc) {

          if (resultTmp < minResult)  {
            this->execute_reallocation(node, neighbor, tourNeighbor);
            this->customersCache.put(node);
            this->customersCache.put(neighbor);
            goto __skip__;
          } else {
            continue;
          }

        }

        if (this->currSol.routesNode[node] != tourNeighbor &&
            resultTmp < - EPSILON &&
            this->currSol.Charge[tourNeighbor][0] + this->demand[node] > this->Q[this->currSol.Category[tourNeighbor]] &&
            this->currSol.Charge[tourNeighbor][0] + this->demand[node] <= this->Q[this->currSol.Category[tourNeighbor]] * 2.0) {

          infeasible_tour = this->start_ec_realloc(node, neighbor, tourNeighbor);

        __skip__:
          resultTmp = this->currSol.SolObjcost - oldCost;

          if (resultTmp < minResult - EPSILON && infeasible_tour == my::dummy) {
            minResult = resultTmp;
            bestSolution = this->currSol;
            bestNode = node;
          }

          if (minResult < -EPSILON) {
            stop = true;
          }

          this->currSol = oldSolution;
          this->update_ksum(this->currSol.get_size());
        }

      }

      if (!stop) {

        for (int k = 0; k < this->Ksum; ++k) {

          double resultTmp = 0.0f;

          if (((this->currSol.routesNode[node] == k && this->currSol.positionNode[node] != 1) || this->currSol.routesNode[node] != k) &&
              this->ij_eval_realloc(node, my::depot, k, resultTmp) ) {

            if (resultTmp < minResult) {
              this->execute_reallocation(node, my::depot, k);
              this->customersCache.put(node);
              minResult = resultTmp;
              bestSolution = this->currSol;
              bestNode = node;

              if (minResult < -EPSILON) {
                stop = true;
              }

            }

          }

        }

      }

      if (!stop){
        this->granularSizeNodes.put(node);
      }

    }

    this->currSol = oldSolution;
    this->update_ksum(this->currSol.get_size());
  }

  void search_ec_path_move(const int &sizeMin, const int &sizePath, double &minResult, int &bestNode, Solution &bestSolution) {  // TODO need special treatment
    const Solution oldSolution = this->currSol;
    const double oldCost = oldSolution.SolObjcost;

    bestSolution = this->currSol;

    EvaluationMemory move;
    move.node = my::dummy;

    minResult = 0.0f;

    bool stop = false;
    bestNode = my::dummy;

    for (int node = this->customersCache.begin(); node != my::dummy && !stop; node = this->customersCache.get_next(node)) {

      move.longPathMove = this->ask_size_path(node, sizePath);

      if (move.longPathMove < sizeMin) {
        continue;
      }

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(node) && !stop; ++j) {

        move.node = node;
        move.neighbor = this->neighborsDistance.at(node, j);

        if (move.node == move.neighbor|| move.neighbor == my::depot) { continue; }

        if (this->taboo.is_any(move.node, move.neighbor)) {
          continue;
        }

        const int tourNeighbor = this->currSol.routesNode[move.neighbor];
        const int categoryNeighbor = this->currSol.Category[tourNeighbor];

        const int tourNode = this->currSol.routesNode[move.node];
        const int categoryNode = this->currSol.Category[tourNode];

        const int posNode = this->currSol.positionNode[move.node];
        double chargeNodes = 0;

        for (int i = 0; i < move.longPathMove; i++) {
          chargeNodes += this->demand[this->currSol.Soltours[tourNode][posNode + i]];
        }

        if (tourNeighbor != tourNode && this->currSol.Charge[tourNeighbor][0] + chargeNodes <= this->Q[categoryNeighbor]) {

          const int posNeighbor = this->currSol.positionNode[move.neighbor];
          const int posTail = posNode + move.longPathMove - 1;
          const int tailNode = this->currSol.Soltours[tourNode][posTail];

          double costReversePath = this->length_reverse_path(tourNode, categoryNeighbor, posTail, posNode);
          double costPath = this->length_path(tourNode, categoryNeighbor, posNode, posTail);

          const int postNeighbor = this->currSol.Soltours[tourNeighbor][posNeighbor + 1];

          double removeTmp = - this->length_path(tourNode, categoryNode, posNode - 1,posNode + move.longPathMove) + this->distance(categoryNode, this->currSol.Soltours[tourNode][posNode - 1], this->currSol.Soltours[tourNode][posNode + move.longPathMove]);

          double d0 = this->distance(categoryNeighbor, move.neighbor, node) + costPath + this->distance(categoryNeighbor, tailNode, postNeighbor) - this->distance(categoryNeighbor, move.neighbor, postNeighbor);
          double d2 = this->distance(categoryNeighbor, move.neighbor, tailNode) + costReversePath + this->distance(categoryNeighbor, node, postNeighbor) - this->distance(categoryNeighbor, move.neighbor, postNeighbor);

          double insertTmp = (d0 < d2) ? d0 : d2;

          int wayTmp = (d0 < d2) ? 0 : 2;

          if (posNeighbor == 1) {
            double d1 = this->distance(categoryNeighbor, my::depot, node) + costPath + this->distance(categoryNeighbor, tailNode, move.neighbor) - this->distance(categoryNeighbor, my::depot, move.neighbor);
            double d3 = this->distance(categoryNeighbor, my::depot, tailNode) + costReversePath + this->distance(categoryNeighbor, node, move.neighbor) - this->distance(categoryNeighbor, my::depot, move.neighbor);

            if (d1 < insertTmp) {
              insertTmp = d1;
              wayTmp = 1;
            }

            if (d3 < insertTmp) {
              insertTmp = d3;
              wayTmp = 3;
            }

          }

          double resultTmp = insertTmp + removeTmp;

          if (resultTmp < minResult - EPSILON) {
            move.addRemove = removeTmp;
            move.addInsert = insertTmp;
            move.way = wayTmp;
            move.result = resultTmp;
            move.charge = this->currSol.Charge[tourNeighbor][0] + chargeNodes;

            this->taboo.put(move.node, move.neighbor);
            this->execute_move(move);
            this->customersCache.put(move.node);
            this->customersCache.put(move.neighbor);
            goto __skip__;
          } else {
            continue;
          }
        }

        if (tourNeighbor != tourNode) {

          this->evaluation_movement(move);

          if (move.result < -EPSILON &&
              !move.is_feasible(this->Q[categoryNeighbor], this->L[categoryNeighbor]) &&
              move.is_feasible(this->Q[categoryNeighbor] * 2.0, this->L[categoryNeighbor] * 2.0) &&
              move.node != my::dummy) {

            this->start_ec_path_move(move, sizePath);

          __skip__:
            double resultTmp = this->currSol.SolObjcost - oldCost;

            if (resultTmp < minResult - EPSILON && this->search_unfeasible_tour() == my::dummy) {
              minResult = resultTmp;
              bestSolution = this->currSol;
              bestNode = node;
            }

            if (minResult < -EPSILON) {
              stop = true;
            }

            this->currSol = oldSolution;
            this->update_ksum(this->currSol.get_size());
          }

        }

      }
      if (!stop){
        this->granularSizeNodes.put(node);
      }

    }

    this->currSol = oldSolution;
    this->update_ksum(this->currSol.get_size());
  }

  int start_ec_realloc(int node, int neighbor, int tourNeighbor) {
    int tour = tourNeighbor;

    if (node != my::dummy) {
      this->execute_reallocation(node, neighbor, tourNeighbor);

      tour = tourNeighbor;
      unsigned int iter = 0u;

      while (tour != my::dummy && node != my::dummy && iter < this->params.get_max_chain_iteration()) {
        this->search_chain_tour_realloc(node, neighbor, tour, tourNeighbor);

        if (node != my::dummy) {
          this->taboo.put(node, neighbor);
          this->execute_reallocation(node, neighbor, tourNeighbor);
          this->customersCache.put(node);
          this->customersCache.put(neighbor);
        }

        if (tour != my::dummy && this->currSol.feasible_tour(tour, this->Q[this->currSol.Category[tour]], this->L[this->currSol.Category[tour]])) {
          tour = this->search_unfeasible_tour();
        }

        ++iter;
      }

    }

    return tour;
  }

  void start_ec_path_move(EvaluationMemory move, const int sizePath) {

    this->execute_move(move);

    int tour = this->currSol.routesNode[move.neighbor];
    unsigned int iter = 0u;

    while (tour != my::dummy && move.node != my::dummy && iter < this->params.get_max_chain_iteration()) {

      this->search_chain_tour_path_move(move, tour, sizePath);

      if (move.node != my::dummy) {
        this->taboo.put(move.node, move.neighbor);
        this->execute_move(move);
        this->customersCache.put(move.node);
        this->customersCache.put(move.neighbor);
      }

      if (this->currSol.feasible_tour(tour, this->Q[this->currSol.Category[tour]], this->L[this->currSol.Category[tour]])) {
        tour = this->search_unfeasible_tour();
      }

      ++iter;
    }
  }

  void search_chain_tour_realloc(int &node, int &neighbor, const int &tour, int &tourNeighbor) {
    const int forbiddenNode = node;
    double minResult = 0.0f;
    bool stop = false;

    const int sizeTour = this->currSol.get_route_size(tour) - 1;

    node = my::dummy;

    for (int i = Solver::get_cust_begin(); i < sizeTour && !stop; ++i) {
      const int nodeTmp = this->currSol.node_at(tour, i);

      if (nodeTmp == forbiddenNode) {
        continue;
      }

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(nodeTmp) && !stop; ++j) {
        const int neighborTmp = this->neighborsDistance.at(nodeTmp, j);

        if (neighborTmp == nodeTmp || neighborTmp == my::depot) { continue; }

        const int tourNeighborTmp = this->currSol.routesNode[neighborTmp];

        if (tourNeighborTmp == tour) {
          continue;
        }

        const double result = this->cost_eval_realloc(nodeTmp, neighborTmp, tourNeighborTmp);
        if (this->currSol.get_route_load(tourNeighborTmp) + this->demand[nodeTmp] <= this->Q[this->currSol.Category[tourNeighbor]] &&
            result < minResult) {

          tourNeighbor = tourNeighborTmp;
          node = nodeTmp;
          neighbor = neighborTmp;
          minResult = result;

          if (minResult < -EPSILON) {
            stop = true;
          }

        }

      }

      if (!stop){
        this->granularSizeNodes.put(nodeTmp);
      }

    }

  }

  void search_chain_tour_path_move(EvaluationMemory &bestMove, const int &tour, const int &sizePath) {
    const int forbiddenNode = bestMove.node;
    double minResult = 0.0f;
    bool stop = false;
    const int sizeTour = this->currSol.Soltours[tour].size() - 1;

    EvaluationMemory moveTmp;

    bestMove.node = my::dummy;

    for (int i = Solver::get_cust_begin(); i < sizeTour && !stop; ++i) {

      moveTmp.node = this->currSol.Soltours[tour][i];

      if (moveTmp.node == forbiddenNode) {
        continue;
      }

      for (int j = 0; j < this->granularSizeNodes.get_neighbor(moveTmp.node) && !stop; ++j) {
        moveTmp.neighbor = this->neighborsDistance.at(moveTmp.node, j);
        if (moveTmp.neighbor == moveTmp.node || moveTmp.neighbor == my::depot) { continue; }

        const int tourNeighbor = this->currSol.routesNode[moveTmp.neighbor];

        if (tourNeighbor == tour) {
          continue;
        }

        moveTmp.longPathMove = this->ask_size_path(moveTmp.node, sizePath);

        const int categoryNeighbor = this->currSol.Category[tourNeighbor];

        this->evaluation_movement(moveTmp);

        if (moveTmp.is_feasible(this->Q[categoryNeighbor],this->L[categoryNeighbor]) && moveTmp.result < minResult - EPSILON) {
          bestMove = moveTmp;
          minResult = moveTmp.result;
          if (minResult < -EPSILON) {
            stop = true;
          }
        }

      }

      if (!stop){
        this->granularSizeNodes.put(moveTmp.node);
      }

    }

  }

  void execute_2opt(const int &tour, const int &pos_begin, const int &pos_end, const double &result) {

    std::vector<int> line;
    for (int h = 0; h <= pos_begin; ++h) {
      line.push_back(this->currSol.Soltours[tour][h]);
    }

    const int category = this->currSol.Category[tour];
    double costTmp = (pos_begin != 0) ? this->currSol.Costours[tour][pos_begin] + this->distance(category, this->currSol.Soltours[tour][pos_begin], this->currSol.Soltours[tour][pos_end - 1]) : this->distance(category, this->currSol.Soltours[tour][pos_begin], this->currSol.Soltours[tour][pos_end - 1]);
    double chargeTmp = (pos_begin != 0) ? this->currSol.Charge[tour][pos_begin] : 0;

    int auxposs = pos_begin + Solver::get_cust_begin();
    for (int h = pos_end - Solver::get_cust_begin(); h > pos_begin; h--) {
      line.push_back(this->currSol.Soltours[tour][h]);
      this->currSol.Costours[tour][auxposs] = costTmp;
      costTmp += this->distance(category, this->currSol.Soltours[tour][h], this->currSol.Soltours[tour][h - 1]);
      chargeTmp += this->demand[this->currSol.Soltours[tour][h]];
      this->currSol.Charge[tour][auxposs] = chargeTmp;
      this->currSol.positionNode[this->currSol.Soltours[tour][h]] = auxposs;
      ++auxposs;
    }

    const int tourSize = this->currSol.Soltours[tour].size();
    for (int h = pos_end; h < tourSize; ++h) {
      line.push_back(this->currSol.Soltours[tour][h]);
      this->currSol.Costours[tour][h] = this->currSol.Costours[tour][h] + result;
    }

    this->currSol.Soltours[tour] = line;
    this->currSol.Costours[tour][0] += result;
    this->currSol.SolObjcost += result;
    line.clear();

    this->currSol.Operational[tour] += result;
  }

  void execute_swap(const int &node, const int &neighbor) {
    int tourNode = this->currSol.routesNode[node];
    const int category = this->currSol.Category[this->currSol.routesNode[node]];
    const int preNode = this->currSol.Soltours[tourNode][this->currSol.positionNode[node] - 1];

    if (preNode == my::depot && this->currSol.Soltours[tourNode].size() == 3) {
      tourNode = my::dummy;
    }

    this->remove_node(node);
    this->insert_node(node, neighbor);
    this->remove_node(neighbor);

    this->currSol.routesNode[neighbor] = tourNode;

    (tourNode == my::dummy) ? this->make_new_tour(neighbor, category) : this->insert_node(neighbor, preNode);
  }

  void execute_reallocation(const int &node, const int &neighbor, int &tourNeighbor) {
    const int firstNode = this->currSol.Soltours[tourNeighbor][1];
    this->remove_node(node);
    tourNeighbor = this->currSol.routesNode[firstNode];

    if (neighbor == my::depot) {
      this->currSol.routesNode[node] = tourNeighbor;
    }

    this->insert_node(node, neighbor);
    tourNeighbor = this->currSol.routesNode[firstNode];
  }

  void remove_node(const int &node) {
    const int tour = this->currSol.routesNode[node];
    const int posNode = this->currSol.positionNode[node];
    const int prevNode = this->currSol.Soltours[tour][posNode - 1];
    const int nextNode = this->currSol.Soltours[tour][posNode + 1];
    const int category = this->currSol.Category[tour];
    const double result = -this->distance(category, prevNode, node) - this->distance(category, node, nextNode) + this->distance(category, prevNode, nextNode);

    this->currSol.Soltours[tour].erase(this->currSol.Soltours[tour].begin() + posNode);

    const int sizeTour = this->currSol.Soltours[tour].size();
    if (this->currSol.Soltours[tour].size() == 2) {
      this->erase_empty_tour(tour);
      this->currSol.SolObjcost += result;
      this->currSol.positionNode[node] = my::dummy;
      this->currSol.routesNode[node] = my::dummy;
    } else {
      auto itaux = this->currSol.Charge.begin() + tour;
      (*itaux).erase((*itaux).begin() + posNode);

      itaux = this->currSol.Costours.begin() + tour;
      (*itaux).erase((*itaux).begin() + posNode);

      for (int i = posNode; i < sizeTour; ++i) {
        this->currSol.Charge[tour][i] -= this->demand[node];
        this->currSol.Costours[tour][i] += result;
      }

      this->currSol.Costours[tour][0] += result;
      this->currSol.Charge[tour][0] -= this->demand[node];
      this->currSol.Operational[tour] += (-this->D[this->currSol.Category[tour]] + result);
      this->currSol.SolObjcost += result;

      this->currSol.routesNode[node] = my::dummy;
      this->currSol.positionNode[node] = my::dummy;

      for (int i = posNode; i < sizeTour; ++i) {
        int a = this->currSol.Soltours[tour][i];
        this->currSol.positionNode[a] = i;
      }

    }

  }

  void insert_node(const int node, const int neighbor) {
    int tour, posNode;
    if (neighbor == my::depot) {
      tour = this->currSol.routesNode[node];
      posNode = 1;
    } else {
      tour = this->currSol.routesNode[neighbor];
      posNode = this->currSol.positionNode[neighbor] + 1;
    }

    const int category = this->currSol.Category[tour];
    if (node == my::depot || tour == my::dummy) {
      my::error_function("insert_node", "tour == -1 for neighbor " + my::itos(neighbor) + " or node== " + my::itos(node));
    }

    const int nextNode = this->currSol.Soltours[tour][posNode];
    const double add = - this->distance(category, neighbor, nextNode) + this->distance(category, neighbor, node) + this->distance(category, node, nextNode);

    this->currSol.Soltours[tour].insert(this->currSol.Soltours[tour].begin() + posNode, node);
    this->currSol.routesNode[node] = tour;
    this->currSol.positionNode[node] = posNode;

    const int tourSize = this->currSol.Soltours[tour].size();
    for (int i = posNode; i < tourSize; ++i) {
      this->currSol.positionNode[this->currSol.Soltours[tour][i]] = i;
    }

    auto itaux = this->currSol.Charge.begin() + tour;
    (posNode > 1) ? (*itaux).insert((*itaux).begin() + posNode, this->demand[node] + (*itaux)[posNode - 1])
                  : (*itaux).insert((*itaux).begin() + posNode, this->demand[node]);

    itaux = this->currSol.Costours.begin() + tour;
    (posNode > 1) ? (*itaux).insert((*itaux).begin() + posNode, distance(category, neighbor, node) + (*itaux)[posNode - 1])
                  : (*itaux).insert((*itaux).begin() + posNode, distance(category, neighbor, node));

    for (int i = posNode + 1; i < tourSize; ++i) {
      this->currSol.Charge[tour][i] += this->demand[node];
      this->currSol.Costours[tour][i] += add;
    }

    this->currSol.Costours[tour][0] += add;
    this->currSol.Charge[tour][0] += this->demand[node];
    this->currSol.SolObjcost += add;
    this->currSol.Operational[tour] += this->D[this->currSol.Category[tour]] + add;
  }

  static int get_cust_begin() { return my::depot + 1; }

  int get_cust_end() { return this->N; }

public:
  int maxGranularSize;
  int maxEliteSolSize;

  double regretIndex;
  Parameters &params;				// Problem parameters
  SplitMix64Engine randEngine;

  std::uniform_real_distribution<float> regretSelect;
  std::uniform_int_distribution<int> booleanDist;
  std::uniform_int_distribution<int> repairOrder;
  std::uniform_int_distribution<int> selectingNode;

  // tabu search data structure
  PairBitMatrix taboo;
  BitMatrix skipEval;

  // pool of elite sols
  EliteSol population;

  // neighborhood structure
  std::vector<std::vector<Saving>> savingTable;
  std::vector<int> unaffectedNodes;
  std::vector<int> remainingFreeNodes;
  GranularManagement granularSizeNodes;

  // double primalIntegral;  // TODO debugging
  // double lastTimeSpend;
  // double lastCostFound;
};

#endif //SOLVER_H