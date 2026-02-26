/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#include "Solver.h"

void Solver::run() {
    this->neighborsDistance.resize(this->N, NEIGHBORHOOD_TABLE_SIZE(this->N-1));
    this->maxGranularSize = this->params.get_granular_size(this->stopsPerRoute);
    this->regretIndex = this->params.get_regret_index(this->stopsPerRoute);
    this->taboo = PairBitMatrix(this->params.get_pair_cache_capacity(), this->N);
    this->skipEval = BitMatrix(this->N);  // there is no capacity for TABU during PR iteratively evaluate neighborhood
    this->granularSizeNodes = GranularManagement(this->N, this->maxGranularSize);

    const int initialGranularSize = this->params.get_initial_granular_size(this->stopsPerRoute);
    const double granularGrowthThreshold = this->params.get_threshold_growth(this->stopsPerRoute);

    this->granularSizeNodes.initialization(initialGranularSize);
    this->nbRun = 0;

    this->generate_solutions();

#ifdef GUIDANCE
    std::cout << "Optimization with feature-based guidance...\n";
    double W = this->reffSol.get_s19(this->Q[CVRP]) - this->reffSol.get_s20(this->Q[CVRP]);
#endif
#ifdef GUIDANCE
    int reStartThreshold = static_cast<int>(std::ceil(THRESHOLD_SETTING * W));
#else
    int reStartThreshold = THRESHOLD_SETTING;
#endif

    double nonImprov = 0.;
    int resetPopulation = 0u;

    double currentBestSol = this->globalBest.SolObjcost;
    while (this->stopping_criteria()) {
        const double perturbIndex = nonImprov/granularGrowthThreshold;

        this->neighborhood_improvement(perturbIndex);
        this->gt_path_relinking(perturbIndex);
        ++this->nbRun;

        if (this->globalBest.SolObjcost < currentBestSol - EPSILON) {

            currentBestSol = this->globalBest.SolObjcost;
            nonImprov = 0.0;
            resetPopulation = 0;
            this->granularSizeNodes.granular_reset(initialGranularSize);

        } else {

            nonImprov += 1.0;

            if (nonImprov >= granularGrowthThreshold) {

                this->granularSizeNodes.growth_neighbors(1);
                resetPopulation += static_cast<int>(nonImprov);
                nonImprov = 0.0;

                if (resetPopulation >= reStartThreshold) {

                    this->customersCache.reset();
                    this->localBest.SolObjcost = std::numeric_limits<float>::max();
                    this->clarke_and_wrigth_solver();

#ifdef GUIDANCE
                    W = this->reffSol.get_s19(this->Q[CVRP]) - this->reffSol.get_s20(this->Q[CVRP]);
                    reStartThreshold = static_cast<int>(std::ceil(THRESHOLD_SETTING * W));
#else
                    reStartThreshold = THRESHOLD_SETTING;
#endif

                    resetPopulation = 0;
                }
            }
        }
    }

    this->currSol = this->globalBest;
    this->update_vertices();
    // this->finalize_primal_integral(); // debugging
}

void Solver::generate_solutions() {
    this->initialize_vertices();

    this->savings();
    this->localBest.SolObjcost = std::numeric_limits<float>::max();
    this->clarke_and_wrigth_solver();
}

void Solver::gt_path_relinking(const double &index) {
    bool isSkip = false;

    // start PR whenever diversity is satisfied
    if (this->population.get_size() < 2) { // adjustable
        isSkip = true;
    }

    if (!isSkip) { // transform VRP solution into giant tour
        this->currSol = this->population.get_sol_worst()->second;

        this->update_ksum(this->currSol.get_size());
        this->randomizeConcat();

        this->initialSol = this->currSol;
        this->currSol = this->population.get_begin()->second;

        this->update_ksum(this->currSol.get_size());
        this->randomizeConcat();

        this->guideSol = this->currSol;
    }

    if (!isSkip && this->initialSol.SolObjcost < this->guideSol.SolObjcost) {
        Solution &tmp = this->initialSol;
        this->initialSol = this->guideSol;
        this->guideSol = tmp;
    }

    if (!isSkip && this->initialSol.SolObjcost == this->guideSol.SolObjcost) {
        isSkip = true;  // same solution
    }

    auto restrictedNeighborhood = std::vector<int>();

    int truncated = my::dummy;
    if (!isSkip) {
        const int nMoves = this->initialSol.n_moves(this->guideSol, restrictedNeighborhood);
        truncated = static_cast<int>(std::round(nMoves * this->params.get_truncated_threshold()));  // truncated path relinking
    }

    if (truncated < this->params.get_min_hamming_size()) {
        isSkip = true;  // too similar
    }

    if (!isSkip) {
        switch (this->booleanDist(this->randEngine)) {
            case 0:
                std::shuffle(restrictedNeighborhood.begin(), restrictedNeighborhood.end(), this->randEngine);
                break;
            case 1:
                std::sort(restrictedNeighborhood.begin(), restrictedNeighborhood.end(), [this](int i, int j) {
                    return this->demand[i] > this->demand[j];
                });
                break;
            default:
                isSkip = true;
                break;
        }

        const double regret = this->regretIndex * (REGRET_CONST - index);
        this->iteratively_evaluate_neighbour(
              this->initialSol,
              this->guideSol,
              truncated,
              restrictedNeighborhood,
              regret
        );

    }
}

void Solver::iteratively_evaluate_neighbour(Solution &initial, Solution &guiding, const int &maxIter, std::vector<int> &neighborhood, const double &regret) {
    this->currSol = initial;
    const double splitThreshold = this->initialSol.SolObjcost;
    for (int prIter = 0; prIter < maxIter; ++prIter) {
        int bestNode = my::dummy;
        int bestNeighbor = my::dummy;
        double minResult = 0u;

        for (auto &p : neighborhood) {

            if (!this->skipEval.is_any(p)) {

                const int posNode = this->currSol.get_position_node(p);
                const int posNeighbor = guiding.get_position_node(p);
                const int node = this->currSol.node_at(my::indexGiantTour, posNode);
                const int neighbor = this->currSol.node_at(my::indexGiantTour, posNeighbor);
                const int tourNeighbor = this->currSol.get_route_index(neighbor);

                double eval = std::numeric_limits<float>::max();

                if (this->currSol.get_route_index(node) == tourNeighbor &&
                    this->ij_eval_swap(node, neighbor, tourNeighbor, eval)) {

                    if (eval < minResult - EPSILON) {
                        minResult = eval;
                        bestNode = node;
                        bestNeighbor = neighbor;
                        break;
                    }

                }

            }

        }

        if (bestNode != my::dummy) {
            this->skipEval.put(bestNode);
            this->execute_swap(bestNode, bestNeighbor);

            Solution tmpSol = this->currSol;

            if (tmpSol.SolObjcost < splitThreshold - EPSILON || this->regretSelect(this->randEngine) < regret) {

                const bool isSplited = this->split_giant_tour(my::indexGiantTour);

                if (isSplited) {

                    const bool isInserted = this->population.set_sols(this->currSol);

                    if (isInserted && this->currSol.SolObjcost < this->globalBest.SolObjcost - EPSILON) {
                        this->update_best();
                    }

                  break;

                }

            }

            this->currSol = tmpSol;
            this->update_ksum(this->currSol.get_size());
        }

    }

}

bool Solver::run_level(const int &level) {
    const double tmpCost = this->currSol.SolObjcost;
    const int lsSize = (level == my::lvlSmall) ? NB_NEIGHBORHOOD : 1;

  __repeat__:

    if (level == my::lvlSmall) {
        std::shuffle(this->permNeighborhood.begin(), this->permNeighborhood.end(), this->randEngine);
    }

    for (int i = 0; i < lsSize; ++i) {

        const int ls = (level == my::lvlSmall) ? this->permNeighborhood[i] : LS_EC_INSERT;

        this->localBest.SolObjcost = std::numeric_limits<float>::max();;

        const bool improv = this->neighborhood_operator(ls);
        if (improv) {
            goto __repeat__;
        }

    }

    return (this->currSol.SolObjcost < tmpCost - EPSILON);
}

bool Solver::neighborhood_operator(const int &index) {
    const double costOld = this->currSol.SolObjcost;
    this->taboo.reset();

    switch (index) {
        case LS_2OPT:
            this->ts_2opt_intra();
            break;

        case LS_INSERT:
            this->ts_realloc();
            break;

        case LS_SWAP:
            this->ts_swap();
            break;

        case LS_CROSSEX:
            this->ts_cross_exchange();
            break;

        case LS_2OPTSTAR:
            this->ts_2optstar();
            break;

        case LS_SINGLE_PATH_MOVE:
            this->ts_path_move(my::pathSize, my::pathSize);
            break;

        case LS_DOUBLE_PATH_MOVE:
            this->ts_path_move(my::pathDoubleSize, my::pathDoubleSize);
            break;

        case LS_EC_INSERT:
            this->ts_ec_realloc();
            break;

        default:
            break; // skip the process
    }

    this->update_best();
    this->current_to_best();
    return (this->currSol.SolObjcost < costOld - EPSILON);
}
