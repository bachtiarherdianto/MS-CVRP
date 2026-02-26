/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef BASE_H
#define BASE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <regex>
#include <iostream>
#include <list>
#include <unordered_map>
#include <string>
#include <vector>
#include <cstdint>
#include <queue>

#define NB_NEIGHBORHOOD (7)  // number of different local searches
#define LS_2OPT (0)
#define LS_INSERT (1)
#define LS_SWAP (2)
#define LS_CROSSEX (3)
#define LS_2OPTSTAR (4)
#define LS_SINGLE_PATH_MOVE (5)
#define LS_DOUBLE_PATH_MOVE (6)
#define LS_EC_INSERT (7)

#ifdef GUIDANCE
#define THRESHOLD_SETTING (4000)
#else
#define THRESHOLD_SETTING (3000)
#endif

#define CVRP (0)
#define REGRET_CONST (1.0)
#define REPAIR_CONST (2)
#define NEIGHBORHOOD_TABLE_SIZE(nb_cust) ((nb_cust < 1500) ? nb_cust : 1500)
#define INDEX_PERTURBATION(index) ((index > 0.5) ? 3 : 2)
#define EPSILON (0.00001)

namespace my {  // group of general function
  static constexpr int dummy = -1;
  static constexpr int depot = 0;
  static constexpr int indexGiantTour = 0;
  static constexpr int lvlSmall = 1;
  static constexpr int pathSize = 2;
  static constexpr int pathDoubleSize = 3;

  enum class ArgType {
    Seed,
    MaxTime,
    GranularSize,
    MaxChainIteration,
    MinTruncated,
    StopThreshold,
    TruncatedSize,
    CacheCapacity,
    PairCacheCapacity,
    MaxRecombination,
    Help,
    Unknown
  };

  void error_function(const std::string &nameFunction, const std::string &comments);          // TODO debugging

  std::string itos(const int & i);

  double round_val(const double &r, const int &places);

  template <typename TT> void display_matrix(const std::vector<std::vector<TT>> &M);  // TODO debugging function

  template <typename TT> void display_vector(const std::vector<TT> &M);

  template <typename TT> void display_list(const std::list<TT> &M);

  double KDComputeDistance(const std::array<double, 2>& a, const std::array<double, 2>& b);

  double KDComputeCoordinateDistance(double a, double b);
}

uint64_t splitmix64(uint64_t &state);

class SplitMix64Engine { // create a custom random engine that uses splitmix64
public:
  using result_type = uint64_t;

  SplitMix64Engine() = default;

  explicit SplitMix64Engine(uint64_t seed);

  ~SplitMix64Engine() = default;

  result_type operator()(); // generate the next random number

  static constexpr result_type min() { return 0; }

  static constexpr result_type max() { return ~static_cast<result_type>(0); }

private:
  uint64_t state;
};

void display_help();

my::ArgType getArgType(const std::string& key);

class Parameters {
public:
  Parameters(const char * instance_path_);

  const char* get_file_input() ;

  double get_truncated_threshold() const;

  double get_min_hamming_size() const ;

  int get_running_time_max() const;

  int get_init_recombination() const;

  int get_max_chain_iteration() const;

  int get_cache_capacity() const;

  int get_pair_cache_capacity() const;

  int get_granular_size(double routeStop) const;

  int get_elite_size(double routeStop) const;

  int get_saving_matrix_size(double routeStop) const;

  int get_initial_granular_size(double routeStop) const;

  double get_threshold_growth(double routeStop) const;

  double get_regret_index(double routeStop) const;

  int get_seed() const;

  void set(const std::string& key, const std::string& value);

private:
  int runningTimeMax = 30;
  int maxChainIteration = 24;
  int granularSize = 25;
  int cacheCapacity = 50;
  int cachePairCapacity = 200;
  int maxInitRecombination = 1000;
  double routeStopThreshold = 20.;
  double truncatedThreshold = 0.2;
  double minTruncated = 4;
  int seed = 0;

  const char * fileInput = "";
};

// Credit to Fabian Giesen at "https://web.archive.org/web/20200912191950/https://fgiesen.wordpress.com/2015/09/24/intervals-in-modular-arithmetic/" for useful implementation tips regarding interval overlaps in modular arithmetics
struct CircleSector {
  int start;
  int end;

  static int positive_mod(const int &i);

  void initialize(const int &point);

  bool isEnclosed(const int &point);

  static bool overlap(const CircleSector & sector1, const CircleSector & sector2);

  void extend(const int & point);

};

template <typename T> class Matrix2DFlat {
public:
  Matrix2DFlat() : rows(0u), cols(0u) {}  // declare the matrix

  void resize(size_t inputRows, size_t inputCols) {  // initially resize the matrix
    this->rows = inputRows;
    this->cols = inputCols;
    this->data.resize(this->rows * this->cols, 0u);
  }

  void at(const int i, const int j, const T value) { // insert value at position i & j
    this->data[i * this->cols + j] = value;
  }

  T at(const int i, const int j)  const  { return this->data[i * this->cols + j]; }  // retrieve (calling out) value at pos i & j

  int get_cols() { return static_cast<int>(this->cols); }

private:
  std::vector<T> data;
  size_t rows;
  size_t cols;
};

class KDTree { // a simple implementation of a kd-tree based on https://github.com/cdalitz/kdtree-cpp
public:
  KDTree(const Matrix2DFlat<double>& coords, const int num_nodes) {
    std::array<double, 2> lowBound = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 2> highBound = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

    for (int i = 0; i < static_cast<int>(num_nodes); ++i) {
      lowBound[0] = std::min(lowBound[0], coords.at(i,0));
      lowBound[1] = std::min(lowBound[1], coords.at(i,1));
      highBound[0] = std::max(highBound[0], coords.at(i,0));
      highBound[1] = std::max(highBound[1], coords.at(i,1));
      this->nodes.emplace_back(i, coords.at(i,0), coords.at(i,1));
    }

    this->root = BuildTree(0, 0, this->nodes.size(), lowBound, highBound);
  }

  ~KDTree() {
    delete root;
  }

  void FIllNeighborhoodTable(double x, double y, int k, int node, Matrix2DFlat<int> &neighbors) const {
    KDTreeHeap heap;
    this->FindNeighbors(this->root, heap, {x, y}, k);

    while (!heap.empty()) {
      const HeapNode& heap_node = heap.top();
      --k;
      neighbors.at(node, k, this->nodes[heap_node.point_index].index);
      heap.pop();
    }

  }

private:
  struct Point {
    Point(int index, double x, double y) : index(index), coords({x, y}) { }
    int index;
    std::array<double, 2> coords; // x and y coordinates
  };

  std::vector<Point> nodes;

  struct Node {
    Node() = default;
    ~Node() {
      delete left;
      delete right;
    }

    int cutdim;
    Node *left, *right;
    std::array<double, 2> lowBound, highBound;
    int point_index;
  };

  Node* root = nullptr;

  struct HeapNode {
    HeapNode(int point_index, double distance) : point_index(point_index), distance(distance) { }
    int point_index;
    double distance;

  };

  struct HeapNodeComparator {
    bool operator()(const HeapNode& a, const HeapNode& b) const {
      return a.distance < b.distance;
    }

  };

  using KDTreeHeap = std::priority_queue<HeapNode, std::vector<HeapNode>, HeapNodeComparator>;

  Node* BuildTree(int depth, int begin, int end, const std::array<double, 2>& lowBound, const std::array<double, 2>& highBound)  {

    const int dimension = depth % 2;

    KDTree::Node* node = new KDTree::Node();
    node->cutdim = dimension;
    node->left = nullptr;
    node->right = nullptr;
    node->lowBound = lowBound;
    node->highBound = highBound;

    if (end - begin <= 1) {
      node->point_index = begin;
    } else {
      int median = (begin + end) / 2;
      std::nth_element(
              this->nodes.begin() + begin,
              this->nodes.begin() + median,
              this->nodes.begin() + end,
              [dimension](const Point& a, const Point& b) {
                return a.coords[dimension] < b.coords[dimension];
              }
      );

      node->point_index = median;

      const int cutVal = this->nodes[median].coords[dimension];

      if (median - begin > 0) {
        std::array<double, 2> nextHighBound = highBound;
        nextHighBound[dimension] = cutVal;
        node->left = BuildTree(depth + 1, begin, median, node->lowBound, nextHighBound);

      }

      if (end - median > 1) {
        std::array<double, 2> nextLowBound = lowBound;
        nextLowBound[dimension] = cutVal;
        node->right = BuildTree(depth + 1, median + 1, end, nextLowBound, highBound);
      }

    }

    return node;
  }

  bool FindNeighbors(KDTree::Node* node, KDTree::KDTreeHeap& heap, const std::array<double, 2>& point, int k) const {
    double currDist = my::KDComputeDistance(point, nodes[node->point_index].coords);
    if (static_cast<int>(heap.size()) < k) {
      heap.push(HeapNode(node->point_index, currDist));
    } else if (currDist < heap.top().distance) {
      heap.pop();
      heap.push(HeapNode(node->point_index, currDist));
    }

    if (point[node->cutdim] < this->nodes[node->point_index].coords[node->cutdim]) {
      if (node->left) {
        if (FindNeighbors(node->left, heap, point, k)) {
          return true;
        }

      }

    } else {
      if (node->right) {
        if (FindNeighbors(node->right, heap, point, k)) {
          return true;
        }

      }

    }

    double dist = static_cast<int>(heap.size()) < k ? std::numeric_limits<double>::max() : heap.top().distance;

    if (point[node->cutdim] < this->nodes[node->point_index].coords[node->cutdim]) {
      if (node->right && BoundsOverlapBall(point, dist, node->right)) {
        if (FindNeighbors(node->right, heap, point, k)) {
          return true;
        }

      }

    } else {
      if (node->left && BoundsOverlapBall(point, dist, node->left)) {
        if (FindNeighbors(node->left, heap, point, k)) {
          return true;
        }

      }

    }

    if (static_cast<int>(heap.size()) == k) {
      dist = heap.top().distance;
    }

    return BallWithinBounds(point, dist, node);
  }

  bool BoundsOverlapBall(const std::array<double, 2>& point, double dist, KDTree::Node* node) const {
    double distsum = 0.;
    for (int i = 0; i < static_cast<int>(point.size()); ++i) {
      if (point[i] < node->lowBound[i]) {
        distsum += my::KDComputeCoordinateDistance(point[i], node->lowBound[i]);
        if (distsum > dist) {
          return false;
        }

      } else if (point[i] > node->highBound[i]) {
        distsum += my::KDComputeCoordinateDistance(point[i], node->highBound[i]);
        if (distsum > dist) {
          return false;
        }

      }

    }

    return true;
  }

  bool BallWithinBounds(const std::array<double, 2>& point, double dist, KDTree::Node* node) const {

    for (int i = 0; i < static_cast<int>(point.size()); ++i) {
      if (my::KDComputeCoordinateDistance(point[i], node->lowBound[i]) <= dist ||
          my::KDComputeCoordinateDistance(point[i], node->highBound[i]) <= dist) {
        return false;
      }

    }

    return true;
  }
};

#endif //BASE_H