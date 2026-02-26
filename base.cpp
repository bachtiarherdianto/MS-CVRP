/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#include "base.h"

void my::error_function(const std::string &nameFunction, const std::string &comments) {         // TODO debugging
 std::cout << "ERROR in function " << nameFunction << ": " << comments << std::endl << std::flush;
 exit(-1);
}

std::string my::itos(const int & i) {
 std::stringstream s;
 s << i;
 return s.str();
}

double my::round_val(const double &r, const int &places) {
 const double off = std::pow(10.0f, places);
 const double x = r * off;
 double fractPart, intPart;  // round double to double with decimal place
 fractPart = std::modf(x, &intPart);

 if (std::abs(fractPart) >= 0.5f) {
  return (x >= 0u) ? std::ceil(x) / off : std::floor(x) / off;
 } else {
  return (x < 0u) ? std::ceil(x) / off : std::floor(x) / off;
 }

}

template <typename TT> void my::display_matrix(const std::vector<std::vector<TT>> &M) {  // TODO debugging function
 int i = 0;
 for (auto itline = M.begin(); itline != M.end(); ++itline) {
  std::cout << i << ":";
  for (auto itercolumn = (*itline).begin(); itercolumn != (*itline).end(); ++itercolumn) {
   std::cout << " " << *itercolumn;
  }

  std::cout << std::endl;
  ++i;
 }

}

template <typename TT> void my::display_vector(const std::vector<TT> &M) {  // TODO debugging function
 static int m_size = static_cast<int>(M.size());
 for (int i = 0; i < m_size; ++i) {
  std::cout << M[i] << " " << std::flush;
 }

 std::cout << std::endl;
}

template <typename TT> void my::display_list(const std::list<TT> &M) {  // TODO debugging function
 for (auto& it : M) {
  std::cout << it << " ";
 }

 std::cout << std::endl;
}

double my::KDComputeDistance(const std::array<double, 2>& a, const std::array<double, 2>& b) {
 return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

double my::KDComputeCoordinateDistance(double a, double b) {
 return (a - b) * (a - b);
}

uint64_t splitmix64(uint64_t &state) { // https://rosettacode.org/wiki/Pseudo-random_numbers/Splitmix64
 uint64_t result = (state += 0x9e3779b97f4a7c15); // update state
 result = (result ^ (result >> 30)) * 0xbf58476d1ce4e5b9;
 result = (result ^ (result >> 27)) * 0x94d049bb133111eb;
 return result ^ (result >> 31);
}

SplitMix64Engine::SplitMix64Engine(uint64_t seed) : state(seed) {}

SplitMix64Engine::result_type SplitMix64Engine::operator()() { // generate the next random number
 return splitmix64(state);
}

void display_help() {
 std::cout << "Usage: MS-CVRP <path-to-instance> [optional]\n\n";
 std::cout << "----------------------------------------------------------------------------------\n";
 std::cout << "Optional Arguments\n";
 std::cout << "[--seed <INT>]\t\t\tSets a fixed random seed.\n";
 std::cout << "[--max-time <INT>]\t\tSets the maximum computation time for the solver.\n";
 std::cout << "[--granular-size <INT>]\tSets the maximum number of granular neighborhoods.\n";
 std::cout << "[--min-truncated <INT>]\tSets the minimum Hamming distance for performing path relinking.\n";
 std::cout << "[--truncated-size <DOUBLE>]\tSets the truncation threshold for path relinking.\n";
 std::cout << "[--cache-capacity <INT>]\tSets the maximum number of customers stored for granular local search.\n";
 std::cout << "[--pair-cache-capacity <INT>]\tSets the maximum number of customer pairs stored for tabu search.\n";
 std::cout << "[--max-chain-iteration <INT>]\tSets the maximum number of iterations within the chained local search.\n";
 std::cout << "[--max-recombination <INT>]\tSets the maximum number of iterations for route recombination when generating the initial solution.\n";
 std::cout << "[--stop-threshold <INT>]\tSets the threshold for the approximate number of routes in a solution.\n";
 std::cout << "----------------------------------------------------------------------------------\n";
}

my::ArgType getArgType(const std::string& key) { // TODO moving to cpp is cleaner and more recommended
 static const std::unordered_map<std::string, my::ArgType> argMap = {
  {"--seed", my::ArgType::Seed},
  {"--max-time", my::ArgType::MaxTime},
  {"--granular-size", my::ArgType::GranularSize},
  {"--max-chain-iteration", my::ArgType::MaxChainIteration},
  {"--min-truncated", my::ArgType::MinTruncated},
  {"--stop-threshold", my::ArgType::StopThreshold},
  {"--truncated-size", my::ArgType::TruncatedSize},
  {"--cache-capacity", my::ArgType::CacheCapacity},
  {"--pair-cache-capacity", my::ArgType::PairCacheCapacity},
  {"--max-recombination", my::ArgType::MaxRecombination},
  {"--help", my::ArgType::Help}
 };

 auto it = argMap.find(key);
 if (it != argMap.end()) {
  return it->second;
 }

 return my::ArgType::Unknown;
}

Parameters::Parameters(const char * instance_path_) : fileInput(instance_path_) {}

const char* Parameters::get_file_input(){
 return this->fileInput;
}

double Parameters::get_truncated_threshold() const  {
 return this->truncatedThreshold;
}

double Parameters::get_min_hamming_size() const {
 return this->minTruncated;
}

int Parameters::get_running_time_max() const {
 return this->runningTimeMax;
}

int Parameters::get_init_recombination() const {
 return this->maxInitRecombination;
}

int Parameters::get_max_chain_iteration() const {
 return this->maxChainIteration;
}

int Parameters::get_cache_capacity() const {
 return this->cacheCapacity;
}

int Parameters::get_pair_cache_capacity() const {
 return this->cachePairCapacity;
}

int Parameters::get_granular_size(double routeStop) const {
 if (this->granularSize <= 50 &&
     routeStop >= this->routeStopThreshold) {
  return (2 * this->granularSize);
     }

 return this->granularSize;

}

int Parameters::get_elite_size(double routeStop) const {
 return (routeStop < this->routeStopThreshold) ? 3 : 2;
}

int Parameters::get_saving_matrix_size(double routeStop) const {
 return (routeStop < this->routeStopThreshold) ? 25 : 80;
}

int Parameters::get_initial_granular_size(double routeStop) const {
 return (routeStop < this->routeStopThreshold) ? 5 : 7;
}

double Parameters::get_threshold_growth(double routeStop) const {
 return (routeStop < this->routeStopThreshold) ? 100. : 150.;
}

double Parameters::get_regret_index(double routeStop) const {
 return (routeStop < this->routeStopThreshold) ? 0.6 : 0.3;
}

int Parameters::get_seed() const {
 return this->seed;
}

void Parameters::set(const std::string& key, const std::string& value) {
 switch (getArgType(key)) {
  case my::ArgType::Seed:
   this->seed = std::stoi(value);
   break;

  case my::ArgType::MaxTime:
   this->runningTimeMax = std::stoi(value);
   break;

  case my::ArgType::GranularSize:
   this->granularSize = std::stoi(value);
   break;

  case my::ArgType::MaxChainIteration:
   this->maxChainIteration = std::stoi(value);
   break;

  case my::ArgType::MinTruncated:
   this->minTruncated = std::stoi(value);
   break;

  case my::ArgType::StopThreshold:
   this->routeStopThreshold = std::stof(value);
   break;

  case my::ArgType::TruncatedSize:
   this->truncatedThreshold = std::stof(value);
   break;

  case my::ArgType::CacheCapacity:
   this->cacheCapacity = std::stoi(value);
   break;

  case my::ArgType::PairCacheCapacity:
   this->cachePairCapacity = std::stoi(value);
   break;

  case my::ArgType::MaxRecombination:
   this->maxInitRecombination = std::stoi(value);
   break;

  case my::ArgType::Help:
   display_help();
   exit(0);

  case my::ArgType::Unknown:
  default:
   std::cout << "UNKNOWN argument '" << key << "'.\n";
   display_help();
   exit(0);
 }

}

int CircleSector::positive_mod(const int &i)  {  // positive modulo 65536
 // 1) using the formula positive_mod(n,x) = (n % x + x) % x
 // 2) remark that "n % 65536" should be automatically compiled in an optimized form as "n & 0xffff" for faster calculations
 return (i % 65536 + 65536) % 65536;
}

void CircleSector::initialize(const int &point) {  // initialize a circle sector from a single point
 this->start = point;
 this->end = point;
}

bool CircleSector::isEnclosed(const int &point) {
 return (positive_mod(point - this->start) <= positive_mod(this->end - this->start)); // tests if a point is enclosed in the circle sector

}

bool CircleSector::overlap(const CircleSector &sector1, const CircleSector &sector2)  { // tests overlap of two circle sectors
 return ((positive_mod(sector2.start - sector1.start) <= positive_mod(sector1.end - sector1.start)) || (positive_mod(sector1.start - sector2.start) <= positive_mod(sector2.end - sector2.start)));
}

void CircleSector::extend(const int &point)  {  // extends the circle sector to include an additional point. done in a "greedy" way, such that the resulting circle sector is the smallest
 if (!isEnclosed(point)) {
  (positive_mod(point - this->end) <= positive_mod(this->start - point)) ? this->end = point : this->start = point;
 }

}