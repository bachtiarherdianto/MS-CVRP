/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef WELFORD_H
#define WELFORD_H

// Welford's online algorithm to compute mean and standard deviation
class Welford {
public:
    Welford() = default;

    void update(const double &x) {
        ++this->count;
        this->sum += x;
        this->sumQ += (x * x);
        const double newMean = this->sum / static_cast<double>(this->count);
        const double newVar = (this->sumQ / static_cast<double>(this->count)) - (newMean * newMean);
        this->mean = newMean;
        this->var = newVar;
    }

    [[nodiscard]] inline double get_mean() const {
        return static_cast<double>(this->mean);
    }

    [[nodiscard]] inline double get_std_dev() const {
        return std::sqrt(this->var);
    }

private:
    unsigned long count = 0ul;
    double sum = 0.0f;
    double sumQ = 0.0f;
    double mean = 0.0f;
    double var = 0.0f;

};

#endif //WELFORD_H