#ifndef STCP_H
#define STCP_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

namespace stcp
{
constexpr double kEps{1e-8};

// Implementation of mixture of E-detectors
class MixE
{
public:
  MixE();
  MixE(const double weight);
  MixE(const double weight, const double log_value);
  MixE(const std::vector<double> &weights);
  MixE(const std::vector<double> &weights, const std::vector<double> &log_values);
  
  std::vector<double> getWeights() { return m_weights; }
  std::vector<double> getLogValues() { return m_log_values; }
  std::vector<double> getLogWeights() { return m_log_weights; }
  
  void print();
  double getLogMixedValue();
  
private:
  std::vector<double> m_weights;
  std::vector<double> m_log_values;
  std::vector<double> m_log_weights;
  
  std::vector<double> validateAndComputeLogWeights(const std::vector<double> &weights);
};

// Constructors
inline MixE::MixE()
  : m_weights{1.0}, m_log_values{0.0}, m_log_weights{0.0}
  {
  }
inline MixE::MixE(const std::vector<double> &weights)
  : m_weights{weights}
  {
    m_log_weights = validateAndComputeLogWeights(weights);
    m_log_values.assign(weights.size(), 0.0);
  }
inline MixE::MixE(const std::vector<double> &weights, const std::vector<double> &log_values)
  : m_weights{weights}, m_log_values{log_values}
  {
    if (weights.size() != log_values.size())
    {
      throw std::runtime_error("Weights and values do not have the same length.");
    }
    m_log_weights = validateAndComputeLogWeights(weights);
  }
inline MixE::MixE(const double weight) : MixE::MixE(std::vector<double>{weight}) {}
inline MixE::MixE(const double weight, const double log_value) : MixE::MixE(std::vector<double>{weight}, std::vector<double>{log_value}) {}

// Public members
inline void MixE::print()
{
  std::cout << "weights: " << std::endl;
  for (auto w : m_weights)
  {
    std::cout << w << ' ';
  }
  std::cout << '\n';
  std::cout << "log of values: " << std::endl;
  for (auto v : m_log_values)
  {
    std::cout << v << ' ';
  }
  std::cout << '\n';
}

inline double MixE::getLogMixedValue()
{
  if (m_log_values.size() == 1)
  {
    // Since weight must be equal to 1 by the construction,
    // we do not need to take account of it.
    return m_log_values[0];
  }
  auto log_wls{m_log_weights};
  for (std::size_t i = 0; i < log_wls.size(); i++)
  {
    log_wls[i] += m_log_values[i];
  }
  
  double max_log_wl = *std::max_element(log_wls.begin(), log_wls.end());
  double sum_exp{0.0};
  for (auto &log_wl : log_wls)
  {
    sum_exp += std::exp(log_wl - max_log_wl);
  }
  
  return log(sum_exp) + max_log_wl;
}

// Private members
inline std::vector<double> MixE::validateAndComputeLogWeights(const std::vector<double> &weights)
{
  double weights_sum{0.0};
  std::vector<double> log_weights;
  
  for (auto &w : weights)
  {
    if (w <= 0)
    {
      throw std::runtime_error("All weights must be strictly positive.");
    }
    weights_sum += w;
    log_weights.push_back(std::log(w));
  }
  if (std::abs(weights_sum - 1.0) > kEps)
  {
    throw std::runtime_error("Sum of weights is not equal to 1.");
  }
  
  return log_weights;
}

} // End of namespace stcp
#endif