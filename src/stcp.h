#ifndef STCP_H
#define STCP_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>

namespace stcp
{
constexpr double kEps { 1e-8 };

// Implementation of mixture of SR-type E-detector
class MixSR
{
public:
  MixSR() = default;
  MixSR(const std::vector<double>& weights)
    : m_weights { weights }
    {
      validateWeights(weights);
      m_log_values.assign(weights.size(), 0.0);
    }
  MixSR(const std::vector<double>& weights, const std::vector<double>& log_values)
    : m_weights { weights }, m_log_values {log_values}
    {
      if (weights.size() != log_values.size())
      {
        throw std::runtime_error("Weights and values do not have the same length.");
      }
      
      validateWeights(weights);
    }
  
  void print()
  {
    std::cout << "weights: " << std::endl;
    for (auto w : m_weights)
    {
      std::cout << w <<  ' ';
    }
    std::cout << '\n';
    std::cout << "log of values: " << std::endl;
    for (auto v : m_log_values)
    {
      std::cout << v <<  ' ';
    }
    std::cout << '\n';
  }
  
  double getLogMixedValue()
  {
    
  }
  
private:
  std::vector<double> m_weights {1};
  std::vector<double> m_log_values {0};
  
  void validateWeights (const std::vector<double>& weights) 
  {
    double weights_sum {0};
    for (auto& w : weights)
    {
      if (w <= 0)
      {
        throw std::runtime_error("All weights must be strictly positive.");
      }
      weights_sum += w;
    }     
    if (std::abs(weights_sum - 1.0) > kEps)
    {
      throw std::runtime_error("Sum of weights is not equal to 1.");
    }
  }
};
}
#endif