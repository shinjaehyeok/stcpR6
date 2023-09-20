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

// Implementation of Baseline increment
class BaselineIncrement
{
public:
  BaselineIncrement()
    : m_lambda{0}, m_s_param{0}, m_v_param{1}
    {
    }
  BaselineIncrement(const double lambda, const double s_param, const double v_param)
    : m_lambda{lambda}, m_s_param{s_param}, m_v_param{v_param}
    {
    }
  
  virtual double computeLogBaseValue(const double x) = 0;
  virtual void reinitialize(const double lambda, const double s_param, const double v_param) = 0;
  
protected:
  double m_lambda;
  double m_s_param;
  double m_v_param;
};

// Implementation of E-value / detector
class ISingleE
{
public:
  virtual double getLogValue() = 0;
  virtual void reinitialize(const double log_value, const double lambda, const double s_param, const double v_param) = 0;
  virtual void updateLogValue(const double x) = 0;
  
  virtual ~ISingleE() {}
};

template <typename L>
class SingleE : public ISingleE
{
  static_assert(std::is_base_of<BaselineIncrement, L>::value, "Type must be derived from BaselineIncrement class.");
  
public:
  SingleE()
    : m_log_value{0.0}, m_base_obj{}
    {
    }
  SingleE(const double log_value, const double lambda, const double s_param, const double v_param)
    : m_log_value{log_value}, m_base_obj{lambda, s_param, v_param}
    {
    }
  
  double getLogValue() override { return m_log_value; }
  void reinitialize(const double log_value, const double lambda, const double s_param, const double v_param) override
  {
    m_log_value = log_value;
    m_base_obj.reinitialize(lambda, s_param, v_param);
  }
  
  virtual void updateLogValue(const double x) = 0;
  
protected:
  double m_log_value;
  L m_base_obj;
};

// Implementation of mixture of E-values / detectors
template <typename E>
class MixE
{
  static_assert(std::is_base_of<ISingleE, E>::value, "Type must be derived from SingleE class.");
  
public:
  MixE();
  MixE(const double log_value, const double lambda, const double s_param, const double v_param);
  MixE(const std::vector<double> &weights,
       const std::vector<double> &log_values,
       const std::vector<double> &lambdas,
       const double s_param,
       const double v_param);
  
  std::vector<double> getWeights() { return m_weights; }
  std::vector<double> getLogWeights() { return m_log_weights; }
  std::vector<double> getLogValues();
  
  void print();
  double getLogMixedValue();
  
private:
  std::vector<double> m_weights;
  std::vector<double> m_log_weights;
  std::vector<E> m_e_objs;
  
  std::vector<double> validateAndComputeLogWeights(const std::vector<double> &weights);
};

// Constructors
template <typename E>
inline MixE<E>::MixE()
  : m_weights{1.0}, m_log_weights{0.0}, m_e_objs(1)
  {
  }
template <typename E>
inline MixE<E>::MixE(const std::vector<double> &weights,
                     const std::vector<double> &log_values,
                     const std::vector<double> &lambdas,
                     const double s_param,
                     const double v_param)
  : m_weights{weights},
    m_log_weights{validateAndComputeLogWeights(weights)},
    m_e_objs(log_values.size())
    {
      if (weights.size() != log_values.size())
      {
        throw std::runtime_error("Weights and values do not have the same length.");
      }
      for (std::size_t i = 0; i < log_values.size(); i++)
      {
        m_e_objs[i].reinitialize(log_values[i], lambdas[i], s_param, v_param);
      }
    }
template <typename E>
inline MixE<E>::MixE(const double log_value, const double lambda, const double s_param, const double v_param)
  : MixE<E>::MixE(std::vector<double>{1.0},
    std::vector<double>{log_value},
    std::vector<double>{lambda},
    s_param, v_param) {}

// Public members
template <typename E>
inline void MixE<E>::print()
{
  std::cout << "weights: " << std::endl;
  for (auto w : m_weights)
  {
    std::cout << w << ' ';
  }
  std::cout << '\n';
  std::cout << "log of values: " << std::endl;
  for (auto v : m_e_objs)
  {
    std::cout << v.getLogValue() << ' ';
  }
  std::cout << '\n';
}

template <typename E>
inline std::vector<double> MixE<E>::getLogValues()
{
  std::vector<double> log_values(m_e_objs.size());
  for (std::size_t i = 0; i < m_e_objs.size(); i++)
  {
    log_values[i] = m_e_objs[i].getLogValue();
  }
  return log_values;
}

template <typename E>
inline double MixE<E>::getLogMixedValue()
{
  if (m_e_objs.size() == 1)
  {
    // Since weight must be equal to 1 by the construction,
    // we do not need to take account of it.
    return m_e_objs[0].getLogValue();
  }
  auto log_wls{m_log_weights};
  for (std::size_t i = 0; i < log_wls.size(); i++)
  {
    log_wls[i] += m_e_objs[i].getLogValue();
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
template <typename E>
inline std::vector<double> MixE<E>::validateAndComputeLogWeights(const std::vector<double> &weights)
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

// Derived classes from BaselineIncrement

// Normal baseline increment
// s_param : mu0
// v_param : sigma
class Normal : public BaselineIncrement
{
public:
  Normal()
    : BaselineIncrement(), m_psi{0}, m_lambda_times_mu_plus_psi{0}
    {
    }
  Normal(const double lambda, const double s_param, const double v_param)
    : BaselineIncrement(lambda, s_param, v_param),
      m_psi{compute_psi(lambda, v_param)},
      m_lambda_times_mu_plus_psi{compute_lambda_times_mu_plus_psi(lambda, s_param, v_param)}
      {
      }
  double computeLogBaseValue(const double x) override
  {
    return m_lambda * x - m_lambda_times_mu_plus_psi;
  }
  void reinitialize(const double lambda, const double s_param, const double v_param) override
  {
    m_lambda = lambda;
    m_s_param = s_param;
    m_v_param = v_param;
    m_psi = compute_psi(lambda, v_param);
    m_lambda_times_mu_plus_psi = compute_lambda_times_mu_plus_psi(lambda, s_param, v_param);
  }
  
private:
  double m_psi;
  double m_lambda_times_mu_plus_psi;
  double compute_psi(const double lambda, const double v_param)
  {
    return lambda * lambda * v_param * v_param * 0.5;
  }
  double compute_lambda_times_mu_plus_psi(const double lambda, const double s_param, const double v_param)
  {
    return lambda * s_param + compute_psi(lambda, v_param);
  }
};

// Derived classes from SingleE
template <typename L>
class ST : public SingleE<L>
{
public:
  using SingleE<L>::SingleE;
  void updateLogValue(const double x) override
  {
    this->m_log_value += this->m_base_obj.computeLogBaseValue(x);
  }
};

} // End of namespace stcp
#endif