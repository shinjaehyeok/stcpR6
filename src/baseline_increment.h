#ifndef BASELINE_INCREMENT_H
#define BASELINE_INCREMENT_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

namespace stcp
{
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
} // End of namespace stcp
#endif