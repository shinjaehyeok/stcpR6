#ifndef BASELINE_INCREMENT_H
#define BASELINE_INCREMENT_H

#include "stcp_interface.h"

namespace stcp
{
    // Implementation of exponential baseline increment
    class ExpBaselineIncrement : public IBaselineIncrement
    {
    public:
        ExpBaselineIncrement()
            : m_lambda{0.0}
        {
        }
        ExpBaselineIncrement(const double &lambda)
            : m_lambda{lambda}
        {
        }
        virtual double computeLogBaseValue(const double &x) override = 0;
        // Exponential baseline can support a batch update
        // by using x_bar = 1/n * sum_{i=1}^n x_i and n values as inputs.
        // Note batch update should take n as a double rather than integer for generality.
        virtual double computeLogBaseValueByAvg(const double &x_bar, const double &n) override = 0;

    protected:
        double m_lambda{0};
    };

    // Derived classes from ExpBaselineIncrement
    // Normal baseline increment
    class Normal : public ExpBaselineIncrement
    {
    public:
        Normal()
            : Normal::Normal(0.0, 0.0, 1.0)
        {
        }
        Normal(const double &lambda)
            : Normal::Normal(lambda, 0.0, 1.0)
        {
        }
        Normal(const double &lambda, const double &mu, const double &sig)
            : ExpBaselineIncrement{lambda}
        {
            setupNormal(lambda, mu, sig);
        }
        double computeLogBaseValue(const double &x) override
        {
            return m_lambda * x - m_lambda_times_mu_plus_psi;
        }
        double computeLogBaseValueByAvg(const double &x_bar, const double &n) override
        {
            return n * Normal::computeLogBaseValue(x_bar);
        }
        

    protected:
        double m_mu{0.0};
        double m_sig{1.0};
        double m_psi{0.0};
        double m_lambda_times_mu_plus_psi{0.0};
        double compute_psi(const double &lambda, const double &sig)
        {
            return lambda * lambda * sig * sig * 0.5;
        }
        double compute_lambda_times_mu_plus_psi(const double &lambda, const double &mu, const double &sig)
        {
            return lambda * mu + compute_psi(lambda, sig);
        }
        void setupNormal(const double &lambda, const double &mu, const double &sig)
        {
            if (sig <= 0)
            {
                throw std::runtime_error("sig must be strictly positive.");
            }
            m_mu = mu;
            m_sig = sig;
            m_psi = compute_psi(lambda, sig);
            m_lambda_times_mu_plus_psi = compute_lambda_times_mu_plus_psi(lambda, mu, sig);
        }
    };

    // Bernoulli baseline increment
    class Ber : public ExpBaselineIncrement
    {
    public:
        Ber()
            : Ber::Ber(0.0, 0.5)
        {
        }
        Ber(const double &lambda)
            : Ber::Ber(lambda, 0.5)
        {
        }
        Ber(const double &lambda, const double &p)
            : ExpBaselineIncrement{lambda}
        {
            setupBer(lambda, p);
        }
        double computeLogBaseValue(const double &x) override
        {
            if (std::abs(x) < kEps)
            {
                return m_log_base_val_x_zero;
            }
            else if (std::abs(x - 1.0) < kEps)
            {
                return m_log_base_val_x_one;
            }
            else
            {
                throw std::runtime_error("Input must be either 0.0 or 1.0 or false or true.");
            }
        }
        double computeLogBaseValueByAvg(const double &x_bar, const double &n) override
        {
            return n * (m_lambda * x_bar + m_log_base_val_x_zero);
        }

    protected:
        double m_p{0.5};
        double m_log_base_val_x_one{0.0};
        double m_log_base_val_x_zero{0.0};
        double compute_psi_uncentered(const double &lambda, const double &p)
        {
            return log(1 - p + p * exp(lambda));
        }
        void check_prob_param_range(const double &p)
        {
            if (p <= 0.0 || p >= 1.0)
            {
                throw std::runtime_error("Probability parameter must be strictly inbetween 0.0 and 1.0.");
            }
        }
        void setupBer(const double &lambda, const double &p)
        {
            check_prob_param_range(p);
            m_p = p,
            m_log_base_val_x_one = lambda - compute_psi_uncentered(lambda, p);
            m_log_base_val_x_zero = -compute_psi_uncentered(lambda, p);
        }
    };

    class Bounded : public ExpBaselineIncrement
    {
    public:
        Bounded()
            : Bounded::Bounded(0.5, 0.5)
        {
        }
        Bounded(const double &lambda)
            : Bounded::Bounded(lambda, 0.5)
        {
        }
        Bounded(const double &lambda, const double &mu)
            : ExpBaselineIncrement{lambda}
        {
            setupBounded(lambda, mu);
        }
        double computeLogBaseValue(const double &x) override
        {
            if (x < 0.0)
            {
                throw std::runtime_error("Input must be non-negative.");
            }

            return log(1.0 + m_lambda * (x / m_mu - 1.0));
        }
        double computeLogBaseValueByAvg(const double &x_bar, const double &n) override
        {
            throw std::runtime_error("computeLogBaseValueByAvg cannot be used for the Bounded case.");
        }

    protected:
        double m_mu{0.5};
        void setupBounded(const double &lambda, const double &mu)
        {
            if (lambda >= 1.0 || lambda <= mu / (mu-1.0))
            {
                throw std::runtime_error("Lambda must be strictly inbetween mu / (mu-1.0) and 1.0.");
            }
            if (mu <= 0.0)
            {
                throw std::runtime_error("The mean parameter must be strictly positive.");
            }
            m_lambda = lambda;
            m_mu = mu;
        }
    };
} // End of namespace stcp
#endif