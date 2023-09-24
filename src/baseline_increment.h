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
        ExpBaselineIncrement(const double lambda)
            : m_lambda{lambda}
        {
        }
        virtual double computeLogBaseValue(const double x) override = 0;

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
        Normal(const double lambda)
            : Normal::Normal(lambda, 0.0, 1.0)
        {
        }
        Normal(const double lambda, const double mu, const double sig)
            : ExpBaselineIncrement{lambda},
              m_mu{mu},
              m_sig{sig},
              m_psi{compute_psi(lambda, sig)},
              m_lambda_times_mu_plus_psi{compute_lambda_times_mu_plus_psi(lambda, mu, sig)}
        {
        }
        double computeLogBaseValue(const double x) override
        {
            return m_lambda * x - m_lambda_times_mu_plus_psi;
        }

    private:
        double m_mu{0.0};
        double m_sig{1.0};
        double m_psi{0.0};
        double m_lambda_times_mu_plus_psi{0.0};
        double compute_psi(const double lambda, const double sig)
        {
            return lambda * lambda * sig * sig * 0.5;
        }
        double compute_lambda_times_mu_plus_psi(const double lambda, const double mu, const double sig)
        {
            return lambda * mu + compute_psi(lambda, sig);
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
        Ber(const double lambda)
            : Ber::Ber(lambda, 0.5)
        {
        }
        Ber(const double lambda, const double p)
            : ExpBaselineIncrement{lambda},
              m_p{p},
              m_log_base_val_x_one{lambda - compute_psi_uncentered(lambda, p)},
              m_log_base_val_x_zero{-compute_psi_uncentered(lambda, p)}
        {
        }
        double computeLogBaseValue(const double x) override
        {
            if (abs(x) < kEps)
            {
                return m_log_base_val_x_zero;
            }
            else if (abs(x - 1.0) < kEps)
            {
                return m_log_base_val_x_one;
            }
            else
            {
                throw std::runtime_error("Input must be either 0.0 or 1.0 or false or true.");
            }
        }

    private:
        double m_p{0.5};
        double m_log_base_val_x_one{0.0};
        double m_log_base_val_x_zero{0.0};
        double compute_psi_uncentered(const double lambda, const double p)
        {
            return log(1 - p + p * exp(lambda));
        }
    };
} // End of namespace stcp
#endif