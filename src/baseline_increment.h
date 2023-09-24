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
} // End of namespace stcp
#endif