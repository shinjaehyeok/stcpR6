#ifndef LOG_LR_INCREMENT_H
#define LOG_LR_INCREMENT_H

#include "stcp_interface.h"

namespace stcp
{
    // Implementation of log likelihood ratio baseline increment
    // Bernoulli
    class BerLR : public ILogLRIncrement
    {
    public:
        BerLR() {}
        BerLR(const double q)
            : BerLR::BerLR(q, 0.5)
        {
        }
        BerLR(const double q, const double p)
            : BerLR::BerLR()
        {
            setupBerLR(q, p);
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
        void updateH1MLE(double &h_1_mle, const double x, const int n) override
        {
            // n: sample size *after" taking account for new observation x
            h_1_mle = (h_1_mle * (n - 1) + x) / n;
        }
        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            if (abs(h_1_mle) < kEps)
            {
                return n * (1 - h_1_mle) * log((1 - h_1_mle) / (1 - this->m_p));
            }
            else if (abs(h_1_mle - 1.0) < kEps)
            {
                return n * h_1_mle * log(h_1_mle / this->m_p);
            }
            return n * (h_1_mle * log(h_1_mle / this->m_p) +
                        (1 - h_1_mle) * log((1 - h_1_mle) / (1 - this->m_p)));
        }

    protected:
        double m_q{0.5};
        double m_p{0.5};
        double m_log_base_val_x_one{0.0};
        double m_log_base_val_x_zero{0.0};
        void check_prob_param_range(const double p)
        {
            if (p <= 0.0 || p >= 1.0)
            {
                throw std::runtime_error("Probability parameter must be strictly inbetween 0.0 and 1.0.");
            }
        }
        void setupBerLR(const double q, const double p)
        {
            check_prob_param_range(q);
            check_prob_param_range(p);
            m_q = q;
            m_p = p;
            m_log_base_val_x_one = log(q / p);
            m_log_base_val_x_zero = log((1 - q) / (1 - p));
        }
    };
    class BerGLR : public BerLR
    {
    public:
        BerGLR()
            : BerLR::BerLR(0.5, 0.5)
        {
        }
        BerGLR(const double p)
            : BerLR::BerLR(p, p)
        {
        }
    };
} // End of namespace stcp
#endif