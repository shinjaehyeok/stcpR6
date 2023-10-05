/*
 * Copyright (c) 2023 Jaehyeok Shin
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LOG_LR_INCREMENT_H
#define LOG_LR_INCREMENT_H

#include "stcp_interface.h"

namespace stcp
{
    // Implementation of log likelihood ratio baseline increment
    // Normal
    class NormalLR : public ILogLRIncrement
    {
    public:
        NormalLR() {}
        NormalLR(const double mu1)
            : NormalLR::NormalLR(mu1, 0.0, 1.0)
        {
        }
        NormalLR(const double mu1, const double mu, const double sig)
            : NormalLR::NormalLR()
        {
            setupNormalLR(mu1, mu, sig);
        }
        double computeLogBaseValue(const double x) override
        {
            return m_mu_delta_by_sig_squared * (x - m_mu1_plus_mu_by_two);
        }
        void updateH1MLE(double &h_1_mle, const double x, const int n) override
        {
            // n: sample size *after" taking account for new observation x
            h_1_mle = (h_1_mle * (n - 1) + x) / n;
        }
        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            return n * ((h_1_mle - m_mu) / m_sig) * ((h_1_mle - m_mu) / m_sig) / 2.0;
        }

    protected:
        double m_mu{0.0};
        double m_sig{1.0};
        double m_mu_delta_by_sig_squared{0.0};
        double m_mu1_plus_mu_by_two{0.0};
        void setupNormalLR(const double mu1, const double mu, const double sig)
        {
            if (sig <= 0)
            {
                throw std::runtime_error("sig must be strictly positive.");
            }
            m_mu = mu;
            m_sig = sig;
            m_mu_delta_by_sig_squared = (mu1 - mu) / sig / sig;
            m_mu1_plus_mu_by_two = (mu1 + mu) / 2.0;
        }
    };
    class NormalGLR : public NormalLR
    {
    public:
        NormalGLR()
            : NormalLR::NormalLR(0.0, 0.0, 1.0)
        {
        }
        NormalGLR(const double mu)
            : NormalLR::NormalLR(mu, mu, 1.0)
        {
        }
        NormalGLR(const double mu, const double sig)
            : NormalLR::NormalLR(mu, mu, sig)
        {
        }
    };

    class NormalGLRGreater : public NormalGLR
    {
    public:
        using NormalGLR::NormalGLR;

        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            double h_1_mle_gt{std::max(h_1_mle, m_mu)};
            return n * ((h_1_mle_gt - m_mu) / m_sig) * ((h_1_mle_gt - m_mu) / m_sig) / 2.0;
        }
    };

    class NormalGLRLess : public NormalGLR
    {
    public:
        using NormalGLR::NormalGLR;

        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            double h_1_mle_ls{std::min(h_1_mle, m_mu)};
            return n * ((h_1_mle_ls - m_mu) / m_sig) * ((h_1_mle_ls - m_mu) / m_sig) / 2.0;
        }
    };

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
            return computeMaxLLRBer(m_p, h_1_mle, n);
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
        double computeMaxLLRBer(const double p, const double x_bar, const int n);
    };

    inline double BerLR::computeMaxLLRBer(const double p, const double x_bar, const int n)
    {
        if (abs(x_bar) < kEps)
        {
            return n * (1 - x_bar) * log((1 - x_bar) / (1 - p));
        }
        else if (abs(x_bar - 1.0) < kEps)
        {
            return n * x_bar * log(x_bar / p);
        }
        return n * (x_bar * log(x_bar / p) +
                    (1 - x_bar) * log((1 - x_bar) / (1 - p)));
    }

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

    class BerGLRGreater : public BerGLR
    {
    public:
        using BerGLR::BerGLR;

        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            return computeMaxLLRBer(this->m_p, std::max(h_1_mle, this->m_p), n);
        }
    };

    class BerGLRLess : public BerGLR
    {
    public:
        using BerGLR::BerGLR;

        double computeMaxLLR(const double h_1_mle, const int n) override
        {
            return computeMaxLLRBer(this->m_p, std::min(h_1_mle, this->m_p), n);
        }
    };

} // End of namespace stcp
#endif