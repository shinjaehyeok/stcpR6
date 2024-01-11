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

#ifndef MIX_E_H
#define MIX_E_H

#include "stcp_interface.h"

namespace stcp
{
    // Implementation of mixture of E-values / detectors
    template <typename E>
    class MixE : public IGeneralE
    {
        static_assert(std::is_base_of<IGeneralE, E>::value, "Type must be derived from IGeneralE class.");

    public:
        MixE();
        MixE(const E &e_obj);
        MixE(const std::vector<E> &e_objs,
             const std::vector<double> &weights);

        double getLogValue() override;
        void reset() override;
        void updateLogValue(const double &x) override;
        void updateLogValueByAvg(const double &x_bar, const double &n) override;

        std::vector<double> getWeights() { return m_weights; }
        std::vector<double> getLogValues();

        void print();

    protected:
        std::vector<E> m_e_objs;
        std::vector<double> m_weights;
        std::vector<double> m_log_weights;
        std::vector<double> validateAndComputeLogWeights(const std::vector<double> &weights);
    };

    // Constructors
    template <typename E>
    inline MixE<E>::MixE()
        : MixE<E>::MixE(std::vector<E>(1),
                        std::vector<double>{1.0})
    {
    }
    template <typename E>
    inline MixE<E>::MixE(const E &e_obj)
        : MixE<E>::MixE(std::vector<E>{e_obj},
                        std::vector<double>{1.0})
    {
    }
    template <typename E>
    inline MixE<E>::MixE(const std::vector<E> &e_objs,
                         const std::vector<double> &weights)
        : m_e_objs{e_objs},
          m_weights{weights},
          m_log_weights{validateAndComputeLogWeights(weights)}

    {
        if (e_objs.size() != weights.size())
        {
            throw std::runtime_error("E objects and Weights do not have the same length.");
        }
    }

    // Public members
    template <typename E>
    inline double MixE<E>::getLogValue()
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

        return logSumExp(log_wls);
    }

    template <typename E>
    inline void MixE<E>::reset()
    {
        for (auto &e_obj : m_e_objs)
        {
            e_obj.reset();
        }
    }

    template <typename E>
    inline void MixE<E>::updateLogValue(const double &x)
    {
        for (auto &e_obj : m_e_objs)
        {
            e_obj.updateLogValue(x);
        }
    }

    template <typename E>
    inline void MixE<E>::updateLogValueByAvg(const double &x_bar, const double &n)
    {
        for (auto &e_obj : m_e_objs)
        {
            e_obj.updateLogValueByAvg(x_bar, n);
        }
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
    inline void MixE<E>::print()
    {
        std::cout << "weights: " << std::endl;
        for (auto &w : m_weights)
        {
            std::cout << w << ' ';
        }
        std::cout << '\n';
        std::cout << "log of values: " << std::endl;
        for (auto &e_obj : m_e_objs)
        {
            std::cout << e_obj.getLogValue() << ' ';
        }
        std::cout << '\n';
    }
    // Private members
    template <typename E>
    inline std::vector<double> MixE<E>::validateAndComputeLogWeights(const std::vector<double> &weights)
    {
        double weights_sum{0.0};
        std::vector<double> log_weights;
        log_weights.reserve(weights.size());
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