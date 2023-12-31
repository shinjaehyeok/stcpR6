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

#ifndef STCP_H
#define STCP_H

#include "stcp_interface.h"
#include "baseline_increment.h"
#include "log_lr_increment.h"
#include "log_lr_e.h"
#include "baseline_e.h"
#include "mix_e.h"

namespace stcp
{
    // Implementation of mixture of E-values / detectors
    template <typename E>
    class Stcp : public IStcp
    {
        static_assert(
            std::is_base_of<IGeneralE, E>::value,
            "Type must be derived from IGeneralE class.");

    public:
        Stcp() {}
        Stcp(const E &e_obj)
            : m_e_obj{e_obj} {}
        Stcp(const E &e_obj, const double &threshold)
            : m_e_obj{e_obj}, m_threshold{threshold} {}

        double getLogValue() override { return m_e_obj.getLogValue(); }
        double getThreshold() override { return m_threshold; };

        bool isStopped() override { return m_is_stopped; };
        int getTime() override { return m_time; };
        int getStoppedTime() override { return m_stopped_time; }

        void reset() override
        {
            m_e_obj.reset();
            m_is_stopped = false;
            m_time = 0;
            m_stopped_time = 0;
        }

        void updateLogValue(const double &x) override;
        void updateLogValues(const std::vector<double> &xs) override;
        void updateLogValuesUntilStop(const std::vector<double> &xs) override;

        double updateAndReturnHistory(const double &x) override;
        std::vector<double> updateAndReturnHistories(const std::vector<double> &xs) override;

    protected:
        E m_e_obj{};
        double m_threshold{log(1.0 / 0.05)}; // Default threshold ues alpha = 0.05.
        int m_time{0};
        bool m_is_stopped{false};
        int m_stopped_time{0};
    };

    // Public members
    template <typename E>
    inline void Stcp<E>::updateLogValue(const double &x)
    {
        m_e_obj.updateLogValue(x);
        m_time++;
        if (this->getLogValue() > m_threshold)
        {
            if (!m_is_stopped)
            {
                // Record the first stopped time only.
                m_stopped_time = m_time;
                m_is_stopped = true;
            }
        }
    }
    template <typename E>
    inline void Stcp<E>::updateLogValues(const std::vector<double> &xs)
    {
        for (auto x : xs)
        {
            this->updateLogValue(x);
        }
    }
    template <typename E>
    inline void Stcp<E>::updateLogValuesUntilStop(const std::vector<double> &xs)
    {
        for (auto x : xs)
        {
            this->updateLogValue(x);
            if (m_is_stopped)
            {
                break;
            }
        }
    }
    template <typename E>
    inline double Stcp<E>::updateAndReturnHistory(const double &x)
    {
        this->updateLogValue(x);
        return this->getLogValue();
    }
    template <typename E>
    inline std::vector<double> Stcp<E>::updateAndReturnHistories(const std::vector<double> &xs)
    {
        std::vector<double> log_values(xs.size());
        for (std::size_t i = 0; i < xs.size(); i++)
        {
            log_values[i] = this->updateAndReturnHistory(xs[i]);
        }
        return log_values;
    }
} // End of namespace stcp
#endif