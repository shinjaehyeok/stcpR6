/*
 Copyright (c) 2023 Jaehyeok Shin

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
    // vector input
    // stoppiong rule
    // combine multiple Es

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
        Stcp(const E &e_obj, const double threshold)
            : m_e_obj{e_obj}, m_threshold{threshold} {}

        double getLogValue() override { return m_e_obj.getLogValue(); }
        double getThreshold() override { return m_threshold; };

        bool isStopped() override { return m_is_stopped; };
        int getTime() override { return m_time; };
        int getStoppedTime() override { return m_stopped_time; }

        void reset() override
        {
            m_e_obj.resetLogValue();
            m_is_stopped = false;
            m_time = 0;
            m_stopped_time = kIntInf;
        }

        void updateLogValue(const double x) override;
        void updateLogValues(const std::vector<double> &xs) override;
        void updateLogValuesUntilStop(const std::vector<double> &xs) override;

        double updateAndReturnHistory(const double x) override;
        std::vector<double> updateAndReturnHistories(const std::vector<double> &xs) override;

    protected:
        E m_e_obj{};
        double m_threshold{log(1.0 / 0.05)}; // Default threshold ues alpha = 0.05.
        int m_time{0};
        bool m_is_stopped{false};
        int m_stopped_time{kIntInf};
    };

    // Public members
    template <typename E>
    inline void Stcp<E>::updateLogValue(const double x)
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
    inline double Stcp<E>::updateAndReturnHistory(const double x)
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