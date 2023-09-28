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

#ifndef LOG_LR_E_H
#define LOG_LR_E_H

#include "stcp_interface.h"

namespace stcp
{
    // Implementation of LR based methods
    template <typename L>
    class LogLRE : public IGeneralE
    {
        static_assert(
            std::is_base_of<ILogLRIncrement, L>::value,
            "Type must be derived from ILogLRIncrement class.");

    public:
        LogLRE()
            : m_log_value{kNegInf}, m_base_obj{}
        {
        }
        LogLRE(const L &base_obj)
            : m_log_value{kNegInf}, m_base_obj{base_obj}
        {
        }

        double getLogValue() override { return m_log_value; }
        void resetLogValue() override { m_log_value = kNegInf; }
        virtual void updateLogValue(const double x) override = 0;

    protected:
        double m_log_value;
        L m_base_obj;
    };

    template <typename L>
    class GLRCU : public LogLRE<L>
    {
    public:
        using LogLRE<L>::LogLRE;
        GLRCU()
            : LogLRE<L>::LogLRE()
        {
        }
        GLRCU(const int window_size)
            : LogLRE<L>::LogLRE(), m_window_size{window_size}
        {
        }
        GLRCU(const L &base_obj, const int window_size)
            : LogLRE<L>::LogLRE(base_obj), m_window_size{window_size}
        {
        }

        void updateLogValue(const double x) override
        {
            if (static_cast<int>(m_h1_mle_history.size()) >= m_window_size)
            {
                m_h1_mle_history.pop_back();
            }

            m_h1_mle_history.push_front(0.0);
            int n{0};
            double max_log_value{kNegInf};
            for (auto &h_1_mle : m_h1_mle_history)
            {
                n++; // n: sample size *after" taking account for new observation x
                this->m_base_obj.updateH1MLE(h_1_mle, x, n);
                double log_value_new{this->m_base_obj.computeMaxLLR(h_1_mle, n)};
                if (max_log_value < log_value_new)
                {
                    max_log_value = log_value_new;
                }
            }
            this->m_log_value = max_log_value;
        }

    private:
        std::deque<double> m_h1_mle_history;
        int m_window_size{100};
    };

} // End of namespace stcp
#endif