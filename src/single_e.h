#ifndef SINGLE_E_H
#define SINGLE_E_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

namespace stcp
{
    // Implementation of E-value / detector
    class ISingleE
    {
    public:
        virtual double getLogValue() = 0;
        virtual void reinitialize(const double log_value, const double lambda, const double s_param, const double v_param) = 0;
        virtual void updateLogValue(const double x) = 0;

        virtual ~ISingleE() {}
    };

    template <typename L>
    class SingleE : public ISingleE
    {
        static_assert(std::is_base_of<BaselineIncrement, L>::value, "Type must be derived from BaselineIncrement class.");

    public:
        SingleE()
            : m_log_value{0.0}, m_base_obj{}
        {
        }
        SingleE(const double log_value, const double lambda, const double s_param, const double v_param)
            : m_log_value{log_value}, m_base_obj{lambda, s_param, v_param}
        {
        }

        double getLogValue() override { return m_log_value; }
        void reinitialize(const double log_value, const double lambda, const double s_param, const double v_param) override
        {
            m_log_value = log_value;
            m_base_obj.reinitialize(lambda, s_param, v_param);
        }

        virtual void updateLogValue(const double x) override = 0;

    protected:
        double m_log_value;
        L m_base_obj;
    };
    // Derived classes from SingleE
    template <typename L>
    class ST : public SingleE<L>
    {
    public:
        using SingleE<L>::SingleE;
        void updateLogValue(const double x) override
        {
            this->m_log_value += this->m_base_obj.computeLogBaseValue(x);
        }
    };

} // End of namespace stcp
#endif