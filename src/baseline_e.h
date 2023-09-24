#ifndef BASELINE_E_H
#define BASELINE_E_H

#include "stcp_interface.h"
#include "baseline_increment.h"

namespace stcp
{
    // Implementation of baseline e-value / e-detector
    template <typename L>
    class BaselineE : public IGeneralE
    {
        static_assert(
            std::is_base_of<IBaselineIncrement, L>::value,
            "Type must be derived from IBaselineIncrement class.");

    public:
        BaselineE()
            : m_log_value{kNegInf}, m_base_obj{}
        {
        }
        BaselineE(const L &base_obj)
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

    // Derived classes from BaselineE
    template <typename L>
    class ST : public BaselineE<L>
    {
    public:
        // ST initalizes log value by 0.0 not -Inf for computational efficiency.
        ST()
            : BaselineE<L>::BaselineE()
        {
            this->m_log_value = 0.0;
        }
        ST(const L &base_obj)
            : BaselineE<L>::BaselineE{base_obj}
        {
            this->m_log_value = 0.0;
        }
        void resetLogValue() override { this->m_log_value = 0.0; }
        void updateLogValue(const double x) override
        {
            this->m_log_value += this->m_base_obj.computeLogBaseValue(x);
        }
    };
    template <typename L>
    class SR : public BaselineE<L>
    {
    public:
        using BaselineE<L>::BaselineE;
        void updateLogValue(const double x) override
        {
            this->m_log_value =
                log(1 + exp(this->m_log_value)) + this->m_base_obj.computeLogBaseValue(x);
        }
    };
    template <typename L>
    class CU : public BaselineE<L>
    {
    public:
        using BaselineE<L>::BaselineE;
        void updateLogValue(const double x) override
        {
            this->m_log_value =
                std::max(0.0, this->m_log_value) + this->m_base_obj.computeLogBaseValue(x);
        }
    };

} // End of namespace stcp
#endif