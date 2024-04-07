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
        double getTime() override { return m_time; };
        double getStoppedTime() override { return m_stopped_time; }

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

        void updateLogValueByAvg(const double &x_bar, const double &n) override;
        void updateLogValuesByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns) override;
        void updateLogValuesUntilStopByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns) override;

        double updateAndReturnHistory(const double &x) override;
        std::vector<double> updateAndReturnHistories(const std::vector<double> &xs) override;

        double updateAndReturnHistoryByAvg(const double &x_bar, const double &n) override;
        std::vector<double> updateAndReturnHistoriesByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns) override;

    protected:
        E m_e_obj{};
        double m_threshold{log(1.0 / 0.05)}; // Default threshold ues alpha = 0.05.
        double m_time{0.0};
        bool m_is_stopped{false};
        double m_stopped_time{0.0};
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
    inline void Stcp<E>::updateLogValueByAvg(const double &x_bar, const double &n)
    {
        m_e_obj.updateLogValueByAvg(x_bar, n);
        m_time += n;
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
    inline void Stcp<E>::updateLogValuesByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns)
    {
        if (x_bars.size() != ns.size())
        {
            throw std::runtime_error("x_bars and ns do not have the same length.");
        }
        for (std::size_t i = 0; i < x_bars.size(); i++)
        {
            this->updateLogValueByAvg(x_bars[i], ns[i]);
        }
    }
    template <typename E>
    inline void Stcp<E>::updateLogValuesUntilStopByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns)
    {
        if (x_bars.size() != ns.size())
        {
            throw std::runtime_error("x_bars and ns do not have the same length.");
        }
        for (std::size_t i = 0; i < x_bars.size(); i++)
        {
            this->updateLogValueByAvg(x_bars[i], ns[i]);
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
    inline double Stcp<E>::updateAndReturnHistoryByAvg(const double &x_bar, const double &n)
    {
        this->updateLogValueByAvg(x_bar, n);
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
    template <typename E>
    inline std::vector<double> Stcp<E>::updateAndReturnHistoriesByAvgs(const std::vector<double> &x_bars, const std::vector<double> &ns)
    {
        if (x_bars.size() != ns.size())
        {
            throw std::runtime_error("x_bars and ns do not have the same length.");
        }
        std::vector<double> log_values(x_bars.size());
        for (std::size_t i = 0; i < x_bars.size(); i++)
        {
            log_values[i] = this->updateAndReturnHistoryByAvg(x_bars[i], ns[i]);
        }
        return log_values;
    }
} // End of namespace stcp
#endif