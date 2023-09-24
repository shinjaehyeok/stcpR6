#ifndef STCP_INTERFACE_H
#define STCP_INTERFACE_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <deque>

namespace stcp
{
    constexpr double kEps{1e-12};
    constexpr double kNegInf{-std::numeric_limits<double>::infinity()};

    class IBaselineIncrement
    {
    public:
        virtual double computeLogBaseValue(const double x) = 0;

        virtual ~IBaselineIncrement() {}
    };

    class ILogLRIncrement
    {
    public:
        virtual double computeLogBaseValue(const double x) = 0;
        virtual void updateH1MLE(double &h_1_mle, const double x, const int n) = 0;
        virtual double computeMaxLLR(const double h_1_mle, const int n) = 0;
        
        virtual ~ILogLRIncrement() {}
    };

    class IGeneralE
    {
    public:
        virtual double getLogValue() = 0;
        virtual void resetLogValue() = 0;
        virtual void updateLogValue(const double x) = 0;

        virtual ~IGeneralE() {}
    };

    class IStcp
    {
    public:
        virtual double getLogValue() = 0;
        virtual double getThreshold() = 0;

        virtual bool isStopped() = 0;
        virtual int getTime() = 0;
        virtual int getStoppedTime() = 0;

        virtual void reset() = 0;

        virtual void updateLogValue(const double x) = 0;
        virtual void updateLogValues(const std::vector<double> &xs) = 0;
        virtual void updateLogValuesUntilStop(const std::vector<double> &xs) = 0;

        // For visualization, IStcp support update and return updated history
        virtual double updateAndReturnHistory(const double x) = 0;
        virtual std::vector<double> updateAndReturnHistories(const std::vector<double> &xs) = 0;

        virtual ~IStcp() {}
    };
} // End of namespace stcp
#endif