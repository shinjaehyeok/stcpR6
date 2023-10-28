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
    // Constants and global helper functions
    constexpr double kEps{1e-12};
    constexpr double kNegInf{-std::numeric_limits<double>::infinity()};

    inline double logSumExp(const std::vector<double> &xs)
    {
        if (xs.empty())
        {
            throw std::runtime_error("Empty vector is not allowed for the logSumExp function.");
        }
        double max_x = *std::max_element(xs.begin(), xs.end());
        double sum_exp{0.0};
        for (auto &x : xs)
        {
            sum_exp += std::exp(x - max_x);
        }

        return log(sum_exp) + max_x;
    }

    // stcp interface classes
    class IBaselineIncrement
    {
    public:
        virtual double computeLogBaseValue(const double &x) = 0;

        virtual ~IBaselineIncrement() {}
    };

    class ILogLRIncrement
    {
    public:
        virtual double computeLogBaseValue(const double &x) = 0;
        // Unrestricted mle for h1 parameter
        virtual void updateH1MLE(double &h_1_mle, const double &x, const int &n) = 0; 
        // MaxLLR must be computed under H1 restriction.
        virtual double computeMaxLLR(const double &h_1_mle, const int &n) = 0; 
        
        virtual ~ILogLRIncrement() {}
    };

    class IGeneralE
    {
    public:
        virtual double getLogValue() = 0;
        virtual void reset() = 0;
        virtual void updateLogValue(const double &x) = 0;

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

        virtual void updateLogValue(const double &x) = 0;
        virtual void updateLogValues(const std::vector<double> &xs) = 0;
        virtual void updateLogValuesUntilStop(const std::vector<double> &xs) = 0;

        // For visualization, IStcp support update and return updated history
        virtual double updateAndReturnHistory(const double &x) = 0;
        virtual std::vector<double> updateAndReturnHistories(const std::vector<double> &xs) = 0;

        virtual ~IStcp() {}
    };
} // End of namespace stcp
#endif