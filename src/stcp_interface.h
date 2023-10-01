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
    constexpr int kIntInf{std::numeric_limits<int>::infinity()};
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
        // Unrestricted mle for h1 parameter
        virtual void updateH1MLE(double &h_1_mle, const double x, const int n) = 0; 
        // MaxLLR must be computed under H1 restriction.
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