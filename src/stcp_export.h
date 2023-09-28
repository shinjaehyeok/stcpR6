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

#ifndef STCP_EXPORT_H
#define STCP_EXPORT_H

#include "stcp.h"

namespace stcp
{
    template <typename E>
    class StcpNormal : public Stcp<MixE<E>>
    {
        static_assert(
            std::is_base_of<ST<Normal>, E>::value ||
                std::is_base_of<SR<Normal>, E>::value ||
                std::is_base_of<CU<Normal>, E>::value,
            "Type must be derived from BaselineE<Normal> class.");

    public:
        StcpNormal()
            : Stcp<MixE<E>>::Stcp()
        {
        }
        StcpNormal(const double threshold,
                   const std::vector<double> &weights,
                   const std::vector<double> &lambdas,
                   const double mu,
                   const double sig)
            : Stcp<MixE<E>>::Stcp()
        {
            this->m_threshold = threshold;
            std::vector<E> e_objs;
            for (auto lambda : lambdas)
            {
                e_objs.push_back(E(Normal(lambda, mu, sig)));
            }
            this->m_e_obj = MixE<E>(e_objs, weights);
        }
    };

    template <typename E>
    class StcpBer : public Stcp<MixE<E>>
    {
        static_assert(
            std::is_base_of<ST<Ber>, E>::value ||
                std::is_base_of<SR<Ber>, E>::value ||
                std::is_base_of<CU<Ber>, E>::value,
            "Type must be derived from BaselineE<Ber> class.");

    public:
        StcpBer()
            : Stcp<MixE<E>>::Stcp()
        {
        }
        StcpBer(
            const double threshold,
            const std::vector<double> &weights,
            const std::vector<double> &lambdas,
            const double p)
            : Stcp<MixE<E>>::Stcp()
        {
            this->m_threshold = threshold;
            std::vector<E> e_objs;
            for (auto lambda : lambdas)
            {
                e_objs.push_back(E(Ber(lambda, p)));
            }
            this->m_e_obj = MixE<E>(e_objs, weights);
        }
    };

    template <typename E>
    class StcpBounded : public Stcp<MixE<E>>
    {
        static_assert(
            std::is_base_of<ST<Bounded>, E>::value ||
                std::is_base_of<SR<Bounded>, E>::value ||
                std::is_base_of<CU<Bounded>, E>::value,
            "Type must be derived from BaselineE<Bounded> class.");

    public:
        StcpBounded()
            : Stcp<MixE<E>>::Stcp()
        {
        }
        StcpBounded(
            const double threshold,
            const std::vector<double> &weights,
            const std::vector<double> &lambdas,
            const double mu)
            : Stcp<MixE<E>>::Stcp()
        {
            this->m_threshold = threshold;
            std::vector<E> e_objs;
            for (auto lambda : lambdas)
            {
                e_objs.push_back(E(Bounded(lambda, mu)));
            }
            this->m_e_obj = MixE<E>(e_objs, weights);
        }
    };

    template <typename L>
    class GLRCUNormal : public Stcp<GLRCU<L>>
    {
        static_assert(
            std::is_base_of<NormalGLR, L>::value ||
                std::is_base_of<NormalGLRGreater, L>::value ||
                std::is_base_of<NormalGLRLess, L>::value,
            "Type must be derived from NormalGLR, NormalGLRGreater or NormalGLRLess.");

    public:
        GLRCUNormal()
            : Stcp<GLRCU<L>>::Stcp()
        {
        }
        GLRCUNormal(
            const double threshold,
            const double mu,
            const double sig,
            const int window_size)
            : Stcp<GLRCU<L>>::Stcp()
        {
            this->m_threshold = threshold;
            this->m_e_obj = GLRCU<L>(L(mu, sig), window_size);
        }
    };

    template <typename L>
    class GLRCUBer : public Stcp<GLRCU<L>>
    {
        static_assert(
            std::is_base_of<BerGLR, L>::value ||
                std::is_base_of<BerGLRGreater, L>::value ||
                std::is_base_of<BerGLRLess, L>::value,
            "Type must be derived from BerGLR, BerGLRGreater or BerGLRLess.");

    public:
        GLRCUBer()
            : Stcp<GLRCU<L>>::Stcp()
        {
        }
        GLRCUBer(
            const double threshold,
            const double p,
            const int window_size)
            : Stcp<GLRCU<L>>::Stcp()
        {
            this->m_threshold = threshold;
            this->m_e_obj = GLRCU<L>(L(p), window_size);
        }
    };

} // End of namespace stcp
#endif