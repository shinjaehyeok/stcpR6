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
        StcpNormal(const double &threshold,
                   const std::vector<double> &weights,
                   const std::vector<double> &lambdas,
                   const double &mu,
                   const double &sig)
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
            const double &threshold,
            const std::vector<double> &weights,
            const std::vector<double> &lambdas,
            const double &p)
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
            const double &threshold,
            const std::vector<double> &weights,
            const std::vector<double> &lambdas,
            const double &mu)
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
            const double &threshold,
            const double &mu,
            const double &sig,
            const int &window_size)
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
            const double &threshold,
            const double &p,
            const int &window_size)
            : Stcp<GLRCU<L>>::Stcp()
        {
            this->m_threshold = threshold;
            this->m_e_obj = GLRCU<L>(L(p), window_size);
        }
    };

} // End of namespace stcp
#endif