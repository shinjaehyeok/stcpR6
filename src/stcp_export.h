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
            "Type must be derived from BaselineE<Normal> class.");

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

    template <>
    class StcpBer<GLRCU<BerGLR>> : public Stcp<GLRCU<BerGLR>>
    {
    public:
        StcpBer<GLRCU<BerGLR>>()
            : Stcp<GLRCU<BerGLR>>::Stcp()
        {
        }
        StcpBer<GLRCU<BerGLR>>(
            const double threshold,
            const double p,
            const int window_size)
            : Stcp<GLRCU<BerGLR>>::Stcp()
        {
            this->m_threshold = threshold;
            this->m_e_obj = GLRCU<BerGLR>(BerGLR(p), window_size);
        }
    };
} // End of namespace stcp
#endif