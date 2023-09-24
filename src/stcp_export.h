#ifndef STCP_EXPORT_H
#define STCP_EXPORT_H

#include "stcp.h"

namespace stcp
{
    template <typename E>
    class MixENormal : public MixE<E>
    {
        static_assert(
            std::is_base_of<ST<Normal>, E>::value ||
                std::is_base_of<SR<Normal>, E>::value ||
                std::is_base_of<CU<Normal>, E>::value,
            "Type must be derived from IGeneralE class.");

    public:
        MixENormal()
            : MixE<E>::MixE()
        {
        }
        MixENormal(const std::vector<double> &weights,
                   const std::vector<double> &lambdas,
                   const double mu,
                   const double sig)
            : MixE<E>::MixE()
        {
            this->m_weights = weights;
            this->m_log_weights = this->validateAndComputeLogWeights(weights);
            this->m_e_objs.clear();
            for (auto weight : weights)
            {
                this->m_weights.push_back(weight);
            }
            for (auto lambda : lambdas)
            {
                this->m_e_objs.push_back(E(Normal(lambda, mu, sig)));
            }
        }
    };

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
        StcpNormal(const std::vector<double> &weights,
                   const std::vector<double> &lambdas,
                   const double mu,
                   const double sig,
                   const double threshold)
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
} // End of namespace stcp
#endif