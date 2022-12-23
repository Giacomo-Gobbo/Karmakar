#include <iostream>

#include "funzioni.hpp"

using namespace boost::numeric::ublas;

template <typename T>
struct LinearConstrainSystem
{
    enum class SolutionType
    {
        BOUNDED,  // trovata una soluzione ottima
        UNBOUNDED // l'insieme delle soluzioni non è limitato superiormente
    };

    enum class ConstrainType
    {
        EQ, // ==
        LE, // <=
        GE, // >=
    };

    enum class OptimizationType
    {
        MIN,
        MAX
    };

private:
    matrix<T> A;
    vector<T> b;
    ConstrainType cstType;

public:
    LinearConstrainSystem();

    LinearConstrainSystem(
        matrix<double> &A,
        vector<double> &b,
        ConstrainType &cstType
    ) : A{A}, b{b}, cstType{cstType} {}

    // aggiungi il vincolo a*x type b, e.g., a*x <= b
    LinearConstrainSystem &add_constrain(const std::vector<T> &a, const T &b, const ConstrainType cstType);

    // testa se il sistema di vincoli è soddisfacibile
    bool is_feasible() const;

    // ottimizza c*x rispetto al sistema di vincoli con x che può anche assumere valori negativi.
    // Il metodo itera fino a quando il criterio di stop indicato da un'istanza della classe
    // STOPPING_CRITERIUM restituisce falso.
    // Il metodo restituisce il tipo di soluzione individuata e, se esiste una soluzione bounded,
    // la scrive in solution

    // template <typename STOPPING_CRITERIUM> da implementare
    SolutionType karmakar(vector<T> &solution, const vector<T> &c, const T gamma, const OptimizationType optType)
    {
        u_int repetitions{10}; // da sostituire con STOPPING_CRITERIUM
        u_int k{0};
        vector<double> v(b.size());
        matrix<double> Dv(b.size(), b.size());
        matrix<double> Dv_inversa(b.size(), b.size());
        matrix<double> mInv(b.size(), b.size());
        vector<double> hx(b.size());
        vector<double> hv(b.size());
        double min;
        double tmp;
        bool isUnbounded{true};
        bool isFirst;

        while (k < repetitions)
        {
            v = b - prod(A, solution) + vector<double>(b.size(), 1e-16);
            std::cout << "v[" << k << "] = " << v << " ____ "
                      << "x[" << k << "] = " << solution << " ____ " << std::endl;
            Dv = diagonale(v);
            if (!invertMatrix(Dv, Dv_inversa))
            {
                throw("Matrice Dv non invertibile");
            }
            Dv = prod(Dv_inversa, Dv_inversa);
            Dv = prod(trans(A), Dv);
            Dv = prod(Dv, A);
            if (!invertMatrix(Dv, Dv_inversa))
            {
                throw("Matrice A^tDv^{-2}A non invertibile");
            }
            hx = prod(Dv_inversa, c);
            hv = prod(-A, hx);
            if ((cstType == ConstrainType::LE && optType == OptimizationType::MAX) ||
                (cstType == ConstrainType::GE && optType == OptimizationType::MIN))
            {
                for (uint i{0}; i < hv.size(); ++i)
                {
                    std::cout << i << std::endl;
                    if (hv(i) < 0)
                    {
                        isUnbounded = false;
                    }
                }
            }
            else
            {
                for (uint i{0}; i < hv.size(); ++i)
                {
                    if (hv(i) > 0)
                    {
                        isUnbounded = false;
                    }
                }
            }
            if (isUnbounded)
            {
                return SolutionType::UNBOUNDED;
            }

            isFirst = true;
            for (uint i{0}; i < b.size(); ++i)
            {
                if ((cstType == ConstrainType::LE && ((optType == OptimizationType::MAX && hv(i) < 0) ||
                                                      (optType == OptimizationType::MIN && hv(i) > 0))) ||
                    (cstType == ConstrainType::GE && ((optType == OptimizationType::MAX && hv(i) > 0) ||
                                                      (optType == OptimizationType::MIN && hv(i) < 0))))
                {
                    if (isFirst)
                    {
                        min = -v(i) / hv(i);
                        isFirst = false;
                    }
                    else
                    {
                        tmp = -v(i) / hv(i);
                        if (min > tmp)
                        {
                            min = tmp;
                        }
                    }
                }
            }
            if (k < repetitions - 1)
            {
                solution = solution + min * gamma * hx;
            }

            ++k;
        }
        double res{0};
        for (uint i{0}; i < c.size(); ++i)
        {
            res += c(i) * solution(i);
        }
        std::cout << "Il massimo è: " << res << std::endl;

        return SolutionType::BOUNDED;
    }
};