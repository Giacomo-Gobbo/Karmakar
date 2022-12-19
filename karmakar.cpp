#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

enum class SolutionType
{
    BOUNDED,  // trovata una soluzione ottima
    UNBOUNDED // l'insieme delle soluzioni non è limitato superiormente
};

/**
 * Funzione per ottenere una matrice diagonale da un vettore
 *
 * La matrice diagonale avrà in posizione (i, i) il valore vector(i)
 *
 * @param vector è il vettore i cui valori devono formare la diagonale della matrice
 * @return è la matrica diagonale
 */
matrix<double> diagonale(vector<double> &vector);

/**
 * Funzione per ottenere il quadrato dell'inversa di una matrice
 *
 * @param input è la matrice su cui si vuole compiere l'inversione
 * @param inverse è la matrice inversa
 * @return è una condizione di verità che è vera se la matrice è invertibile
 */
bool invertMatrix(const matrix<double> &input, matrix<double> &inverse);

/**
 * Versione del algoritmo di Karmakar con l'utilizzo delle trasformazioni affini
 * Non è un algoritmo di complessità polinomiale
 *
 * @param A matrice di Ax <= b
 * @param b
 * @param c vettore di cTx (cT è il vettore riga di c)
 * @param x0
 * @param repetitions è il numero di iterazioni consentite per ottenere il risultato
 * @return è il vettore massimo
 */
SolutionType affineScaling(
    matrix<double> &A,
    vector<double> &b,
    vector<double> &c,
    vector<double> &x0,
    uint &repetitions);

int main()
{
    using namespace boost::numeric::ublas;

    matrix<double> A{identity_matrix<double>(2)};
    vector<double> b(2, 1);
    vector<double> c(2, 2);
    c(1) = 1;
    vector<double> x0(2, 0.1);
    uint repetitions{10};
    if (affineScaling(A, b, c, x0, repetitions) == SolutionType::BOUNDED)
    {
        std::cout << "BOUNDED solution" << std::endl;
    }
    else
    {
        std::cout << "UNBOUNDED solution" << std::endl;
    }

    return 0;
}

matrix<double> diagonale(vector<double> &vector)
{
    matrix<double> m(vector.size(), vector.size());

    for (size_t i{0}; i < vector.size(); ++i)
    {
        m(i, i) = vector(i);
    }

    return m;
}

bool invertMatrix(const matrix<double> &input, matrix<double> &inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<double> A(input);             // copia di lavoro della matrice input
    pmatrix pm(A.size1());               // crea una matrice di permutazione per la fattorizzazione LU
    const int res = lu_factorize(A, pm); // esegue la fattorizzazione LU
    if (res != 0)
    {
        return false;
    }
    inverse.assign(identity_matrix<double>(A.size1())); // crea la matrice identità
    lu_substitute(A, pm, inverse);                      // backsubstitute to get the inverse
    return true;
}

SolutionType affineScaling(matrix<double> &A, vector<double> &b, vector<double> &c, vector<double> &x0, uint &repetitions)
{
    u_int k{0};
    vector<double> v[repetitions];
    vector<double> x[repetitions];
    matrix<double> Dv(b.size(), b.size());
    matrix<double> Dv_inversa(b.size(), b.size());
    matrix<double> mInv(b.size(), b.size());
    vector<double> hx(b.size());
    vector<double> hv(b.size());
    // double alpha;
    double min{0};
    double tmp;
    // double gamma{1};

    x[0] = x0;

    while (k < repetitions)
    {
        v[k] = b - prod(A, x[k]);
        Dv = diagonale(b);
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
        for (uint i{0}; i < hv.size(); ++i)
        {
            if (hv(i) > 0)
            {
                return SolutionType::UNBOUNDED;
            }
        }
        for (uint i{0}; i < b.size() - 1; ++i)
        {
            if (hv(i) < 0)
            {
                do
                {
                    min = -v[k](i) / hv(i);
                } while (false);
                tmp = -v[k](i + 1) / hv(i + 1);
                if (min > tmp)
                {
                    min = tmp;
                }
            }
        }
        // alpha = gamma*min;
        x[k + 1] = x[k] + min * hx;

        ++k;
    }
    double res{0};
    for (uint i{0}; i < c.size(); ++i) {
        res += c(i)*x[repetitions - 1](i);
    }
    std::cout << "Il massimo è: " << res << std::endl;

    return SolutionType::BOUNDED;
}
