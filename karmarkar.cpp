#include <iostream>
#include <exception>
#include <cstdlib>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/qvm/mat_operations.hpp>

using namespace boost::numeric::ublas;

struct LinearConstraintSystem{

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

enum class ConstraintType{
    EQ, // ==
    LE, // <=
    GE, // >=
};

enum class OptimizationType { MIN, MAX };

LinearConstraintSystem();

LinearConstraintSystem& add_constraint(const vector<double>& a, const double& b, const ConstraintType type);

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
    unsigned int &repetitions,
    OptimizationType opt,
    ConstraintType con);

};

LinearConstraintSystem& LinearConstraintSystem::add_constraint(const vector<double>& a, const double& b, const LinearConstraintSystem::ConstraintType type){

}

matrix<double> LinearConstraintSystem::diagonale(vector<double> &vector)
{
    matrix<double> m(vector.size(), vector.size());

    for (size_t i{0}; i < vector.size(); ++i){
        for (size_t j{0}; j < vector.size(); ++j){
            m(i,j) = 0;
        }
    }
    for (size_t i{0}; i < vector.size(); ++i)
    {
        m(i, i) = vector(i);
    }

    return m;
}

bool LinearConstraintSystem::invertMatrix(const matrix<double> &input, matrix<double> &inverse)
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

LinearConstraintSystem::SolutionType LinearConstraintSystem::affineScaling(matrix<double> &A, vector<double> &b, vector<double> &c, vector<double> &x0, unsigned int &repetitions, LinearConstraintSystem::OptimizationType opt, LinearConstraintSystem::ConstraintType con){
    unsigned int k{0};
    vector<double> v[repetitions];
    vector<double> x[repetitions];
    matrix<double> Dv(b.size(), b.size());
    matrix<double> Dv_inversa(b.size(), b.size());
    matrix<double> Dv_inversa_squared(b.size(), b.size());
    matrix<double> Dv_toinvert(b.size(), b.size());
    matrix<double> Dvinvert(b.size(), b.size());
    matrix<double> mInv(b.size(), b.size());
    vector<double> hx(b.size());
    vector<double> hv(b.size());
    double alpha;
    double min{0};
    double tmp;
    double gamma{0.5};

    x[0] = x0;

    while (k < repetitions)
    {
        
        std::cout << "\n\nSTARTING ITERATION n. " << k << std::endl;

        v[k] = b - prod(A, x[k]);

        std::cout << "x[k]: " << x[k] << std::endl;

        Dv = diagonale(v[k]);


        if (!invertMatrix(Dv, Dv_inversa))
        {
            std::cout << "Matrice non invertibile --> Matrice: " << Dv << std::endl;
            throw("Matrice A^tDv^{-2}A non invertibile");
        };

        Dv_inversa_squared = prod(Dv_inversa, Dv_inversa);
        Dv_inversa_squared = prod(trans(A), Dv_inversa_squared);
        Dv_inversa_squared = prod(Dv_inversa_squared, A);

        if (!invertMatrix(Dv_inversa_squared, Dv_inversa))
        {
            std::cout << "Matrice non invertibile --> Matrice: " << Dv << std::endl;
            throw("Matrice A^tDv^{-2}A non invertibile");
        };


        hx = prod(Dv_inversa, c);
        
        hv = prod(-A, hx);
        
        unsigned int counter{0};

        for (unsigned int i{0}; i < hv.size(); ++i)
        {
            if (con==ConstraintType::LE){
                if (opt==OptimizationType::MIN){

                    if (hv(i) < 0) // min: <, max: >
                    {   
                        ++counter;
                    }
                }

                else if (opt==OptimizationType::MAX){
                    if (hv(i) > 0) // min: <, max: >
                    {   
                        ++counter;
                    }
                }
            }

            else if (con==ConstraintType::GE){
                if (opt==OptimizationType::MIN){

                    if (hv(i) > 0) // min: <, max: >
                    {   
                        ++counter;
                    }
                }

                else if (opt==OptimizationType::MAX){
                    if (hv(i) < 0) // min: <, max: >
                    {   
                        ++counter;
                    }
                }
            }

        }

        if (counter==hv.size()){
            std::cout << "UNBOUNDED" << std::endl;
            return SolutionType::UNBOUNDED;
        }


        // TROVARE IL MINIMO:


        bool isFirst = true;
        for (unsigned int i{0}; i < b.size(); ++i)
        {
            if (con==ConstraintType::LE){
                if (opt==OptimizationType::MIN){
                    if (hv(i) > 0) // min: >, max: <
                    {
                        if (isFirst)
                        {
                            min = -v[k][i] / hv[i];
                            isFirst = false;
                        }
                        else
                        {
                            tmp = -v[k][i] / hv[i];
                            if (min > tmp)
                            {
                                min = tmp;
                            }
                        }
                    }
                }
                else if (opt==OptimizationType::MAX){
                    if (hv(i) < 0) // min: >, max: <
                    {
                        if (isFirst)
                        {
                            min = -v[k][i] / hv[i];
                            isFirst = false;
                        }
                        else
                        {
                            tmp = -v[k][i] / hv[i];
                            if (min > tmp)
                            {
                                min = tmp;
                            }
                        }
                    }
                }
            }
            else if (con==ConstraintType::GE){
                if (opt==OptimizationType::MIN){
                    if (hv(i) < 0) // min: >, max: <
                    {
                        if (isFirst)
                        {
                            min = -v[k][i] / hv[i];
                            isFirst = false;
                        }
                        else
                        {
                            tmp = -v[k][i] / hv[i];
                            if (min > tmp)
                            {
                                min = tmp;
                            }
                        }
                    }
                }
                else if (opt==OptimizationType::MAX){
                    if (hv(i) > 0) // min: >, max: <
                    {
                        if (isFirst)
                        {
                            min = -v[k][i] / hv[i];
                            isFirst = false;
                        }
                        else
                        {
                            tmp = -v[k][i] / hv[i];
                            if (min > tmp)
                            {
                                min = tmp;
                            }
                        }
                    }
                }
            }
        }
        x[k + 1] = x[k] + gamma * min * hx;

        double temp_sum{0};
        for(unsigned int i{0}; i<c.size(); ++i){
            temp_sum += c(i)*x[k+1](i);
        }
        std::cout << "\nGuess: " << temp_sum;
        ++k;

    }
    double res{0};
    for (unsigned int i{0}; i < c.size(); ++i) {
        res += c(i)*x[repetitions - 1](i);
    }
    std::cout << "\n\nOptimum: " << res << std::endl;

    return SolutionType::BOUNDED;
}

int main(){

    using namespace boost::numeric::ublas;

    unsigned int repetitions{24};

    // inizializzazione:
    matrix<double> A(2, 2);
    vector<double> b(2);
    vector<double> c(2);
    vector<double> x0(2);
    LinearConstraintSystem::OptimizationType opt{LinearConstraintSystem::OptimizationType::MIN};
    LinearConstraintSystem::ConstraintType con{LinearConstraintSystem::ConstraintType::GE};

    // TEST 1:
    A(0,0) = 1;
    A(0,1) = 0;
    A(1,0) = 0;
    A(1,1) = 1;
    b[0] = 5;
    b[1] = 5;
    c[0] = 1;
    c[1] = 1;
    x0[0] = 10;
    x0[1] = 10;

    LinearConstraintSystem lcs;

    if (lcs.affineScaling(A, b, c, x0, repetitions, opt, con) == LinearConstraintSystem::SolutionType::BOUNDED)
    {
        std::cout << "BOUNDED solution" << std::endl;
    }
    else
    {
        std::cout << "UNBOUNDED solution" << std::endl;
    }

    return 0;
}
