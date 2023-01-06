#include <iostream>
#include <string>
#include <exception>
#include <cstdlib>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/qvm/mat_operations.hpp>

using namespace boost::numeric::ublas;

template <typename T>
class STOPPING_CRITERIUM{ // definizione della classe StoppingCriterium
    public:
    std::string tol_type{"max_iter"}; // tol_type deve essere 'max_iter' o 'eps'
    T tol{10}; // di default è max_iter con 10 iterazioni 

    void check_criterium(const std::string tol_type){ // funzione che controlla se l'input del tipo del criterio è corretto
        if (tol_type != "max_iter" && tol_type != "eps"){
            std::cout << "Errore nella definizione del criterio di stop: tol_type deve essere 'max_iter' o 'eps'" << std::endl;
            throw("");
        }
    }
};

template <typename T> // metodo per controllare se è rispettato il criterio max_iter oppure il criterio eps
bool check_criterium(const STOPPING_CRITERIUM<T>& criterium, const unsigned int& to_check = 0, const double& par_1 = 0, const double& par_2 = 0){
    if (criterium.tol_type == "max_iter"){
        if (criterium.tol-1 < to_check){
            return false;
        }
        return true;
    } 
 
    else if (criterium.tol_type == "eps"){
        if (abs(par_1-par_2) < criterium.tol){
            return false;
        }
        return true;
    }

   return false;
}

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

template <typename T>
SolutionType affineScaling(
    matrix<double> &A,
    vector<double> &b,
    vector<double> &c,
    vector<double> &x0,
    unsigned int &repetitions,
    OptimizationType opt,
    ConstraintType con,
    STOPPING_CRITERIUM<T> stop
    );

};

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

template <typename T>
LinearConstraintSystem::SolutionType LinearConstraintSystem::affineScaling(matrix<double> &A, vector<double> &b, vector<double> &c, vector<double> &x0, unsigned int &repetitions, LinearConstraintSystem::OptimizationType opt, LinearConstraintSystem::ConstraintType con, STOPPING_CRITERIUM<T> stop){
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
    vector<double> for_storing(repetitions); // vettore ausiliario per il criterio 'eps'
    double alpha;
    double min{0};
    double tmp;
    double gamma{0.5};

    x[0] = x0;

    double par_1 = 0;
    double par_2 = 10e5;

    while (k < repetitions && check_criterium(stop/*istanza della classe StoppingCriterium*/, k, par_1/*x[k]*/, par_2/*x[k-1]*/)==true)
    {
        std::cout << "\n\nSTARTING ITERATION n. " << k+1 << std::endl;

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
            temp_sum += c(i)*x[k](i);
        }
        std::cout << "\nGuess: " << temp_sum << std::endl;
        for_storing[k] = temp_sum;
        if(k>=1){
            par_1 = for_storing[k];
            par_2 = for_storing[k-1];
        }
        ++k;
    }

    return SolutionType::BOUNDED;
}

int main(){
    using namespace boost::numeric::ublas;

    unsigned int repetitions{24};

    // criterio di stop:
    STOPPING_CRITERIUM <int> crit;
    crit.tol_type = "eps";
    crit.tol = 1;

    check_criterium(crit);

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

    if (lcs.affineScaling(A, b, c, x0, repetitions, opt, con, crit) == LinearConstraintSystem::SolutionType::BOUNDED)
    {
        std::cout << "BOUNDED solution" << std::endl;
    }
    else
    {
        std::cout << "UNBOUNDED solution" << std::endl;
    }

    return 0;
}
