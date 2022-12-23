#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

/**
 * @brief Classe template che rappresenta il sistema lineare di vincoli per un problema di programmazione lineare 
 *  
 * @tparam T: Tipo dei coefficienti del sistema lineare
*/
template<typename T>
struct LinearConstrainSystem {
    matrix<T> A;        //!< Matrice dei coefficienti dei vincoli
    vector<T> b;        //!< Condizioni dei vincoli
    uint k;             //!< Numero di ripetizioni

    /**
     * @brief Classe di enumerazione per definire se il sistema ha una soluzione ottima o meno
    */
    enum class SolutionType
    {
        BOUNDED,  //!< Indica che il sistema ha una soluzione ottima 
        UNBOUNDED //!< Indica che il sistema NON ha una soluzione ottima 
    };

    /**
     * @brief Classe di enumerazione per definire se il sistema è utilizzato per un problema di massimizzazione o minimizzazione
    */
    enum class OptimizationType { 
        MIN, //!< Indica che il problema è di minimizzazione
        MAX  //!< Indica che il problema è di massimizzazione
    };

    LinearConstrainSystem()
        : A{nullptr}, b{nullptr}, k{0}
    {}
    
    LinearConstrainSystem(matrix<T> A, vector<T> b)
        : A{A}, b{b}, k{0}
    {}

    /**
     * @brief Metodo che crea una matrice diagonale con gli elementi del vettore passato in input
     * 
     * @param vector: Vettore da inserire come diagonale principale
     * @return m: Matrice diagonale
    */
    matrix<T> diagonale(vector<T> &vector){
        // Creo una matrice quadrata di dimensione pari a quella del vettore in input con tutti gli elementi posti a zero
        matrix<T> m{zero_matrix<T>(vector.size(), vector.size())};

        // Inserisco nella diagonale principale il gli elementi del vettore
        for (size_t i{0}; i < vector.size(); ++i)
        {
            m(i, i) = vector(i);
        }

        return m;
    }

    /**
     * @brief Metodo che inverte la matrice passata in input
     * 
     * @param input: Matrice da invertire
     * @param inverse: Matrice dove inserire l'inversa
     * @return booleano che indica se la matrice è invertibile o meno
    */
    bool invertMatrix(const matrix<T> &input, matrix<T> &inverse){
        // Definisco una matrice di permutazione
        typedef permutation_matrix<std::size_t> pmatrix;
        // Inizializzo una matrice di lavoro
        matrix<T> A(input);
        // Crea una matrice di permutazione per la fattorizzazione LU
        pmatrix pm(A.size1());     
        // Esegue la fattorizzazione LU          
        const int res = lu_factorize(A, pm); 
        // Se la fattorizzazione è nulla allora non è invertibile
        if (res != 0)
        {
            return false;
        }
        // Assegnamo alla matrice inversa il valore della matrice identità
        inverse.assign(identity_matrix<T>(A.size1()));
        // Inverto la matrice
        lu_substitute(A, pm, inverse);
        return true;
    }

    /**
     * @brief Metodo che calcola la differenza fra due array
     * 
     * @param a: primo vettore
     * @param b: secondo vettore
     * @return difference: La differenza tra i due array
    */
    T diff(vector<T> &a, vector<T> &b){
        // Se la dimensione dei due array differisce usciamo
        if (a.size() != b.size()){
                std::cout << "Array di dimensione diversa, uscita..." << std::endl;
                exit(1);
        }
        
        // Calcolo la difference come la somma dei valori assoluti della differenza dei due array
        T difference{0};
        for(uint i{0}; i < a.size(); ++i){
            difference += std::abs(a[i]-b[i]);
        }

        return difference;
    }

    /**
     * @brief Metodo che controlla se un problema ha soluzione ottimale o meno
     * 
     * @param hv: Array dal quale effettuare il controllo
     * @param opt: Tipo di ottimizzazione selezionato
     * @return unbounded: Booleano che indica se il sistema NON ha soluzione
    */
    bool isUnbounded(vector<T> hv, const OptimizationType &opt){
            bool unbounded{true};

            // Controllo se tutti gli elementi...
            for (uint i{0}; i < hv.size(); ++i){
                // ...Sono negativi nel caso della massimizzazione
                if (hv(i) < 0 && opt==OptimizationType::MAX){
                    unbounded = false;
                }
                // ...Sono positivi nel caso della minimizzazione
                else if (hv(i) > 0 && opt==OptimizationType::MIN){
                    unbounded = false;
                }
            }

            return unbounded;
    }

    /**
     * @brief Calcolo la direzione da seguire per il nuovo punto ottimale
     * 
     * @param v: array v
     * @param hv: array hv
     * @param opt: Tipo di soluzione selezionata
     * @return dir: la direzione da seguire
    */
    T direction(vector<T> v, vector<T> hv, const OptimizationType &opt){
        T dir;
        bool first{true};

        // Per ogni elemento dei vettori...
        for (uint i{0}; i < hv.size(); ++i){
            // ...Se il valore di hv è negativo e stiamo massimizzando...
            if (hv(i) < 0 && opt==OptimizationType::MAX){
                T tmp = -v(i) / hv(i);
                // ...Ed il valore trovato è il primo o è minore di quello precedente lo salviamo
                if (dir > tmp || first){
                    if (first){
                        first = false;
                    }
                    dir = tmp;
                } 
            }
            // ...Altrimenti se il valore di hv è positivo e stiamo minimizzando...
            else if (hv(i) > 0 && opt==OptimizationType::MIN){
                T tmp = -v(i) / hv(i);
                // ...Ed il valore trovato è il primo o è maggiore di quello precedente lo salviamo
                if (dir < tmp || first){
                    if (first){
                        first = false;
                    }
                    dir = tmp;
                } 
            }
        }

        return dir;
    }

    /**
     * @brief Metodo che effettua l'algoritmo di karmakar
     * 
     * @tparam STOP_CRITERION: criterio di Stop selezionato
     * @param solution: Array che conterrà la soluzione
     * @param c: Array contenente i costi del sistema
     * @param gamma: Step-size
     * @param type: Tipo di ottimizzazione scelto
     * @return Se il sistema ammette una soluzione ottima o meno
    */
    template <typename STOP_CRITERION>
    SolutionType karmakar(vector<T>& solution, const vector<T>& c, const T gamma, const OptimizationType type, uint repetitions){
        vector<T> xprev{solution*2};
        
        while (((this->k) < repetitions-1) && (diff(xprev, solution) > 1e-16)){
            ++(this->k);
            
            // Calcolo il vettore v che contiene la differenza tra l'array b ed il prodotto della matrice A con il vettore dei punti x
            vector<T> v((this->b).size());
            v = (this->b) - prod((this->A), solution) + vector<T>((this->b).size(), 1e-16);
            std::cout << "v[" << (this->k) << "] = " << v << " ____ "
                    << "x[" << (this->k) << "] = " << solution << " ____ " << std::endl;

            // Calcolo la matrice diagonale Dv
            matrix<T> Dv((this->b).size(), (this->b).size());
            Dv = diagonale(v);
            // Se non è invertibile usciamo
            matrix<T> Dv_inversa((this->b).size(), (this->b).size());
            if (!invertMatrix(Dv, Dv_inversa))
            {
                std::cout << "Matrice Dv non invertibile" << std::endl;
                exit(1);
            }

            matrix<T> m((this->b).size(), (this->b).size());
            matrix<T> mInv((this->b).size(), (this->b).size());
            // Calcoliamo il quadrato della matrice inversa di Dv
            m = prod(Dv_inversa, Dv_inversa);
            // Effettuo il prodotto tra la matrice A trasposta e Dv^-2
            m = prod(trans((this->A)), m);
            // Effettuo il prodotto tra la matrice T(A)Dv^-2 e A
            m = prod(m, (this->A));
            // Se non è invertibile usciamo
            if (!invertMatrix(m, mInv))
            {
                std::cout << "Matrice T(A)Dv^{-2}A non invertibile" << std::endl;
                exit(1);
            }

            // Effettuiamo il prodotto tra la matrice inversa di T(A)Dv^{-2}A e c
            vector<T> hx((this->b).size());
            hx = prod(mInv, c);
            // Effettuiamo il prodotto tra la matrice A e la matrice hx
            vector<T> hv((this->b).size());
            hv = prod(-(this->A), hx);
            // Se il sistema non ammette una soluzione ottima usciamo
            if (isUnbounded(hv, type))
            {
                return SolutionType::UNBOUNDED;
            }

            // Calcoliamo la direzione da seguire per aggiornare la soluzione ottimale
            T dir = direction(v, hv, type);
            
            // Calcoliamo la nuova soluzione ottimale salvando quella precedente
            xprev = solution;
            solution = solution + dir * gamma * hx;
        }
        
        // Se arriviamo a questo punto abbiamo trovato una soluzione ottimale
        T res{0};
        for (uint i{0}; i < c.size(); ++i)
        {
            res += c(i) * solution(i);
        }
        std::cout << "Il massimo è: " << res << std::endl;

        return SolutionType::BOUNDED;
    }
};

int main()
{

    /*    matrix<T> m(2, 2, 0.9);
    m(0, 1) = 0;
    m(1, 0) = 0;
    matrix<T> inv_m(2, 2, 0);
    bool isInvertible{invertMatrix(m, inv_m)};
    std::cout << "m" << m << " ___ inv_m" << inv_m << std::endl;
    */
    matrix<double> A{identity_matrix<double>(2)};
    A(0, 1) = -1;
    A(1, 0) = 1;
    std::cout << "A:" << A << std::endl;
    vector<double> b(2, 2);
    b(1) = 6;
    vector<double> c(2, 1);
    c(1) = 0.5;
    vector<double> x0(2, 0.5);
    x0(1) = 0;
    uint repetitions{10};
    LinearConstrainSystem<double> obj(A, b);
    auto solution = obj.karmakar<int>(x0, c, 0.5, obj.OptimizationType::MAX, 50);
    if (solution == obj.SolutionType::BOUNDED){
        std::cout << "Bounded" << std::endl;
        std::cout << "Soluzione: " << x0 << std::endl;
    } else {
        std::cout << "Unbounded" << std::endl;
    }
    // if (affineScaling(A, b, c, x0, OptimizationType::MIN, repetitions, 0.5) == SolutionType::BOUNDED)
    // {
    //     std::cout << "BOUNDED solution" << std::endl;
    // }
    // else
    // {
    //     std::cout << "UNBOUNDED solution" << std::endl;
    // }

    return 0;
}