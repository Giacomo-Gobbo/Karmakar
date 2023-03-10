#include <iostream>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

template <typename T>
class Stopping{
    /**
     * @brief Classe per lo Stopping Criterium
     * 
     * @tparam T: Tipo dei valori dei vettori utilizzati
    */
    std::string type;   //!< Criterio di stop scelto (max_iter o eps)
    long double tol;    //!< Valore della tolleranza

    public:
        /**
         * @brief Costruttore Vuoto
        */
        Stopping()
            : type{"max_iter"}, tol{10}
        {}

        /**
         * @brief Costruttore con parametri
         * 
         * @param type: Criterio di stop scelto
         * @param tol: Valore della tolleranza
        */
        Stopping(std::string type, long double tol)
            : type{type}, tol{tol}
        {
            if (this->type != "max_iter" && this->type != "eps"){
                throw(std::domain_error("Errore nella definizione del criterio di stop: tol_type deve essere 'max_iter' o 'eps'"));
            }
        }

        /**
         * @brief Metodo che controlla se il criterio di stop è verificato
         * 
         * @param par_1: Primo vettore utilizzato per il criterio eps
         * @param par_2: Secondo vettore utilizzato per il criterio eps
         * @param to_check: Numero di iterazioni da usare per il criterio max_iter
        */
        bool check_criterium(vector<T> par_1, vector<T> par_2, const unsigned int& to_check = 0){
            // Se il criterio di stop è max_iter
            if (this->type == "max_iter"){
                // Controllo se il numero di iterazioni effettuate super il numero di iterazioni massime
                if (this->tol <= to_check){
                    return true;
                }
            // Se il criterio di stop è eps
            } else {
                // Controllo se la differenza dei due vettori è minore della tolleranza
                if (norm_1(par_1-par_2) <= this->tol){
                    return true;
                }
            }

            return false;
        }
};

/**
 * @brief Classe template che rappresenta il sistema lineare di vincoli per un problema di programmazione lineare 
 *  
 * @tparam T: Tipo dei coefficienti del sistema lineare
*/
template<typename T>
struct LinearConstrainSystem {
    /**
     * @brief Classe di enumerazione per definire se il sistema ha una soluzione ottima o meno
    */
    enum class SolutionType
    {
        BOUNDED,  //!< Indica che il sistema ha una soluzione ottima 
        UNBOUNDED //!< Indica che il sistema NON ha una soluzione ottima 
    };

    enum class ConstrainType {
        LE, // <=
        EQ, // ==
        GE, // >=
    };

    /**
     * @brief Classe di enumerazione per definire se il sistema è utilizzato per un problema di massimizzazione o minimizzazione
    */
    enum class OptimizationType { 
        MIN, //!< Indica che il problema è di minimizzazione
        MAX  //!< Indica che il problema è di massimizzazione
    };

    matrix<T> A;            //!< Matrice dei coefficienti dei vincoli
    vector<T> b;            //!< Condizioni dei vincoli
    vector<T> slack;        //!< Vettore di variabili di slack
    vector<T> solution;     //!< Vettore soluzione
    uint k;                 //!< Numero di ripetizioni

    /**
     * @brief Costruttore Vuoto
    */
    LinearConstrainSystem()
        : k{0}
    {}
    
    /**
     * @brief Costruttore
     * 
     * @param A: Matrice dei vincoli
     * @param b: Vettore delle condizioni dei vincoli
     * @param slack: Vettore delle variabili di slack iniziali
    */
    LinearConstrainSystem(matrix<T> A, vector<T> b, vector<ConstrainType> slack)
        : A{A}, b{b}, slack{slack}, k{0}
    {}

    /**
     * @brief Metodo che aggiunge un vincolo al sistema
     * 
     * @param a: Vettore dei coefficienti del nuovo vincolo
     * @param b: Condizione del nuovo vincolo
     * @param type: Tipo del nuovo vincolo
    */
    LinearConstrainSystem& add_constrain(const vector<T>& a, const T& b, const ConstrainType type){
        // Calcolo il numero di righe e colonne della matrice dei vincoli
        unsigned long int nrow{A.size1()};
        unsigned long int ncol{A.size2()};

        if(ncol == 0){
            ncol = a.size();
        }
        vector<T> tmp(ncol);
        // Se il nuovo vincolo contiene un numero di variabili maggiore o minore a quelle disponibili elimino quelle in eccesso ed aggiungo quelle in meno impostandole
        // a zero
        if(a.size() != ncol){
            std::cout << "Array a di dimensione errata (" << a.size() << "), imposto a dimensione corretta (" << ncol << ")" << std::endl;

            // Copio il nuovo vettore dei vincoli in quello di dimensione corretta
            for(uint i{0}; i < ncol && i < a.size(); ++i){
                tmp[i] = a[i];
            }
        } else {
            tmp = a;
        }

        // Calcolo le dimensioni della nuova matrice dei vincoli aggiungendo una riga
        unsigned long int newRow{nrow+1};

        // Creo una nuova matrice dei vincoli
        matrix<T> newA(newRow, ncol);
        for (uint i{0}; i < newRow; ++i){
            for (uint j{0}; j < ncol; ++j){
                // Se il coefficiente da copiare appartiene alla matrice precedente, lo copio dalla matrice iniziale
                if (i < nrow){
                    newA(i, j) = A(i, j);
                // Se il coefficiente da copiare corrisponde ad uno del nuovo vincolo, allora lo copio dal vettore passato in input
                } else {
                    newA(i, j) = tmp[j];
                }
            }
        }

        // Aggiungo alle condizioni dei vincoli la nuova condizione
        // E la nuova tipologia
        vector<T> newb(newRow);
        for (uint i{0}; i < newRow; ++i){
            if (i < nrow){
                newb[i] = this->b[i];
            } else if (i >= nrow) {
                newb[i] = b;
            }
        }

        // Aggiungo un valore al vettore di slack
        vector<T> newslack(slack.size()+1);
        for (uint i{0}; i<slack.size()+1; ++i){
            // Copio i valori precedenti
            if (i < slack.size()){
                newslack(i) = slack[i];
            // Se il nuovo vincolo è di minore uguale pongo il nuovo valore ad 1
            } else if(type == ConstrainType::LE){
                newslack(i) = 1;
            // Se il nuovo vincolo è di uguaglianza pongo il nuovo valore a 0
            } else if(type == ConstrainType::EQ){
                newslack(i) = 0;
            // Se il nuovo vincolo è di maggiore uguale pongo il nuovo valore a -1
            } else {
                newslack(i) = -1;
            }
        }

        // Salvo le nuove matrici/vettori
        A = newA;
        this->b = newb;
        slack = newslack;

        return *this;
    }

    /**
     * @brief Metodo che trasforma la matrice A in una matrice che tiene conto dei vincoli
    */
    matrix<T> unione() const{
        // Calcolo il numero di vincoli di uguaglianza
        uint n{0};
        for (uint i{0}; i<slack.size(); ++i) {
            if (slack[i] == 0){
                ++n;
            }
        }

        // La nuova matrice avrà lo stesso numero di righe della matrice originale
        // ed l'unione numero di colonne della matrice A e della matrice diagonale 
        // delle variabili di slack.
        // (a meno delle colonne dei vincoli di uguaglianza poiché vettori nulli)
        unsigned long int Row{A.size1()};
        unsigned long int Col{A.size2()+slack.size()-n};

        // Copio la matrice A
        matrix<T> workA(Row, Col, 0);
        for (uint i{0}; i < A.size1(); ++i){
            for (uint j{0}; j < A.size2(); ++j){
                workA(i, j) = A(i,j);
            }
        }
 
        // Copio le variabili di slack nella diagonale principale
        // (a meno dei vincoli di uguaglianza)
        uint j{0};
        for (uint i{0}; i < slack.size(); ++i){
            if (slack[i] != 0){
                workA(i, j+A.size2()) = slack[i];
                j++;
            }
        }

        return workA;
    }

    /**
     * @brief Controllo se una soluzione trovata soddisfa il sistema o meno
    */
    bool is_feasible() const {
        // Se abbiamo una soluzione
        if(solution.size() != 0){
            // Consideriamo la matrice contenente le variabili di slack
            matrix<T> workA{unione()};
            // Per ogni vincolo presente
            for (int i = 0; i < A.size1(); ++i){
                // Calcoliamo la somma dipendentemente dal vincolo
                T sum{0};
                for (int j = 0; j < A.size2(); ++j){
                    sum += A(i, j)*solution(j); 
                }

                // Se la differenza col vincolo è diversa da zero la soluzione non soddisfa i vincoli
                if (abs(sum - b[i]) >= 1e-10){
                    return false;
                }
            }
            // Se arriviamo alla fine del ciclo i vincoli sono stati rispettati
            return true;
        } else {
            std::cout << "Calcola una soluzione prima!" << std::endl; 
            return false;
        }
    }

    /**
     * @brief Metodo che crea una matrice diagonale con gli elementi del vettore passato in input
     * 
     * @param vector: Vettore da inserire come diagonale principale
     * @return m: Matrice diagonale
    */
    matrix<T> diagonale(vector<T>& vector) const{
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
     * @ref: https://gist.github.com/lilac/2464434
     * 
     * @param input: Matrice da invertire
     * @param inverse: Matrice dove inserire l'inversa
     * @return booleano che indica se la matrice è invertibile o meno
    */
    static bool invertMatrix(const matrix<T>& input, matrix<T>& inverse){
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
        inverse = identity_matrix<T>(A.size1());
        // Inverto la matrice
        lu_substitute(A, pm, inverse);
        return true;
    }

    /**
     * @brief Metodo che controlla se un problema ha soluzione ottimale o meno
     * 
     * @param hv: Array dal quale effettuare il controllo
     * @param opt: Tipo di ottimizzazione selezionato
     * @return unbounded: Booleano che indica se il sistema NON ha soluzione
    */
    bool isUnbounded(const vector<T>& hv, const OptimizationType& opt) const{
        // Controllo se tutti gli elementi sono negativi nel caso della massimizzazione
        if (opt == OptimizationType::MAX){
            auto P = [](T i){ return i<0; };
            return std::all_of(hv.begin(),hv.end(),P);
        // Controllo se tutti gli elementi sono positivi nel caso della minimizzazione
        } else {
            auto P = [](T i){ return i>0; };                
            return std::all_of(hv.begin(),hv.end(),P);
        }
    }

    /**
     * @brief Calcolo la direzione da seguire per il nuovo punto ottimale
     * 
     * @param v: array v
     * @param hv: array hv
     * @param opt: Tipo di soluzione selezionata
     * @return dir: la direzione da seguire
    */
    T direction(const vector<T>& v, const vector<T>& hv, const OptimizationType& opt) const{
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
     * @param solution: Array che conterrà la soluzione
     * @param c: Array contenente i costi del sistema
     * @param gamma: Step-size
     * @param type: Tipo di ottimizzazione scelto
     * @param stop: criterio di Stop selezionato
     * @return Se il sistema ammette una soluzione ottima o meno
    */
    SolutionType karmakar(vector<T>& solution, const vector<T>& c, const T gamma, const OptimizationType type, Stopping<T>& stop){
        // Calcolo la matrice di lavoro A ed aggiungo al vettore dei costi
        // i costi delle variabili di slack (poste a zero)
        matrix<T> workA{unione()};
        vector<T> workc(workA.size2(), 0);
        for (uint i{0}; i<c.size(); ++i){
            workc[i] = c[i];
        }

        std::cout << "workA: " << workA << std::endl;
        std::cout << "workc: " << workc << std::endl << std::endl;

        // Calcolo la soluzione iniziale
        solution = norm_1(b) / norm_1(prod(workA, workc)) * workc;
        vector<T> xprev{solution*2};
    
        while (!stop.check_criterium(xprev, solution, k)) {  
            ++(this->k);
            
            // Calcolo il vettore v che contiene la differenza tra l'array b ed il prodotto della matrice A con il vettore dei punti x
            vector<T> v{b - prod(workA, solution)};
            std::cout << "x[" << (this->k) << "] = " << solution << std::endl;

            // Calcolo la matrice diagonale Dv
            matrix<T> Dv(diagonale(v));
            
            // Se non è invertibile usciamo
            matrix<T> Dv_inversa;
            if (!invertMatrix(Dv, Dv_inversa))
            {
                throw(std::domain_error("Matrice Dv non invertibile"));
            }

            matrix<T> m(trans(workA));
            // (A.T)Dv^-1
            m = prod(m, Dv_inversa);
            // (A.T)Dv^-2
            m = prod(m, Dv_inversa);
            // (A.T)Dv^-2(A)
            m = prod(m, workA);

            matrix<T> mInv;
            // Se non è invertibile usciamo
            if (!invertMatrix(m, mInv))
            {
                throw(std::domain_error("Matrice T(A)Dv^{-2}A non invertibile"));
            }

            // Effettuiamo il prodotto tra la matrice inversa di ((A.T)Dv^-2(A))^-1 e c
            vector<T> hx(prod(mInv, workc));
            // Effettuiamo il prodotto tra la matrice A e la matrice hx
            vector<T> hv(-prod(workA, hx));
            // Se il sistema non ammette una soluzione ottima usciamo

            if (isUnbounded(hv, type))
            {
                return SolutionType::UNBOUNDED;
            }

            // Calcoliamo la direzione da seguire per aggiornare la soluzione ottimale
            T dir = direction(v, hv, type);
            
            // Calcoliamo la nuova soluzione ottimale salvando quella precedente
            xprev = solution;
            solution += dir * gamma * hx;
        }
        this->solution = solution;
        
        // Se arriviamo a questo punto abbiamo trovato una soluzione ottimale
        // Calcoliamo la funzione costo ottenuto
        T res{0};
        for (uint i{0}; i < workc.size(); ++i)
        {
            res += workc(i) * solution(i);
        }

        if (type == OptimizationType::MAX){
            std::cout << std::endl << "Il massimo è: ";
        } else {
            std::cout << std::endl << "Il minimo è: ";
        }
        std::cout << res << std::endl;

        return SolutionType::BOUNDED;
    }
};

int main()
{    
    
    LinearConstrainSystem<long double> obj;
    vector<long double> vett1(2, 1);
    vector<long double> vett2(2, -1);
    vett2[1] = 1;

    obj.add_constrain(vett1, 5, obj.ConstrainType::EQ);
    obj.add_constrain(vett2, 3, obj.ConstrainType::EQ);

    vector<long double> c(2, 1);
    c[0] = 2;

    std::cout << "A: " << obj.A << std::endl;
    std::cout << "b: " << obj.b << std::endl;
    std::cout << "slack: " << obj.diagonale(obj.slack) << std::endl;
    std::cout << "c: " << c << std::endl;

    vector<long double> x0;

    Stopping<long double> stop;
    Stopping<long double> stop2("eps", 1e-11);

    try {

        auto solution = obj.karmakar(x0, c, 0.5, obj.OptimizationType::MAX, stop2);
        if (solution == obj.SolutionType::BOUNDED){
            std::cout << "Bounded" << std::endl;
            std::cout << "Soluzione: " << x0 << std::endl;
        
            if (obj.is_feasible()){
                std::cout << "FEASIBLE" << std::endl;
            } else {
                std::cout << "NOT FEASIBLE" << std::endl;
            }

        } else {
            std::cout << "Unbounded" << std::endl;
        }
    } catch(std::exception& e){
        std::cout << typeid(e).name() << ": " << e.what() << std::endl;
    }

    return 0;
}
