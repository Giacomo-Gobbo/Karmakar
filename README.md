# Struttura Stopping
Una struttura template che rappresenta il criterio di stop usato per l'algoritmo di karmakar.

## Attributi
- **type**: Stringa con criterio di stop scelto. (max_iter o eps)
- **tol**: Valore per la tolleranza

#### Metodo Costruttore Vuoto
`Input: None`  
`Output: None`
Seleziona il criterio di stop avente un numero di iterazioni massimo pari a 10.

#### Metodo Costruttore
`Input: String type, longo double tol`  
`Output: None`  
Seleziona il criterio di stop e relativo vincolo. Se il criterio di stop non è tra quelli previsti lancia un errore.

### add_costrain
`Input: Vettore par_1, Vettore par_2, unsigned int to_check`  
`Output: Boolean`  
Controlla se il criterio di stop è verificato.


# Struttura LinearCostraintSystem
Una struttura template che rappresenta un sistema di vincoli di un problema di programmazione lineare.
Utilizza l'algoritmo di **Karmakar**, più precisamente la variante dell'**Affine Scaling**, per ottimizzare il problema. È possibile trovare la soluzione ottima, se esiste, che massimizza o minimizza il sistema.

## Attributi
- **A**: Matrice dei coefficienti dei vincoli, ogni riga rappresenta un vincolo ed ogni colonna la singola variabile (ordinate, **coefficiente = 0** se non presente).
- **b**: Vettore contenente le condizioni dei vincoli.
- **slack**: Vettore contenente le variabili di slack.
- **solution**: Vettore della soluzione.
- **k**: Numero di iterazioni effettuate per trovare la soluzione ottima.

## Enum Class
### SolutionType
Indica il tipo di soluzione trovata:
- **BOUNDED**: Il sistema ha soluzione ottima.
- **UNBOUNDED**: Il sistema è illimitato.

### OptimizationType
Indica il tipo di ottimizzazione applicare:
- **MIN**: Ricerca della soluzione che minimizza il sistema.
- **MAX**: Ricerca della soluzione che massimizza il sistema.



## Metodi
### LinearConstrainSystem
#### Metodo Costruttore Vuoto
`Input: None`  
`Output: None`  
Crea una matrice ed un vettore dei vincoli vuoti ed inizializza a zero il numero di passi effettuati per l'ottimizzazione a zero.

#### Metodo Costruttore
`Input: matrice A, Vettore b, Vettore slack`  
`Output: None`  
Inizializza la matrice ed i vettori dei vincoli e di slack. Inizializza a zero il numero di passi effettuati per l'ottimizzazione a zero.

### add_costrain
`Input: Vettore a, T b, ConstrainType type`  
`Output: LinearConstrainSystem`  
Aggiunge un nuovo vincolo nella matrice dei vincoli, aggiungendo anche le variabili di slack (se necessario).

### unione
`Input: None`  
`Output: boolean`  
Controllo se la soluzione è soddisfacibile controllando se la differenza tra l'elemento b del vincolo e della combinazione lineare dei vettori a\*x è inferiore ad una certa tolleranza.

### is_feasible
`Input: None`  
`Output: Matrice workA`  
Aggiungo alla matrice dei vincoli le colonne della matrice degli slack (non aggiungo vincoli di uguaglianza perché colonne formate da vettori nulle).

### diagonale
`Input: Vettore vector`  
`Output: Matrice m`  
Crea una matrice diagonale avente nella diagonale principale il vettore passato in input.

### invertMatrix
`Input: Matrice input, Matrice inverse`  
`Output: Boolean, Matrice inverse`  
Controlla se la matrice input è invertibile e, nel caso, la inserisce nella matrice inverse. Restituisce un booleano che indica se la matrice è invertibile o meno.

### isUnbounded
`Input: Vettore hv, OptimizationType opt`  
`Output: Boolean unbounded`  
Controlla se un sistema è limitato i meno dipendentemente dal tipo di ottimizzazione scelto. Restituisce True se il valore è illimitato, False se ammette una soluzione ottima.

### direction
`Input: Vettore v, Vettore hv, OptimizationType opt`  
`Output: T dir`  
Calcola la direzione verso la quale aggiornare la soluzione ottimale dipendentemente dal tipo di ottimizzazione scelto.

### karmakar
`Input: Vettore solution, Vettore c, T gamma, OptimizationType type, Stopping stop`  
`Output: SolutionType, Vettore solution`  
Effettua l'algoritmo di Karmakar mediante l'affine scaling restituendo se il sistema ammette una soluzione ottimale o meno ed, in caso, il vettore della soluzione. Si arresta se il criterio di stop è soddisfatto
