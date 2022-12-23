# Struttura LinearCostraintSystem
Una struttura template che rappresenta un sistema di vincoli di un problema di programmazione lineare.
Utilizza l'algoritmo di **Karmakar**, più precisamente la variante dell'**Affine Scaling**, per ottimizzare il problema. È possibile trovare la soluzione ottima, se esiste, che massimizza o minimizza il sistema.
## Attributi
- **A**: Matrice dei coefficienti dei vincoli, ogni riga rappresenta un vincolo ed ogni colonna la singola variabile (ordinate, **coefficiente = 0** se non presente).
- **b**: Vettore contenente le condizioni dei vincoli.
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
`Crea una matrice ed un vettore dei vincoli vuoti ed inizializza a zero il numero di passi effettuati per l'ottimizzazione a zero.`
#### Metodo Costruttore
`Input: matrice A, Vettore b`  
`Output: None`  
`Inizializza la matrice ed il vettore dei vincoli ed inizializza a zero il numero di passi effettuati per l'ottimizzazione a zero.`

### diagonale
`Input: Vettore vector`  
`Output: Matrice m`  
`Crea una matrice diagonale avente nella diagonale principale il vettore passato in input.`

### invertMatrix
`Input: Matrice input, Matrice inverse`  
`Output: Boolean, Matrice inverse`  
`Controlla se la matrice input è invertibile e, nel caso, la inserisce nella matrice inverse. Restituisce un booleano che indica se la matrice è invertibile o meno.`

### diff
`Input: Vettore a, Vettore b`  
`Output: Double differece`  
`Calcola la differenza di due vettori come la somma dei valori assoluti delle differenze dei singoli elementi.`
### isUnbounded
`Input: Vettore hv, OptimizationType opt`  
`Output: Boolean unbounded`  
`Controlla se un sistema è limitato i meno dipendentemente dal tipo di ottimizzazione scelto. Restituisce True se il valore è illimitato, False se ammette una soluzione ottima.`
### direction
`Input: Vettore v, Vettore hv, OptimizationType opt`  
`Output: T dir`  
`Calcola la direzione verso la quale aggiornare la soluzione ottimale dipendentemente dal tipo di ottimizzazione scelto.`
### karmakar
`Input: Vettore solution, Vettore c, T gamma, OptimizationType type`  
`Output: SolutionType, Vettore solution`  
`Effettua l'algoritmo di Karmakar mediante l'affine scaling restituendo se il sistema ammette una soluzione ottimale o meno ed, in caso, il vettore della soluzione.`