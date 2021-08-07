# [Schelling's-Model-of-Segregation](http://nifty.stanford.edu/2014/mccown-schelling-model-segregation/)

| **Petrazzuolo Lorenzo** | **0522500894** | **23/08/2021** |
| --- | --- | --- |

Progetto per l'esame di _Programmazione Concorrente, Parallela e sul Cloud_ dell'anno di corso _2020/2021_. 
Laurea magistrale in Computer Science, curriculum in Cloud Computing.
___
## **Indice**

- [**Introduzione**](#introduzione)
  - [_Il modello di segregazione di Schelling_](#il-modello-di-segregazione-di-Schelling)
  - [_Descrizione del problema_](#descrizione-del-problema)

- [**Dettagli implementativi**](#dettagli-implementativi)
  - [_Creazione Matrice_](#creazione-matrice)
  - [_Distribuzione Matrice_](#distribuzione-matrice)
  - [_Calcolo soddisfazione_](#calcolo-soddisfazione)
  - [_Riposizionamento insoddisfatti_](#riposizionamento-insoddisfatti)
  - [_Ricomposizione Matrice_](#ricomposizione-matrice)
  - [_Stampa dei risultati_](#stampa-dei-risultati)

 - [**Note sull'implementazione**](#note_sullimplementazione)
  - [_Compilazione_](#compilazione)
  - [_Esecuzione_](#esecuzione)

- [**Benchmarks**](#benchmarks)

- [**Correttezza**](#correttezza)

- [**Conclusioni**](#conclusioni)


___

## Introduzione 

### _Il modello di segregazione di Schelling_

Nel 1971, l'economista americano Thomas Schelling creò un modello basato su agenti che suggeriva che anche il comportamento involontario potesse contribuire alla segregazione. Il suo modello di segregazione ha mostrato che anche quando gli individui (o "agenti") non si preoccupavano di essere circondati o vivere da agenti di una razza o un background economico diverso, avrebbero comunque scelto di separarsi dagli altri agenti nel tempo. Sebbene il modello sia abbastanza semplice, fornisce uno sguardo affascinante su come gli individui potrebbero auto-segregarsi, anche quando non hanno un desiderio esplicito di farlo.
> _[Riferimenti](https://en.wikipedia.org/wiki/Schelling's_model_of_segregation)_


### _Descrizione del problema_

Nel seguente progetto si propone un'implementazione di una simulazione del modello di Shelling. Una volta impostati una serie di parametri, il programma dovrà essere in grado di descrivere graficamente le caratteristiche del modello. Più precisamente si suppone di avere due tipi di agenti _**X**_ e _**O**_ che andranno a comporre la popolazione di una matrice _**NxM**_. Fissata una soglia di soddisfazione, ogni agente si dirà soddisfatto se la percentuale delle celle limitrofe è del suo stesso tipo, altrimenti sarà un agente insoddisfatto e dovrà essere spostato randomicamente in un'altra cella vuota della matrice. Questo avviene ad ogni round fino a che:

- Tutti gli agenti della matrice sono soddisfatti

oppure

- Il numero di round ha oltrepassato il limite massimo di round

___
## Dettagli implementativi

Per la risoluzione del problema è stato scelto di usare **```C```** come linguaggio di programmazione e la libreria **```MPI (Message Passing Interface)```** per lavorare in logica distribuita. Dal seguente diagramma di flusso è possibile evidenziare i passi principali di computazione che hanno permesso di trovare una soluzione al quesito posto. Successivamente verranno illustrati i vari punti nel dettaglio.

![diagramma](./media/diagramma_di_flusso.png)

### _Creazione Matrice_
Per effettuare la generazione della matrice è importante andare a definire il numero preciso di _**X**_, _**O**_ e **'  '** (celle vuote). Questo viene calcolato tramite il prodotto tra righe, colonne e percentuale del valore.

```c
int E = (ROWS*COLUMNS*PERC_E);
int O = (ROWS*COLUMNS-E)*PERC_O, X = (ROWS*COLUMNS-E)*PERC_X;
if(ROWS*COLUMNS != (O+E+X)) E+=(ROWS*COLUMNS)-(O+E+X);
```

Mediante la funzione _**create_matrix**_ è possibile generare una matrice di dimensioni fissate _**ROWSxCOLUMS**_ dove i valori delle celle sono scelti randomicamente tramite la funzione _**rand()**_ il cui seed viene fissato precedentemente in funzione del numero di secondi dell'ora locale. Ogniqualvolta viene inserito un valore, viene decrementato il suo contatore fino a riempire la matrice con le diverse categorie di agenti. 

```c
void create_matrix(char* mat, int E, int O, int X){

    int i, val_rand=0;
    //Seed per funzione di rand()
    time_t now = time(NULL);
    struct tm *tm_struct = localtime(&now);
    srand ( localtime(&now)->tm_sec );

    for(i=0; i<ROWS*COLUMNS; i++){
        val_rand=rand()%3;
        if(val_rand==0 && O>0){
            mat[i] = 'O';
            O--;
        }else if(val_rand==1 && X>0){
            mat[i] = 'X';
            X--;
        }else if(val_rand==2 && E>0){
            mat[i] = ' ';
            E--;
        }else{
            
            if(E>0){
                mat[i] = ' ';
                E--;
            }else if((val_rand==0 || val_rand==2) && (E==0 || O==0)){
                mat[i] = 'X';
                X--;
            }else if((val_rand==1 || val_rand==2) && (E==0 || X==0)){
                mat[i] = 'O';
                O--;
            }
        }
    }
}
```

### _Distribuzione Matrice_
Affinchè si distribuisca la matrice in maniera intelligente è importante andare a definire una struttura dati dedicata che permetta ad ogni processo di salvare la propria sottomatrice, gli indici e il numero di valori della sottomatrice con e senza righe supplementari utili per il calcolo della soddisfazione degli agenti.


```c
typedef struct{
    //Sottomatrice assegnata ad ogni thread 
    char *submatrix;
    
    //Valori per la suddivisione e ricomposizione matrice in funzione di una Scatter
    int *scounts_scatter;
    int *displ_scatter;

    //Valori per la suddivisione e ricomposizione matrice in funzione di una Gather
    int *scounts_gather;
    int *displ_gather;

} Info_submatrix;
```

La distribuzione della matrice di partenza avviene in base al numero di processi coinvolti e non è altro che la sua suddivisione in righe. Al fine di calcolare in un secondo momento la soddisfazione degli agenti, si è scelto di ampliare la sezione di righe dedicata ai singoli thread di una o due righe supplementari, rispettivamente una per il primo (MASTER) e l'ultimo processo, due per i processi a cui sono dedicate sezioni "centrali" della matrice. In merito a ciò, _**scount_scatter**_ e _**displs_scatter**_ rappresentano rispettivamente il numero di elementi e gli indici delle sottomatrici con righe supplementari di ogni thread, _**scount_gather**_ e _**displs_gather**_ invece rappresentano il numero di elementi e gli indici delle sottomatrici effettive quindi senza righe supplementari dedicate ai singoli processi. I primi due array di valori utili per il calcolo della soddisfazione degli agenti, gli altri due utili per la ricomposizione della matrice.

```c
//Distribuzione matrice in sottomatrici per thread: calcolo size e displacements di tutte le sottomatrici dei vari thread
void distribute_matrix(Info_submatrix t_mat, int numproc){
        
        if(numproc>1){
            int i=0;
            int size_r_int = ROWS / numproc;    
            int size_r_rest = ROWS % numproc;

            for(i=0;i<numproc;i++){

                if(i<size_r_rest){
                    t_mat.scounts_scatter[i] = (size_r_int+1) * COLUMNS;
                    t_mat.scounts_gather[i] = (size_r_int+1) * COLUMNS;
                }
                else{
                    t_mat.scounts_scatter[i] = size_r_int * COLUMNS;
                    t_mat.scounts_gather[i] = size_r_int * COLUMNS;
                }

                if(i==0){
                    t_mat.scounts_scatter[i] += COLUMNS;
                    t_mat.displ_scatter[i] = 0;
                }else if(i==(numproc-1)){
                    t_mat.scounts_scatter[i] += COLUMNS;
                    t_mat.displ_scatter[i] = ROWS*COLUMNS - t_mat.scounts_scatter[i];
                }else{
                    t_mat.scounts_scatter[i] += COLUMNS*2;
                    t_mat.displ_scatter[i] = t_mat.displ_scatter[i-1] + t_mat.scounts_scatter[i-1] - COLUMNS*2;
                }
                t_mat.displ_gather[i] = i==0 ? 0 : t_mat.displ_gather[i-1] + t_mat.scounts_gather[i-1];   
            }
        }else{  //Un solo processo
            t_mat.scounts_scatter[0] = ROWS * COLUMNS;
            t_mat.displ_scatter[0] = 0;
            t_mat.scounts_gather[0] = ROWS*COLUMNS;
            t_mat.displ_gather[0] = 0;
        }
}
```

A questo punto, in base ai valori precedentemente calcolati, è possibile suddividere la matrice in sottomatrici e inviare ogni sezione al processo dedicato tramite una _**Scatterv**_

```c
//Suddivisione effettiva della matrice per i vari thread
MPI_Scatterv(mat, t_mat.scounts_scatter, t_mat.displ_scatter, MPI_CHAR, t_mat.submatrix, t_mat.scounts_scatter[myrank], MPI_CHAR, 0, MPI_COMM_WORLD);       
```

### _Calcolo soddisfazione_
Il calcolo della soddisfazione degli agenti viene effettuato da _**satisfaction_step**_ che restituisce sia il numero di celle soddisfatte che una struttura dati dedicata alla memorizzazione di celle insoddisfatte e celle vuote, di seguito presentata.

```c
typedef struct{
    //Array composto dagli indici delle celle insoddisfatte
    char *unsatisfied;  
    //Numero di celle insoddisfatte         
    int n_unsatisfied;
    //Array composto dagli indici delle celle vuote            
    int *freeslots;           
    //Numero di celle vuote   
    int n_freeslots;             

} Info_cellpositions;
```

Per calcolare la soddisfazione di un agente è importante verificare dapprima il rango del processo per capire se si tratta di un thread a cui è stata dedicata una sottomatrice "centrale" oppure no. Una volta verificato ciò, se la cella è piena allora, in base all'indice, verrà chiamata un'apposita funzione che permetterà di calcolare i valori delle celle vicine:
1. per una cella angolo  
    - _**calc_firstangle**_ 
    - _**calc_secondangle**_ 
    - _**calc_thirdangle**_
    - _**calc_fourthangle**_

1. per una cella sulla cornice
    - _**calc_Ledge**_ = _Cornice sinistra_
    - _**calc_Redge**_ = _Cornice destra_
    - _**calc_Aedge**_ = _Cornice in alto_
    - _**calc_Bedge**_ = _Cornice in basso_

1. per una cella centrale
    - _**calc_center**_

Ad ognuna di queste funzioni verrà passata sia la propria sottomatrice che l'indice della cella in esame. Dopo aver calcolato il numero di celle vicine aventi come agente quello della cella passata, il calcolo della soddisfazione sarà dato dal confronto di due variabili: 
```c
return (similar_cells>=similarity);
```
dove
```c
float similarity = PERC_SIM*N; //Grado di soddisfazione 
```

> Il valore **PERC_SIM** è la percentuale di soddisfazione fissata prima della fase di compilazione ([**Note sull'implementazione**](#note_implementazione))
> Il valore **N** varia a seconda della posizione della cella: 
> - 3 per cella angolo 
> - 5 per cella cornice
> - 8 per cella centrale

e 
```c
float similar_cells=0.0; //Numero di celle vicine simili alla cella in esame
```

Se il numero di celle vicine è maggiore o uguale al valore di similarity allora la cella in esame è soddisfatta e verrà restituito 1 come valore. 

___

```c
if(i==0) {
    if(calc_firstangle(t_mat)) local_satisfied++; 
    else{
        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
        cellpos.freeslots[cellpos.n_freeslots]=i;
        cellpos.n_freeslots++;
        cellpos.n_unsatisfied++;
        t_mat.submatrix[i]=' ';
    }
}
```

Prendendo come esempio il calcolo di _**first_angle**_, se il valore restituito è 1 allora viene incrementata la variabile di soddisfatti locali altrimenti:

1. L'indice della cella insoddisfatta viene salvato nell'array delle celle insoddisfatte;
1. L'indice della cella insoddisfatta viene salvato nell'array di celle vuote;
1. Il numero di celle vuote viene incrementato;
1. Il numero di celle insoddisfatte viene incrementato;
1. La cella insoddisfatta viene sostituita con una cella vuota.

Se invece la cella è vuota:
1. L'indice della cella vuota viene salvato nell'array di celle vuote;
1. Il numero di celle vuote viene incrementato.

Nel momento in cui ogni thread ha calcolato il numero delle celle soddisfatte della proprio sottomatrice, tramite una _**MPI_Allreduce**_ ogni processo determinerà il valore della soddisfazione globale, cioè quella dell'intera matrice

```c
MPI_Allreduce(&local_satisfied, &global_satisfied, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
```

### _Riposizionamento insoddisfatti_

Il riposizionamento degli insoddisfatti viene effettuato dalla funzione _displacements_step_. L'obiettivo qui è quello di aggiornare tutti i processi circa l'insieme delle celle insoddisfatte e l'insieme di celle vuote di tutti gli altri thread. Questo permette di tenere traccia sia del numero e degli indici delle celle vuote sia del numero e degli indici delle celle insoddisfatte così da costruire un array di valori che, una volta randomizzato, dia la possibilità di ricollocare gli agenti insoddisfatti in posizioni randomiche dell'intera matrice.

Più precisamente verranno fatte eseguire due _**Allgather**_ rispettivamente per il calcolo del numero di celle insoddisfatte per ogni thread e il calcolo del numero di celle vuote per ogni thread.

```c
MPI_Allgather(&cellpos.n_unsatisfied, 1, MPI_INT, number_of_unsatisfied_per_th, 1, MPI_INT, MPI_COMM_WORLD);

MPI_Allgather(&cellpos.n_freeslots, 1, MPI_INT, number_of_freeslots_per_th, 1, MPI_INT, MPI_COMM_WORLD);
```

Grazie alle informazioni contenute in _**number_of_unsatisfied_per_th**_ e _**number_of_freeslots_per_th**_ è possibile calcolare i displacements delle celle insoddisfatte e quelli delle celle vuote insieme ai quali si vanno a definire due array

- **all_freeslots** = Array di indici di celle vuote dell'intera matrice
- **total_unsat_freeslots** = Array di indici di celle insoddisfatte dell'intera matrice

tramite due _**Allgatherv**_

```c
MPI_Allgatherv(cellpos.freeslots, cellpos.n_freeslots, MPI_INT, all_freeslots, number_of_freeslots_per_th, displs_freeslots, MPI_INT, MPI_COMM_WORLD);

MPI_Allgatherv(cellpos.unsatisfied, cellpos.n_unsatisfied, MPI_CHAR, total_unsat_freeslots, number_of_unsatisfied_per_th, displs_unsatisfied, MPI_CHAR, MPI_COMM_WORLD);
```

Siccome il numero di celle vuote è maggiore o uguale al numero di celle insoddisfatte è importante riempire _**total_unsat_freeslots**_, che al momento contiene solo valori di agenti insoddisfatti, di celle vuote così da poterlo suddividere in sezioni, ognuna dedicata ad un processo. 

```c
//Prima di randomizzare l'array di celle insoddisfatte, viene riempito con ' ' (celle vuote) per far coincidere la size dell'array con il numero di celle vuote totali
for(i=unsatisfied; i<sum_number_freeslots; i++) total_unsat_freeslots[i]=' '; 
```

A questo punto è possibile randomizzare l'array, 

```c
void randomize ( char *arr, int n )
{
    for (int i = n-1; i > 0; i--)
    {
        int j = rand() % (i+1);
        swap(&arr[i], &arr[j]);
    }
}
```

sostituire i valori della sottomatrice con quelli presenti nell'array randomizzato

```c
//Vengono sostituiti i valori nei freeslots delle varie sottomatrici con i valori delle sezioni di total_unsatisfied in funzione del proprio rank 
for(i=0;i<cellpos.n_freeslots;i++) t_mat.submatrix[cellpos.freeslots[i]]=total_unsat_freeslots[displs_freeslots[myrank]+i];
```

e aggiornare le righe supplementari con i nuovi valori ricollocati nella loro nuova posizione all'interno della matrice. 

```c
//Prima di assegnare i nuovi valori delle celle insoddisfatti alle matrici dei vari thread, vengono aggiornate le righe supplementari utili al calcolo della soddisfazione al successivo round

if(myrank==0)   //Rank==0 -> MASTER quindi bisogna aggiornare solamente una riga supplementare cioè l'ultima della sottomatrice
    for(i=0; i<number_of_freeslots_per_th[1];i++){
        if(all_freeslots[displs_freeslots[1]+i]>=(COLUMNS*2)) break;
        else if(all_freeslots[displs_freeslots[1]+i]>=COLUMNS){
            t_mat.submatrix[t_mat.scounts_scatter[myrank]-COLUMNS+(all_freeslots[displs_freeslots[1]+i]%COLUMNS)] = total_unsat_freeslots[displs_freeslots[1]+i];
        }
    }
else if(myrank==(numproc-1)){   //Rank==(totale processi - 1) -> ULTIMO THREAD quindi bisogna aggiornare solamente una riga cioè la prima della sottomatrice
    for(i=0; i<number_of_freeslots_per_th[numproc-2];i++)
        if(all_freeslots[displs_freeslots[numproc-2]+i]>=(t_mat.scounts_scatter[numproc-2]-(COLUMNS*2))){
            t_mat.submatrix[all_freeslots[(displs_freeslots[numproc-2]+i)]%COLUMNS] = total_unsat_freeslots[displs_freeslots[numproc-2]+i];
        }
}
else{   //Rank compreso tra 0 e l'ultimo thread quindi bisogna aggiornare due righe supplementari, la prima e l'ultima della sottomatrice
    //PRIMA RIGA
    for(i=0; i<number_of_freeslots_per_th[myrank-1];i++)
        if(all_freeslots[displs_freeslots[myrank-1]+i]>=(t_mat.scounts_scatter[myrank-1]-(COLUMNS*2)) &&  all_freeslots[displs_freeslots[myrank-1]+i]<(t_mat.scounts_scatter[myrank-1]-COLUMNS)){ 
            t_mat.submatrix[all_freeslots[displs_freeslots[myrank-1]+i]%COLUMNS] = total_unsat_freeslots[displs_freeslots[myrank-1]+i];
        }
    //ULTIMA RIGA
    for(i=0; i<number_of_freeslots_per_th[myrank+1];i++)
        if(all_freeslots[displs_freeslots[myrank+1]+i]>=COLUMNS && all_freeslots[displs_freeslots[myrank+1]+i]<(COLUMNS*2)){
            t_mat.submatrix[t_mat.scounts_scatter[myrank]-COLUMNS+(all_freeslots[displs_freeslots[myrank+1]+i]%COLUMNS)] = total_unsat_freeslots[displs_freeslots[myrank+1]+i];
        }
}
```

### _Ricomposizione Matrice_
Se tutti gli agenti sono soddisfatti o è stato superato il numero di round massimo, allora si esce dal ciclo, viene ricomposta la matrice finale e vengono stampati i risultati.

La ricomposizione della matrice avviene tramite la funzione _**recompose_mat**_ che non fa altro che generare un array con i valori della sottomatrice assegnata al processo e inviarlo tramite una _**Gatherv**_ al master.

```c
MPI_Gatherv(submatGather, t_mat.scounts_gather[myrank], MPI_CHAR, final_mat, t_mat.scounts_gather, t_mat.displ_gather, MPI_CHAR, 0, MPI_COMM_WORLD);
```
### _Stampa dei risultati_
La stampa dei risultati è l'ultimo passo del flusso di operazioni del programma. I dati stampati a video riguardano:

- Matrice iniziale e matrice finale
- Risultato ottenuto dalla computazione
    - Numero di celle soddisfatte
    - Percentuale di soddisfazione ottenuta
    - Tempo di computazione 
    - Numero di round raggiunto
- Parametri in input  
    - Taglia matrice 
    - Percentuale valori X, O e celle vuote
    - Percentuale di soddisfazione

___
## Note sull'implementazione 
Per rendere l'implementazione grafica più accattivante è stato pensato di dare la possibilità di visualizzare la matrice aggiornata ad ogni round, ovviamente a discapito delle prestazioni. Per questo motivo si consiglia di optare per questa scelta qualora non si vogliano misurare le performance del programma.

Per fare ciò è fondamentale la modifica dei valori **DELAY**, **PRINT_ROUNDS** e **PERFORMANCE**. Il primo per sospendere il processo master per secondi dopo la stampa della matrice aggiornata, il secondo valore per consentire la visualizzazione a video della matrice risultante mentre il terzo per autorizzare la stampa a sfavore delle prestazioni. 
```c
//Delay di attesa per round
#define DELAY 0

//'1' se si vuole stampare la matrice risultante ad ogni round
#define PRINT_ROUNDS 1          

//'1' se si vuole calcolare matrice risultante senza stampe a favore delle prestazioni
#define PERFORMANCE 0           
```

Inoltre, prima di effettuare la compilazione, è importante settare i valori utili per la costruzione della matrice, per il calcolo della soddisfazione e per il numero di round di computazione massimi.

```c
//Numero Colonne
#define COLUMNS 30

//Numero Righe
#define ROWS 30

//Numero Round Massimi
#define N_ROUND_MAX 300         

//Percentuale di 'O'
#define PERC_O 0.5              

//Percentuale di 'X'
#define PERC_X (1-PERC_O)       

//Percentuale di celle vuote
#define PERC_E 0.3              

//Percentuale di soddisfazione
#define PERC_SIM 0.3            
```

### _Compilazione_

### _Esecuzione_

___
## Benchmarks
bench
___

## Correttezza
Correttezza: dimostrare che due esecuzioni con 4 processori e 4 processori con lo stesso seed per la funzione di random nella costruzione della matrice, generano lo stesso output. 
___

## Conclusioni
Conclusioni