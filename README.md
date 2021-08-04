# [Schelling's-Model-of-Segregation](http://nifty.stanford.edu/2014/mccown-schelling-model-segregation/)

| **Petrazzuolo Lorenzo** | **0522500894** | **23/08/2021** |
| --- | --- | --- |

Progetto per l'esame di _Programmazione Concorrente, Parallela e sul Cloud_ dell'anno di corso _2020/2021_. 
Laurea magistrale in Computer Science, curriculum in Cloud Computing.

## Introduzione 
___
### _Il modello di segregazione di Shelling: Cos'è_

Nel 1971, l'economista americano Thomas Schelling creò un modello basato su agenti che suggeriva che anche il comportamento involontario potesse contribuire alla segregazione. Il suo modello di segregazione ha mostrato che anche quando gli individui (o "agenti") non si preoccupavano di essere circondati o vivere da agenti di una razza o un background economico diverso, avrebbero comunque scelto di separarsi dagli altri agenti nel tempo. Sebbene il modello sia abbastanza semplice, fornisce uno sguardo affascinante su come gli individui potrebbero auto-segregarsi, anche quando non hanno un desiderio esplicito di farlo.
> _[Riferimenti](https://en.wikipedia.org/wiki/Schelling's_model_of_segregation)_


### _Descrizione del problema_

Nel seguente progetto si propone un'implementazione di una simulazione del modello di Shelling. Una volta impostati una serie di parametri, il programma dovrà essere in grado di descrivere graficamente le caratteristiche del modello. Più precisamente si suppone di avere due tipi di agenti _**X**_ e _**O**_ che andranno a comporre la popolazione di una matrice _**NxM**_. Fissata una soglia di soddisfazione, ogni agente si dirà soddisfatto se la percentuale delle celle limitrofe è del suo stesso tipo, altrimenti sarà un agente insoddisfatto e dovrà essere spostato randomicamente in un'altra cella vuota della matrice. Questo avviene ad ogni round fino a che:

- Tutti gli agenti della matrice sono soddisfatti

oppure

- Il numero di round ha oltrepassato il limite massimo di round

## Dettagli implementativi
___

Per la risoluzione del problema è stato scelto di usare **C** come linguaggio di programmazione e la libreria **MPI** per lavorare in logica distribuita. Dal seguente diagramma di flusso è possibile evidenziare i quattro passi principali di computazione che hanno permesso di trovare una soluzione al quesito posto. Successivamente verranno illustrati i vari punti nel dettaglio.

![diagramma](./media/diagramma_di_flusso.png)








## Note sull'implementazione 
___





### Compilazione
### Esecuzione

## Benchmarks
___

## Correttezza
___

## Conclusioni
___