#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#include "mpi.h"
#include <time.h>

#define COLUMNS 25              //Numero Colonne              
#define ROWS 35                 //Numero Righe
#define N_ROUND_MAX 300         //Numero Round Massimi
#define PERC_O 0.5              //Percentuale di 'O'
#define PERC_X (1-PERC_O)       //Percentuale di 'X'
#define PERC_E 0.3              //Percentuale di celle vuote
#define PERC_SIM 0.10           //Percentuale di soddisfazione
#define DELAY 1                 //Delay di attesa per round
#define PRINT_ROUNDS 1          //'1' se si vuole stampare la matrice risultante ad ogni round
#define PERFORMANCE 0           // 1-> Calcolo matrice risultante senza stampe a favore delle prestazioni

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

typedef struct{
    
    char *unsatisfied;           //Array composto dagli indici delle celle insoddisfatte
    int n_unsatisfied;           //Numero di celle insoddisfatte
    int *freeslots;              //Array composto dagli indici delle celle vuote
    int n_freeslots;             //Numero di celle vuote

} Info_cellpositions; 

_Bool calc_firstangle(Info_submatrix mat);
_Bool calc_secondangle(Info_submatrix mat);
_Bool calc_thirdangle(Info_submatrix mat, int);
_Bool calc_fourthangle(Info_submatrix mat, int);
_Bool calc_Ledge(Info_submatrix mat, int, int myrank);
_Bool calc_Redge(Info_submatrix mat, int, int myrank);
_Bool calc_Aedge(Info_submatrix mat, int);
_Bool calc_Bedge(Info_submatrix mat, int, int myrank);
_Bool calc_center(Info_submatrix mat, int, int myrank); 
int sum(int *, int);
void swap(char *, char *);
void randomize(char *, int);
void create_matrix(char* mat, int E, int O, int X);
void print_matrix(char *i_mat);
void distribute_matrix(Info_submatrix t_mat, int numproc);
Info_cellpositions satisfaction_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions, int *);
void displacements_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions, int, int);
void recompose_mat(int myrank, Info_submatrix t_mat, char *final_mat, int numproc);

int main(int argc, char *argv[]){
	system("clear");
    
	int numproc, myrank, source, *arr, i, value, round=1, satisfied=0;
    double start, end;
    char *mat, *final_mat;
    Info_submatrix t_mat;
    Info_cellpositions cellpos;

    //Definizione seed per creazione randomica della matrice
    time_t now = time(NULL);
    struct tm *tm_struct = localtime(&now);
    int hour = tm_struct->tm_hour;

    //Calcolo numero di valori 'O', 'X' e ' ' per costruzione matrice
    int E = (ROWS*COLUMNS*PERC_E);
    int O = (ROWS*COLUMNS-E)*PERC_O, X = (ROWS*COLUMNS-E)*PERC_X;
    if(ROWS*COLUMNS != (O+E+X)) E+=(ROWS*COLUMNS)-(O+E+X);

    //Inizializzazione MPI
	MPI_Init(&argc, &argv);

	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    //Definizione valori di Info_subamtrix
    t_mat.scounts_scatter = (int *)malloc(sizeof(int)*numproc);
    t_mat.scounts_gather = (int *)malloc(sizeof(int)*numproc);
    t_mat.displ_scatter = (int *)malloc(sizeof(int)*numproc);
    t_mat.displ_gather = (int *)malloc(sizeof(int)*numproc);

    if(numproc<=ROWS){
        if(myrank == 0){
            int val=1;
            mat = (char *)malloc(ROWS * COLUMNS * sizeof(char));
            final_mat = (char *)malloc(ROWS * COLUMNS * sizeof(char));
            //Creazione matrice con valori randomici
            create_matrix(mat, E, O, X);    
            if(!PERFORMANCE) print_matrix(mat);
            sleep(DELAY);
            printf("\n");
        }

        //Distribuzione della matrice per il numero di thread
        distribute_matrix(t_mat, numproc);
        t_mat.submatrix = (char *)malloc(sizeof(char)*t_mat.scounts_scatter[myrank]);

        //Suddivisione effettiva della matrice per i vari thread
        MPI_Scatterv(mat, t_mat.scounts_scatter, t_mat.displ_scatter, MPI_CHAR, t_mat.submatrix, t_mat.scounts_scatter[myrank], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        start = MPI_Wtime();

        while(round<=N_ROUND_MAX){

            cellpos.unsatisfied = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]); 
            cellpos.freeslots = (int *)malloc(sizeof(int)*t_mat.scounts_gather[myrank]);
            
            //1. Calcolo insoddisfatti la cui posizione viene sostituita con ' ';
            //2. Calcolo freeslots;
            //3. Calcolo soddisfazione totale
            cellpos = satisfaction_step(t_mat, myrank, numproc, cellpos, &satisfied);

            //Se tutte le celle con valore 'X' e 'O' sono soddisfatte allora si esce dal while
            if(satisfied==(ROWS*COLUMNS-E)){
                free(cellpos.freeslots);
                free(cellpos.unsatisfied);
                break;
            }
           
            //Riempio array di freeslots con valori di insoddisfatti e randomizzo l'array. 
            //Divido l'array in sezioni, ognuna corrispondente a quella del rank. 
            //Infine sostituisco i valori dell'array con i valori dei freeslots di ogni processo.
            displacements_step(t_mat, myrank, numproc, cellpos, (ROWS*COLUMNS-E-satisfied), round);
            
            if(PRINT_ROUNDS){
                //printf("\033[%d;%dH",0,0);
                printf("\e[1;1H\e[2J\n");
                if(numproc>1){
                    //Ricomposizione e stampa matrice ad ogni round se PRINT_ROUNDS==1
                    recompose_mat(myrank, t_mat, final_mat, numproc);
                    if(myrank==0) print_matrix(final_mat);
                    free(final_mat);
                    final_mat = (char *)malloc(sizeof(char)*ROWS*COLUMNS);
                }
                else{
                    //Caso numero di processi è uguale a 1
                    print_matrix(t_mat.submatrix);
                }
            }

            //Stampa aggiornamento soddisfatti per ogni round
            if(myrank == 0){
                if(!PERFORMANCE) printf("\nRound: %d/%d\nSatisfied: %.1f %% (%d/%d)\n", round, N_ROUND_MAX,((float)satisfied*100)/(float)(ROWS*COLUMNS-E), satisfied, ROWS*COLUMNS-E);
                sleep(DELAY);
            }
            
            free(cellpos.freeslots);
            free(cellpos.unsatisfied);
            round++;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        system("clear");
        
        end = MPI_Wtime();

        if(numproc>1) recompose_mat(myrank, t_mat, final_mat, numproc);
        else final_mat=t_mat.submatrix;

        //Stampa dei risultati 
        if(myrank==0){
            if(!PERFORMANCE){
                printf("\n\tSTARTING MATRIX\n");
                print_matrix(mat);
                printf("\n\tFINAL MATRIX\n");
                print_matrix(final_mat);
            }
            if(satisfied==(ROWS*COLUMNS-E))
                printf("\n\n \033[32m ALL AGENTS WERE SATISFIED !!!\x1b[0m \n\n");
            else
                printf("\n\n \x1b[31m  NOT ALL AGENTS ARE SATISFIED..\x1b[0m \n\n");
            printf("\n\n\tMATRIX SIZE: %dx%d \n\n\tSIMILAR: %d%%\n\n\tRED/BLUE: %d%%/%d%%\n\n\tEMPTY: %d%% \n\n\tTIME PASSED: %f \n\n\tSATISFIED: %d/%d\n\n\tPERC SATISFIED: %.1f %%\n\n\tROUNDS: %d/%d\n\n\n", ROWS, COLUMNS, (int)(PERC_SIM*100),(int)(PERC_O*100), (int)(PERC_X*100), (int)(PERC_E*100),end-start, satisfied, (ROWS*COLUMNS-E),((float)satisfied*100)/(float)(ROWS*COLUMNS-E), round, N_ROUND_MAX);
        }
	}else{
        if(myrank==0) printf("\n\n\n\x1b[31mPlease enter fewer processes than %d (ROWS NUMBER) \x1b[0m\n\n\n\n\n", ROWS);
    }

	MPI_Finalize();

    free(t_mat.displ_gather);
    free(t_mat.displ_scatter);
    free(t_mat.scounts_gather);
    free(t_mat.scounts_scatter);
    free(t_mat.submatrix);
    free(mat);
    free(final_mat);

	return 0;
}

//Ricomposizione matrice finale: ogni thread, tramite una Gather, comunica la propria sottomatrice 
//al master il quale costruirà la matrice finale
void recompose_mat(int myrank, Info_submatrix t_mat, char *final_mat, int numproc){

        int i=0, j=0;
        char *submatGather = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]);
        if(numproc>1){

            if(myrank==0) for(i=0; i<t_mat.scounts_gather[myrank];i++) submatGather[i]=t_mat.submatrix[i];
            else if(myrank==(numproc-1)){
                
                for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
                    submatGather[j]=t_mat.submatrix[i];
                    j++;
                }
            }
            else{
                for(i=COLUMNS; i<(t_mat.scounts_scatter[myrank]-COLUMNS);i++){
                    submatGather[j]=t_mat.submatrix[i];
                    j++;
                }
            }
            MPI_Gatherv(submatGather, t_mat.scounts_gather[myrank], MPI_CHAR, final_mat, t_mat.scounts_gather, t_mat.displ_gather, MPI_CHAR, 0, MPI_COMM_WORLD);
            
        }
        free(submatGather);
}

//Somma valori contenuti in un array
int sum(int *arr, int length){
    int sum=0;
    for(int i=0;i<length;i++) sum+=arr[i];
    return sum;
}

//Funzione dedicata allo spostamento delle celle insoddisfatte e all'aggiornamento delle righe supplementari di ogni sottomatrice di ogni thread
void displacements_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int unsatisfied, int round){
    
    int sum_number_freeslots=0,                 //Somma totale del numero di celle vuote per thread
        displs_unsatisfied[numproc],            //Displacements degli insoddisfatti per thread
        displs_freeslots[numproc],              //Displacements delle celle vuote per thread
        number_of_unsatisfied_per_th[numproc],  //Numero di insoddisfatti per thread
        number_of_freeslots_per_th[numproc],    //Numero di celle vuote per thread
        i=0;
    
    if(numproc>1){
        
        MPI_Allgather(&cellpos.n_unsatisfied, 1, MPI_INT, number_of_unsatisfied_per_th, 1, MPI_INT, MPI_COMM_WORLD);

        MPI_Allgather(&cellpos.n_freeslots, 1, MPI_INT, number_of_freeslots_per_th, 1, MPI_INT, MPI_COMM_WORLD);

        sum_number_freeslots = sum(number_of_freeslots_per_th, numproc);

        int all_freeslots[sum_number_freeslots];            //Array di indici di tutte le celle vuote dell'intera matrice
        char total_unsat_freeslots[sum_number_freeslots];   //Array di indici di celle insoddisfatte dell'intera matrice da randomizzare
        
        for(i=0; i<numproc; i++){
            if(i==0){
                displs_unsatisfied[i]=0;
                displs_freeslots[i]=0;
            }
            else{
                displs_unsatisfied[i]=displs_unsatisfied[i-1]+number_of_unsatisfied_per_th[i-1];
                displs_freeslots[i]=displs_freeslots[i-1]+number_of_freeslots_per_th[i-1];
            }
        }

        MPI_Allgatherv(cellpos.freeslots, cellpos.n_freeslots, MPI_INT, all_freeslots, number_of_freeslots_per_th, displs_freeslots, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(cellpos.unsatisfied, cellpos.n_unsatisfied, MPI_CHAR, total_unsat_freeslots, number_of_unsatisfied_per_th, displs_unsatisfied, MPI_CHAR, MPI_COMM_WORLD);

        //Prima di randomizzare l'array di celle insoddisfatte, viene riempito con ' ' (celle vuote) per far coincidere la size dell'array 
        //con il numero di celle vuote totali
        for(i=unsatisfied; i<sum_number_freeslots; i++) total_unsat_freeslots[i]=' ';
        
        srand ( round );    //Seed per la funzione rand()
        
        randomize(total_unsat_freeslots, sum_number_freeslots);     //Randomizzazione array di insoddisfatti + celle vuote rimanenti 

        //Prima di assegnare i nuovi valori delle celle insoddisfatti alle matrici dei vari thread,
        //vengono aggiornate le righe supplementari utili al calcolo della soddisfazione al successivo round
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

        //Vengono sostituiti i valori nei freeslots delle varie sottomatrici con i valori delle sezioni di total_unsatisfied in funzione del proprio rank 
        for(i=0;i<cellpos.n_freeslots;i++) t_mat.submatrix[cellpos.freeslots[i]]=total_unsat_freeslots[displs_freeslots[myrank]+i];
    }
    else{
        //UN SOLO PROCESSO
        char total_unsat_freeslots[cellpos.n_freeslots];
        for(i=0;i<cellpos.n_freeslots;i++) 
            if(i<cellpos.n_unsatisfied) total_unsat_freeslots[i]=cellpos.unsatisfied[i];
            else total_unsat_freeslots[i]=' ';

        randomize(total_unsat_freeslots, cellpos.n_freeslots);

        for(i=0;i<cellpos.n_freeslots;i++)
            t_mat.submatrix[cellpos.freeslots[i]]=total_unsat_freeslots[i];

    }
}

void swap (char *a, char *b)
{
    char temp = *a;
    *a = *b;
    *b = temp;
}

//Funzione di randomizzazione 
void randomize ( char *arr, int n )
{
    
    for (int i = n-1; i > 0; i--)
    {
        int j = rand() % (i+1);
 
        swap(&arr[i], &arr[j]);
    }
}

//Identificazione di celle insoddisfatte e celle vuote. Ogni cella insoddisfatta verrà sostituita con una cella vuota
Info_cellpositions satisfaction_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int *satisfied){

    int i=0, local_satisfied=0, global_satisfied=0;
    cellpos.n_unsatisfied=0;
    cellpos.n_freeslots=0;

    if(numproc>1){
        if(myrank==0){

            for(i=0; i<(t_mat.scounts_scatter[myrank]-COLUMNS);i++){
                if(t_mat.submatrix[i]!=' '){
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
                    else if(i==(COLUMNS-1)){
                        if(calc_secondangle(t_mat)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i%COLUMNS==0){
                        if(calc_Ledge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i%COLUMNS==(COLUMNS-1)){
                        if(calc_Redge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i>0 && i<COLUMNS){
                        if(calc_Aedge(t_mat, i)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else{
                        if(calc_center(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                }
                else{
                    cellpos.freeslots[cellpos.n_freeslots]=i;
                    cellpos.n_freeslots++;
                }
            }

        }else if(myrank==(numproc-1)){

            for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
                if(t_mat.submatrix[i]!=' '){
                    if(i%COLUMNS==0 && i!=(t_mat.scounts_scatter[myrank]-COLUMNS+1)){
                        if(calc_Ledge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i%COLUMNS==(COLUMNS-1) && i!=(t_mat.scounts_scatter[myrank]-1)){
                        if(calc_Redge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i>(t_mat.scounts_scatter[myrank]-COLUMNS+1) && i<(t_mat.scounts_scatter[myrank]-1)){
                        if(calc_Bedge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i==(t_mat.scounts_scatter[myrank]-COLUMNS+1)){
                        if(calc_thirdangle(t_mat, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i==(t_mat.scounts_scatter[myrank]-1)){
                        if(calc_fourthangle(t_mat, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else{
                        if(calc_center(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                }
                else{
                    cellpos.freeslots[cellpos.n_freeslots]=i;
                    cellpos.n_freeslots++;
                }
            }

        }else{
            
            for(i=COLUMNS; i<(t_mat.scounts_scatter[myrank]-COLUMNS);i++){
                if(t_mat.submatrix[i]!=' '){
                    if(i%COLUMNS==0){
                        if(calc_Ledge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else if(i%COLUMNS==(COLUMNS-1)){
                        if(calc_Redge(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                    else{ 
                        if(calc_center(t_mat, i, myrank)) local_satisfied++; 
                        else{
                            cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                            cellpos.freeslots[cellpos.n_freeslots]=i;
                            cellpos.n_freeslots++;
                            cellpos.n_unsatisfied++;
                            t_mat.submatrix[i]=' ';
                        }
                    }
                }
                else{
                    cellpos.freeslots[cellpos.n_freeslots]=i;
                    cellpos.n_freeslots++;
                }
            }

        }
        //Calcolo soddisfazione totale
        MPI_Allreduce(&local_satisfied, &global_satisfied, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        *satisfied = global_satisfied;

    }
    else{

        // UN PROCESSO
        for(i=0; i<t_mat.scounts_scatter[myrank]; i++){
            if(t_mat.submatrix[i]!=' '){
                if(i==0) 
                    if(calc_firstangle(t_mat)) local_satisfied++;
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                else if(i==(COLUMNS-1)){
                    if(calc_secondangle(t_mat)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i>0 && i<COLUMNS){
                    if(calc_Aedge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i%COLUMNS==0 && i!=(t_mat.scounts_scatter[myrank]-COLUMNS+1)){
                    if(calc_Ledge(t_mat, i, myrank)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i%COLUMNS==(COLUMNS-1) && i!=(t_mat.scounts_scatter[myrank]-1)){
                    if(calc_Redge(t_mat, i, myrank)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i==(t_mat.scounts_scatter[myrank]-COLUMNS+1)){
                    if(calc_thirdangle(t_mat, myrank)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i==(t_mat.scounts_scatter[myrank]-1)){
                    if(calc_fourthangle(t_mat, myrank)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else{
                    if(calc_center(t_mat, i, myrank)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
            }else{
                cellpos.freeslots[cellpos.n_freeslots]=i;
                cellpos.n_freeslots++;
            }
        }
        *satisfied = local_satisfied;
    }
    

    return cellpos;
}

//Calcolo soddisfazione cella angolo 
_Bool calc_firstangle(Info_submatrix mat){
    float similarity = PERC_SIM*3;
    float similar_cells=0.0;
    if(mat.submatrix[0]==mat.submatrix[1]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS+1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella angolo 
_Bool calc_secondangle(Info_submatrix mat){
    float similarity = PERC_SIM*3;
    float similar_cells=0.0;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[COLUMNS-2]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-1]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-2]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella angolo 
_Bool calc_thirdangle(Info_submatrix mat, int myrank){
    float similarity = PERC_SIM*3;
    float similar_cells=0.0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+1]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+2]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella angolo 
_Bool calc_fourthangle(Info_submatrix mat, int myrank){
    float similarity = PERC_SIM*3;
    float similar_cells=0.0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-1-COLUMNS]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2-COLUMNS]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella cornice sinistra 
_Bool calc_Ledge(Info_submatrix mat, int index, int myrank){
    float similarity = PERC_SIM*5;
    float similar_cells=0.0;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella cornice destra
_Bool calc_Redge(Info_submatrix mat, int index, int myrank){
    float similarity = PERC_SIM*5;
    float similar_cells=0.0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella cornice alta
_Bool calc_Aedge(Info_submatrix mat, int index){
    float similarity = PERC_SIM*5;
    float similar_cells=0.0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella cornice bassa
_Bool calc_Bedge(Info_submatrix mat, int index, int myrank){
    float similarity = PERC_SIM*5;
    float similar_cells=0.0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Calcolo soddisfazione cella centrale
_Bool calc_center(Info_submatrix mat, int index, int myrank){
    float similarity = PERC_SIM*8;
    float similar_cells=0.0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    return (similar_cells>=similarity);
}

//Creazione matrice: il numero di valori da inserire nella matrice sono già precalcolati nella main function
//e vengono decrementati ogniqualvolta ne viene inserito uno nuovo
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

//Stampa matrice
void print_matrix(char *mat) {

    int i=0;
    printf("\n");
    for (i=0; i<(ROWS*COLUMNS); i++) {
        printf("|");
        switch(mat[i])
        {
        case 'O':
            printf("\x1b[31m %c \x1b[0m", mat[i]);
            break;
        case 'X':
            printf("\033[0;34m %c \x1b[0m", mat[i]);
            break;
        case ' ':
            printf("   ");
            break;
        }
        
        if (((i+1)%(COLUMNS)==0) && (i!=0))
            printf("|\n");
        if ((ROWS*COLUMNS)==1)
            printf("|\n");
        if (COLUMNS==1 && i==0)
            printf("|\n");
    }
    printf("\n");
}

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