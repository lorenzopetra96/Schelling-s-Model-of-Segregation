#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#include "mpi.h"
#include <time.h>


#define COLUMNS 10
#define ROWS 10
#define N_ROUND_MAX 4
#define PERC_O 0.5
#define PERC_X (1-PERC_O)
#define PERC_E 0.2
#define PERC_SIM 0.3
#define DELAY 1
#define PRINT_ROUNDS 0
#define SEED 

typedef struct{

    char *submatrix;

    int *scounts_scatter;
    int *displ_scatter;

    int *scounts_gather;
    int *displ_gather;


} Info_submatrix;

typedef struct{

    char *unsatisfied; 
    int n_unsatisfied;
    int *freeslots;
    int n_freeslots;

} Info_cellpositions; 


_Bool calc_firstangle(Info_submatrix mat);
_Bool calc_secondangle(Info_submatrix mat);
_Bool calc_thirdangle(Info_submatrix mat, int);
_Bool calc_fourthangle(Info_submatrix mat, int);
_Bool calc_Ledge(Info_submatrix mat, int);
_Bool calc_Redge(Info_submatrix mat, int);
_Bool calc_Aedge(Info_submatrix mat, int);
_Bool calc_Bedge(Info_submatrix mat, int);
_Bool calc_center(Info_submatrix mat, int); 
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

    time_t now = time(NULL);
    struct tm *tm_struct = localtime(&now);
    int hour = tm_struct->tm_hour;
    
    srand ( 1 );
    
    //Per il controllo della satisfaction
    int E = (ROWS*COLUMNS*PERC_E);
    int O = (ROWS*COLUMNS-E)*PERC_O, X = (ROWS*COLUMNS-E)*PERC_X;
    if(ROWS*COLUMNS != (O+E+X)) E+=(ROWS*COLUMNS)-(O+E+X);
    

	//srand(time(NULL)+myrank);

	MPI_Init(&argc, &argv);

	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    t_mat.scounts_scatter = (int *)malloc(sizeof(int)*numproc);
    t_mat.scounts_gather = (int *)malloc(sizeof(int)*numproc);
    t_mat.displ_scatter = (int *)malloc(sizeof(int)*numproc);
    t_mat.displ_gather = (int *)malloc(sizeof(int)*numproc);
    

    if(numproc<=ROWS){
        if(myrank == 0){
            int val=1;
            mat = (char *)malloc(ROWS * COLUMNS * sizeof(char));
            final_mat = (char *)malloc(ROWS * COLUMNS * sizeof(char));
            create_matrix(mat, E, O, X);
            printf("MATRICE INIZIALE\n"); val++;
            print_matrix(mat);
            sleep(2);
            printf("\n");
        }

        distribute_matrix(t_mat, numproc);
        t_mat.submatrix = (char *)malloc(sizeof(char)*t_mat.scounts_scatter[myrank]);

        
        
        MPI_Scatterv(mat, t_mat.scounts_scatter, t_mat.displ_scatter, MPI_CHAR, t_mat.submatrix, t_mat.scounts_scatter[myrank], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while(round<=N_ROUND_MAX){

            cellpos.unsatisfied = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]); 
            cellpos.freeslots = (int *)malloc(sizeof(int)*t_mat.scounts_gather[myrank]);

            
            //Calcolo insoddisfatti la cui posizione viene sostituita con ' '. Calcolo freeslots.
            cellpos = satisfaction_step(t_mat, myrank, numproc, cellpos, &satisfied);

            

            if(satisfied==(ROWS*COLUMNS-E)){
                free(cellpos.freeslots);
                free(cellpos.unsatisfied);
                break;
            }
            //if(myrank == 0) printf("\nNEL WHILE-------RANK: %d - FREESLOTS: %d\nSATISFIED: %d\n", myrank, cellpos.n_freeslots, satisfied);
            
            //if(myrank==2) for(i=0;i<3;i++) printf("\nELEM: %c", cellpos.unsatisfied[i]);
            
            //Riempio array di freeslots con valori di insoddisfatti e randomizzo l'array. Divido l'array
            //in sezioni, ognuna corrispondente a quella del rank. Infine sostituisco i valori dell'array
            //con i valori dei freeslots di ogni processo 
            displacements_step(t_mat, myrank, numproc, cellpos, (ROWS*COLUMNS-E-satisfied), round);
            
            if(PRINT_ROUNDS){
                //RICOMPONI E STAMPA MATRICE
                recompose_mat(myrank, t_mat, final_mat, numproc);
                print_matrix(final_mat);
            }

            if(myrank == 0){
                
                
                printf("\nNumero round: %d\nSoddisfatti: %1f %%\n", round, ((float)satisfied*100)/(float)(ROWS*COLUMNS-E));
                sleep(2);
                
            }
            
            free(cellpos.freeslots);
            free(cellpos.unsatisfied);
            round++;
        }


        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();

        if(myrank == 0 && satisfied==(ROWS*COLUMNS-E)){
            printf("\n\n\n\n___________________USCITO PER SODDISFAZIONE___________________\n\n\n\n\n\n");
        }

        recompose_mat(myrank, t_mat, final_mat, numproc);
        if(myrank==0){
            printf("\n\tMATRICE INIZIALE\n");
            print_matrix(mat);
            printf("\n\tMATRICE FINALE\n");
            print_matrix(final_mat);
            printf("\n\n\tTAGLIA MATRICE: %dx%d \n\n\tRED/BLUE: %d%%/%d%%\n\n\tEMPTY: %d%% \n\n\tTEMPO IMPIEGATO: %f \n\n\tNUMERO SODDISFATTI: %d/%d\n\n\tPERC SODDISFATTI: %.1f %%\n\n\n\n", ROWS, COLUMNS, (int)PERC_O*100, (int)PERC_X*100, (int)PERC_E*100,end-start, satisfied, (ROWS*COLUMNS-E),((float)satisfied*100)/(float)(ROWS*COLUMNS-E));
        }
	}else{
        //aggiungere controllo se numproc>ROWS
        printf("numproc>ROWS... Ripetere");
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

void recompose_mat(int myrank, Info_submatrix t_mat, char *final_mat, int numproc){

        char *submatGather = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]);
        int i=0, j=0;
        
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
        free(submatGather);
}

int sum(int *arr, int length){
    int sum=0;
    for(int i=0;i<length;i++) sum+=arr[i];
    return sum;
}

void displacements_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int unsatisfied, int round){

    
    int sum_number_freeslots=0, displs_unsatisfied[numproc], displs_freeslots[numproc], number_of_unsatisfied_per_th[numproc], number_of_freeslots_per_th[numproc], i=0;
    
    MPI_Allgather(&cellpos.n_unsatisfied, 1, MPI_INT, number_of_unsatisfied_per_th, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Allgather(&cellpos.n_freeslots, 1, MPI_INT, number_of_freeslots_per_th, 1, MPI_INT, MPI_COMM_WORLD);

    //if(myrank==0) printf("\nSASAJDOASJODAS: %d ", unsatisfied);

    //MPI_Allreduce(&cellpos.n_freeslots, &sum_number_freeslots, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    sum_number_freeslots = sum(number_of_freeslots_per_th, numproc);

    //if(myrank==0) printf("\nSOMMA SLOTS: %d\nUNSATISF: %d\n", sum_number_freeslots, unsatisfied);
    
    char total_unsat_freeslots[sum_number_freeslots];
    
    for(i=0; i<numproc; i++){
        if(i==0){
            displs_unsatisfied[i]=0;
            displs_freeslots[i]=0;
        }
        else{
            displs_unsatisfied[i]=displs_unsatisfied[i-1]+number_of_unsatisfied_per_th[i-1];
            displs_freeslots[i]=displs_freeslots[i-1]+number_of_freeslots_per_th[i-1];
            //if(myrank==0) printf("\n%d) N_INSODD: %d\nN_FREESLOTS %d\n", i-1, number_of_unsatisfied_per_th[i-1], number_of_freeslots_per_th[i-1]);
        }
    }

    //if(myrank==1) for(i=0;i<3;i++) printf("\nN_UNSAT: %d - N_FREES: %d - DISPLS_UNSAT: %d - DISPLS_FREES: %d\n", number_of_unsatisfied_per_th[i], number_of_freeslots_per_th[i], displs_unsatisfied[i], displs_freeslots[i]);
     

    //printf("\nRANK: %d - N_UNSAT: %d\n", myrank, cellpos.n_unsatisfied);

    MPI_Allgatherv(cellpos.unsatisfied, cellpos.n_unsatisfied, MPI_CHAR, total_unsat_freeslots, number_of_unsatisfied_per_th, displs_unsatisfied, MPI_CHAR, MPI_COMM_WORLD);

    //if(myrank==0) for(i=0;i<unsatisfied;i++) printf("\nUNSAT: %c", total_unsat_freeslots[i]);

    //Riempio il resto dell'array degli insoddisfatti con ' ' per far arrivare la size alla stessa misura dei number_of_freeslots_per_th
    for(i=unsatisfied; i<sum_number_freeslots; i++) total_unsat_freeslots[i]=' ';

    //if(myrank==0) for(i=0;i<sum_number_freeslots;i++) printf("\nVALUE: %c", total_unsat_freeslots[i]);
    //if(myrank == 0) for(i=0;i<sum_number_freeslots;i++) printf("%c ",total_unsat_freeslots[i]);
    srand ( round );
    randomize(total_unsat_freeslots, sum_number_freeslots);

    //if(myrank == 0) { printf("\nRANK: %d - POST RANDOMIZE\n", myrank); for(i=0;i<sum_number_freeslots;i++) printf("%c ",total_unsat_freeslots[i]); printf("\n\n");}
    //if(myrank == 1) { printf("\nRANK: %d - POST RANDOMIZE\n", myrank); for(i=0;i<sum_number_freeslots;i++) printf("%c ",total_unsat_freeslots[i]); printf("\n\n");}
    //if(myrank == 2) { printf("\nRANK: %d - POST RANDOMIZE\n", myrank); for(i=0;i<sum_number_freeslots;i++) printf("%c ",total_unsat_freeslots[i]); printf("\n\n");}

    //if(myrank == 0) for(i=0;i<sum_number_freeslots;i++) printf("%c ",total_unsat_freeslots[i]);
     
    //Sostituisco i valori nei freeslots delle varie submatrix con i valori delle sezioni di total_unsatisfied in funzione del proprio rank 
    for(i=0;i<cellpos.n_freeslots;i++) t_mat.submatrix[cellpos.freeslots[i]]=total_unsat_freeslots[displs_freeslots[myrank]+i];
    
    //if(myrank==0) { printf("\nRANK: %d     |",myrank ); for(i=0;i<number_of_freeslots_per_th[myrank];i++) printf("%c ",total_unsat_freeslots[(displs_freeslots[myrank]+i)]); }
    //if(myrank==1) { printf("\nRANK: %d     |",myrank ); for(i=0;i<number_of_freeslots_per_th[myrank];i++) printf("%c ",total_unsat_freeslots[(displs_freeslots[myrank]+i)]); }
    //if(myrank==2) { printf("\nRANK: %d     |",myrank ); for(i=0;i<number_of_freeslots_per_th[myrank];i++) printf("%c ",total_unsat_freeslots[(displs_freeslots[myrank]+i)]); }


}

void swap (char *a, char *b)
{
    char temp = *a;
    *a = *b;
    *b = temp;
}

// A function to generate a random permutation of arr[]
void randomize ( char *arr, int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    
 
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}

Info_cellpositions satisfaction_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int *satisfied){

    int i=0, local_satisfied=0, global_satisfied=0, end_submatrix=0, empty_cells=0;
    cellpos.n_unsatisfied=0;
    cellpos.n_freeslots=0;

    i = (myrank!=0) ? COLUMNS : 0;
    end_submatrix = (myrank==(numproc-1)) ? t_mat.scounts_scatter[myrank] : (t_mat.scounts_scatter[myrank]-COLUMNS);

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
                    if(calc_Ledge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i%COLUMNS==(COLUMNS-1)){
                    if(calc_Redge(t_mat, i)) local_satisfied++; 
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
                    if(calc_center(t_mat, i)) local_satisfied++; 
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
                if(i%COLUMNS==0){
                    if(calc_Ledge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i%COLUMNS==(COLUMNS-1)){
                    if(calc_Redge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i>(t_mat.scounts_scatter[myrank]-COLUMNS+1) && i<(t_mat.scounts_scatter[myrank]-1)){
                    if(calc_Bedge(t_mat, i)) local_satisfied++; 
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
                    if(calc_center(t_mat, i)) local_satisfied++; 
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
                    if(calc_Ledge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else if(i%COLUMNS==(COLUMNS-1)){
                    if(calc_Redge(t_mat, i)) local_satisfied++; 
                    else{
                        cellpos.unsatisfied[cellpos.n_unsatisfied]=t_mat.submatrix[i];
                        cellpos.freeslots[cellpos.n_freeslots]=i;
                        cellpos.n_freeslots++;
                        cellpos.n_unsatisfied++;
                        t_mat.submatrix[i]=' ';
                    }
                }
                else{ 
                    if(calc_center(t_mat, i)) local_satisfied++; 
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

    cellpos.n_freeslots += empty_cells;

    MPI_Allreduce(&local_satisfied, &global_satisfied, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    *satisfied = global_satisfied;

    return cellpos;
}

_Bool calc_firstangle(Info_submatrix mat){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[0]==mat.submatrix[1]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS+1]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_secondangle(Info_submatrix mat){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[COLUMNS-2]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-1]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-2]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_thirdangle(Info_submatrix mat, int myrank){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+1]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+2]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_fourthangle(Info_submatrix mat, int myrank){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-1-COLUMNS]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2-COLUMNS]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_Ledge(Info_submatrix mat, int index){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_Redge(Info_submatrix mat, int index){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_Aedge(Info_submatrix mat, int index){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_Bedge(Info_submatrix mat, int index){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    return (similar_cells>=similarity);
}

_Bool calc_center(Info_submatrix mat, int index){
    float similarity = PERC_SIM*8;
    int similar_cells=0;
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

void create_matrix(char* mat, int E, int O, int X){

    int i=0, val_rand=0;

    //printf("\nPRIMA\nO:%d - X:%d - E:%d\n", O, X, E);

    for(; i<ROWS*COLUMNS; i++){
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
    //printf("\nDOPO\nO:%d - X:%d - E:%d\n", O, X, E);

}

void print_matrix(char *mat) {

    printf("\n");
    for (int i = 0; i < (ROWS * COLUMNS); i++) {
        printf("|");
        switch (mat[i])
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
        
        if (((i + 1) % (COLUMNS) == 0) && (i != 0))
            printf("|\n");
        if ((ROWS * COLUMNS) == 1)
            printf("|\n");
        if (COLUMNS == 1 && i == 0)
            printf("|\n");
    }
    printf("\n");
}

void distribute_matrix(Info_submatrix t_mat, int numproc){
        
        if(numproc>1){

            int size_r_int = ROWS / numproc, i=0;
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
        }else{
            t_mat.scounts_scatter[0] = ROWS * COLUMNS;
            t_mat.displ_scatter[0] = 0;
        }
        
}