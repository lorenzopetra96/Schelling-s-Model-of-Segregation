#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include <time.h>
#include <unistd.h>

#define COLUMNS 10
#define ROWS 10
#define N_ROUND_MAX 10
#define PERC_O 0.5
#define PERC_X (1-PERC_O)
#define PERC_E 0.3
#define PERC_SIM 0.3
#define DELAY 1
#define PRINT_ROUNDS 0

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


_Bool calc_firstangle(Info_submatrix mat, Info_cellpositions cellpos);
_Bool calc_secondangle(Info_submatrix mat, Info_cellpositions cellpos);
_Bool calc_thirdangle(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_fourthangle(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_Ledge(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_Redge(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_Aedge(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_Bedge(Info_submatrix mat, int, Info_cellpositions cellpos);
_Bool calc_center(Info_submatrix mat, int, Info_cellpositions cellpos); 
int sum(int *, int);
void swap(char *, char *);
void randomize(char *, int);
void create_matrix(char* mat, int E, int O, int X);
void print_matrix(char *i_mat);
void distribute_matrix(Info_submatrix t_mat, int numproc);
int satisfaction_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos);
void displacements_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int);


int main(int argc, char *argv[]){
	
	int numproc, myrank, source, *arr, i, value, round=1, satisfied=0;
    double start, end;
    char *mat;
    Info_submatrix t_mat;
    Info_cellpositions cellpos;
    
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
            create_matrix(mat, E, O, X);
            printf("MATRICE INIZIALE\n"); val++;
            print_matrix(mat);
            printf("\n\n");
        }

        distribute_matrix(t_mat, numproc);
        t_mat.submatrix = (char *)malloc(sizeof(char)*t_mat.scounts_scatter[myrank]);

        
        
        MPI_Scatterv(mat, t_mat.scounts_scatter, t_mat.displ_scatter, MPI_CHAR, t_mat.submatrix, t_mat.scounts_scatter[myrank], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while(round<N_ROUND_MAX){

            if(myrank == 0){
                system("clear");
                printf("Numero round: %d\n", round);
            }

            //Calcolo insoddisfatti la cui posizione viene sostituita con ' '. Calcolo freeslots.
            if((satisfied = satisfaction_step(t_mat, myrank, numproc, cellpos))==(ROWS*COLUMNS-E)) break;

            //Riempio array di freeslots con valori di insoddisfatti e randomizzo l'array. Divido l'array
            //in sezioni, ognuna corrispondente a quella del rank. Infine sostituisco i valori dell'array
            //con i valori dei freeslots di ogni processo 
            displacements_step(t_mat, myrank, numproc, cellpos, (ROWS*COLUMNS-E-satisfied));

            if(PRINT_ROUNDS){
                //RICOMPONI E STAMPA MATRICE

            }

            if(myrank == 0){

                system("clear");
                printf("Numero round: %d\nSoddisfatti: %1f %%\n", round, ((float)satisfied*100)/(float)(ROWS*COLUMNS-E));

            }


            round++;
        }


        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();



	}else{
        //aggiungere controllo se numproc>ROWS
        printf("numproc>ROWS... Ripetere");
    }
	
	 
	MPI_Finalize();



	return 0;
}

int sum(int *arr, int length){
    int sum=0;
    for(int i=0;i<length;i++) sum+=arr[i];
    return sum;
}

void displacements_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos, int unsatisfied){

    int n_total_unsat_freeslots=0, displs_unsatisfied[numproc], displs_freeslots[numproc], number_of_unsatisfied_per_th[numproc], number_of_freeslots_per_th[numproc], i=0;
    
    MPI_Allgather(&cellpos.n_unsatisfied, 1, MPI_INT, number_of_unsatisfied_per_th, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Allgather(&cellpos.n_freeslots, 1, MPI_INT, number_of_freeslots_per_th, 1, MPI_INT, MPI_COMM_WORLD);

    n_total_unsat_freeslots = sum(number_of_freeslots_per_th, numproc);
    
    char total_unsat_freeslots[n_total_unsat_freeslots];
    
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

    MPI_Allgatherv(cellpos.unsatisfied, cellpos.n_unsatisfied, MPI_CHAR, total_unsat_freeslots, number_of_unsatisfied_per_th, displs_unsatisfied, MPI_CHAR, MPI_COMM_WORLD);

    //Riempio il resto dell'array degli insoddisfatti con ' ' per far arrivare la size alla stessa misura dei number_of_freeslots_per_th
    for(i=unsatisfied; i<n_total_unsat_freeslots; i++) total_unsat_freeslots[i]=' ';

    randomize(total_unsat_freeslots, n_total_unsat_freeslots);

    //Sostituisco i valori nei freeslots delle varie submatrix con i valori delle sezioni di total_unsatisfied in funzione del proprio rank 
    for(i=0;i<cellpos.n_freeslots;i++) t_mat.submatrix[cellpos.freeslots[i]]=total_unsat_freeslots[displs_freeslots[myrank]+i];
    
    
    free(cellpos.freeslots);
    free(cellpos.unsatisfied);

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
    srand(1);
 
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

int satisfaction_step(Info_submatrix t_mat, int myrank, int numproc, Info_cellpositions cellpos){

    int i=0, local_satisfied=0, global_satisfied=0, end_submatrix=0;
    cellpos.unsatisfied = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]); 
    cellpos.freeslots = (int *)calloc(t_mat.scounts_gather[myrank], sizeof(int));
    cellpos.n_unsatisfied=0;
    cellpos.n_freeslots=0;
    

    i = (myrank!=0) ? COLUMNS : 0;
    end_submatrix = (myrank==(numproc-1)) ? t_mat.scounts_scatter[myrank] : (t_mat.scounts_scatter[myrank]-COLUMNS);

    if(myrank==0){

        for(i=0; i<(t_mat.scounts_scatter[myrank]-COLUMNS);i++){
            if(t_mat.submatrix[i]!=' '){
                if(i==0) local_satisfied += calc_firstangle(t_mat, cellpos);
                else if(i==(COLUMNS-1)) local_satisfied += calc_secondangle(t_mat, cellpos);
                else if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, cellpos);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, cellpos);
                else if(i>0 && i<COLUMNS) local_satisfied += calc_Aedge(t_mat, i, cellpos);
                else local_satisfied += calc_center(t_mat, i, cellpos);
            }
            else{
                cellpos.freeslots[cellpos.n_freeslots]=i;
                cellpos.n_freeslots++;
            }
        }

    }else if(myrank==(numproc-1)){

        for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
            if(t_mat.submatrix[i]!=' '){
                if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, cellpos);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, cellpos);
                else if(i>(t_mat.scounts_scatter[myrank]-COLUMNS+1) && i<(t_mat.scounts_scatter[myrank]-1)) local_satisfied += calc_Bedge(t_mat, i, cellpos);
                else if(i==(t_mat.scounts_scatter[myrank]-COLUMNS+1)) local_satisfied += calc_thirdangle(t_mat, myrank, cellpos);
                else if(i==(t_mat.scounts_scatter[myrank]-1)) local_satisfied += calc_fourthangle(t_mat, myrank, cellpos);
                else local_satisfied += calc_center(t_mat, i, cellpos);
            }
            else{
                cellpos.freeslots[cellpos.n_freeslots]=i;
                cellpos.n_freeslots++;
            }
        }

    }else{
        
        for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
            if(t_mat.submatrix[i]!=' '){
                if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, cellpos);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, cellpos);
                else local_satisfied += calc_center(t_mat, i, cellpos);
            }
            else{
                cellpos.freeslots[cellpos.n_freeslots]=i;
                cellpos.n_freeslots++;
            }
        }

    }
    cellpos.n_unsatisfied=local_satisfied;
    MPI_Allreduce(&local_satisfied, &global_satisfied, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return global_satisfied;
}

_Bool calc_firstangle(Info_submatrix mat, Info_cellpositions cellpos){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[0]==mat.submatrix[1]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS+1]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[0];
        cellpos.freeslots[cellpos.n_freeslots]=0;
        cellpos.n_freeslots++;
        mat.submatrix[0]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_secondangle(Info_submatrix mat, Info_cellpositions cellpos){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[COLUMNS-2]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-1]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-2]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[COLUMNS-1];
        cellpos.freeslots[cellpos.n_freeslots]=COLUMNS-1;
        cellpos.n_freeslots++;
        mat.submatrix[COLUMNS-1]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_thirdangle(Info_submatrix mat, int myrank, Info_cellpositions cellpos){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+1]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+2]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1];
        cellpos.freeslots[cellpos.n_freeslots]=mat.scounts_scatter[myrank]-COLUMNS+1;
        cellpos.n_freeslots++;
        mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_fourthangle(Info_submatrix mat, int myrank, Info_cellpositions cellpos){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-1-COLUMNS]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2-COLUMNS]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[mat.scounts_scatter[myrank]-1];
        cellpos.freeslots[cellpos.n_freeslots]=mat.scounts_scatter[myrank]-1;
        cellpos.n_freeslots++;
        mat.submatrix[mat.scounts_scatter[myrank]-1]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_Ledge(Info_submatrix mat, int index, Info_cellpositions cellpos){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[index];
        cellpos.freeslots[cellpos.n_freeslots]=index;
        cellpos.n_freeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_Redge(Info_submatrix mat, int index, Info_cellpositions cellpos){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[index];
        cellpos.freeslots[cellpos.n_freeslots]=index;
        cellpos.n_freeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_Aedge(Info_submatrix mat, int index, Info_cellpositions cellpos){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[index];
        cellpos.freeslots[cellpos.n_freeslots]=index;
        cellpos.n_freeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_Bedge(Info_submatrix mat, int index, Info_cellpositions cellpos){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[index];
        cellpos.freeslots[cellpos.n_freeslots]=index;
        cellpos.n_freeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

_Bool calc_center(Info_submatrix mat, int index, Info_cellpositions cellpos){
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
    if(similar_cells<similarity){
        cellpos.unsatisfied[cellpos.n_freeslots]=mat.submatrix[index];
        cellpos.freeslots[cellpos.n_freeslots]=index;
        cellpos.n_freeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

void create_matrix(char* mat, int E, int O, int X){

    int i=0, val_rand=0;

    printf("\nPRIMA\nO:%d - X:%d - E:%d\n", O, X, E);

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
    printf("\nDOPO\nO:%d - X:%d - E:%d\n", O, X, E);

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

            int size_r_int = ROWS / numproc;
            int size_r_rest = ROWS % numproc;

            for(int i=0;i<numproc;i++){

                if(i<size_r_rest) t_mat.scounts_scatter[i] = t_mat.scounts_gather[i] = (size_r_int+1) * COLUMNS;
                else t_mat.scounts_scatter[i] = size_r_int * COLUMNS;
                

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
                
                //??????????????????????????????????
                //t_mat.inizio_submatrix = i==0 ? 0 : t_mat.displ_gather[i];
                //t_mat.fine_submatrix = (t_mat.displ_gather[i]+t_mat.scounts_gather[i])/COLUMNS;
                
            }
        }else{
            t_mat.scounts_scatter[0] = ROWS * COLUMNS;
            t_mat.displ_scatter[0] = 0;
        }
}