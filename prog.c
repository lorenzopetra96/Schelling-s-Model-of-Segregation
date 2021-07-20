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

typedef struct{
    //Salvo numero riga inizio e numero riga fine
    int inizio_submatrix;
    int fine_submatrix;

    char *submatrix;



    int *scounts_scatter;
    int *displ_scatter;

    int *scounts_gather;
    int *displ_gather;


} Info_submatrix;

struct info_cellpositions{

    char *unsatisfied; 
    int *freeslots;
    int nfreeslots;

}; 


int calc_firstangle(Info_submatrix mat);
int calc_secondangle(Info_submatrix mat);
int calc_thirdangle(Info_submatrix mat, int);
int calc_fourthangle(Info_submatrix mat, int);
int calc_Ledge(Info_submatrix mat, int);
int calc_Redge(Info_submatrix mat, int);
int calc_Aedge(Info_submatrix mat, int);
int calc_Bedge(Info_submatrix mat, int);
int calc_center(Info_submatrix mat, int); 
void create_matrix(char* mat);
void print_matrix(char *i_mat);
void distribute_matrix(Info_submatrix t_mat, int numproc);
int calculate_satisfaction(Info_submatrix t_mat, int myrank, int numproc);

int main(int argc, char *argv[]){
	
	int numproc, myrank, source, *arr, i, value, round=1, satisfied=0;
    double start, end;
    Info_submatrix t_mat;
    
    char *mat;
	srand(time(NULL)+myrank);

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
            create_matrix(mat);
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

            if((satisfied = calculate_satisfaction(t_mat, myrank, numproc))==100) break;



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

int calculate_satisfaction(Info_submatrix t_mat, int myrank, int numproc){

    int i=0, local_satisfied=0, global_satisfied=0, end_submatrix=0;
    char *unsatisfied = (char *)malloc(sizeof(char)*t_mat.scounts_gather[myrank]); 
    int *freeslots = (int *)calloc(t_mat.scounts_gather[myrank], sizeof(int)), nfreeslots=0;

    i = (myrank!=0) ? COLUMNS : 0;
    end_submatrix = (myrank==(numproc-1)) ? t_mat.scounts_scatter[myrank] : (t_mat.scounts_scatter[myrank]-COLUMNS);

    if(myrank==0){

        for(i=0; i<(t_mat.scounts_scatter[myrank]-COLUMNS);i++){
            if(t_mat.submatrix[i]!=' '){
                if(i==0) local_satisfied += calc_firstangle(t_mat, unsatisfied, freeslots, nfreeslots);
                else if(i==(COLUMNS-1)) local_satisfied += calc_secondangle(t_mat, unsatisfied, freeslots, nfreeslots);
                else if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i>0 && i<COLUMNS) local_satisfied += calc_Aedge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else local_satisfied += calc_center(t_mat, i, unsatisfied, freeslots, nfreeslots);
            }
            else{
                freeslots[nfreeslots]=i;
                nfreeslots++;
            }
        }

    }else if(myrank==(numproc-1)){

        for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
            if(t_mat.submatrix[i]!=' '){
                if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i>(t_mat.scounts_scatter[myrank]-COLUMNS+1) && i<(t_mat.scounts_scatter[myrank]-1)) local_satisfied += calc_Bedge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i==(t_mat.scounts_scatter[myrank]-COLUMNS+1)) local_satisfied += calc_thirdangle(t_mat, myrank, unsatisfied, freeslots, nfreeslots);
                else if(i==(t_mat.scounts_scatter[myrank]-1)) local_satisfied += calc_fourthangle(t_mat, myrank, unsatisfied, freeslots, nfreeslots);
                else local_satisfied += calc_center(t_mat, i, unsatisfied, freeslots, nfreeslots);
            }
            else{
                freeslots[nfreeslots]=i;
                nfreeslots++;
            }
        }

    }else{
        
        for(i=COLUMNS; i<t_mat.scounts_scatter[myrank];i++){
            if(t_mat.submatrix[i]!=' '){
                if(i%COLUMNS==0) local_satisfied += calc_Ledge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else if(i%COLUMNS==(COLUMNS-1)) local_satisfied += calc_Redge(t_mat, i, unsatisfied, freeslots, nfreeslots);
                else local_satisfied += calc_center(t_mat, i, *unsatisfied, *freeslots, nfreeslots);
            }
            else{
                freeslots[nfreeslots]=i;
                nfreeslots++;
            }
        }

    }

    MPI_Allreduce(&local_satisfied, &global_satisfied, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    return global_satisfied;
}

int calc_firstangle(Info_submatrix mat, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[0]==mat.submatrix[1]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS]) similar_cells++;
    if(mat.submatrix[0]==mat.submatrix[COLUMNS+1]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[0];
        freeslots[nfreeslots]=0;
        nfreeslots++;
        mat.submatrix[0]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_secondangle(Info_submatrix mat, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[COLUMNS-2]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-1]) similar_cells++;
    if(mat.submatrix[COLUMNS-1]==mat.submatrix[(COLUMNS*2)-2]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[COLUMNS-1];
        freeslots[nfreeslots]=COLUMNS-1;
        nfreeslots++;
        mat.submatrix[COLUMNS-1]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_thirdangle(Info_submatrix mat, int myrank, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+1]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]==mat.submatrix[mat.scounts_scatter[myrank]-(COLUMNS*2)+2]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1];
        freeslots[nfreeslots]=mat.scounts_scatter[myrank]-COLUMNS+1;
        nfreeslots++;
        mat.submatrix[mat.scounts_scatter[myrank]-COLUMNS+1]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_fourthangle(Info_submatrix mat, int myrank, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*3;
    int similar_cells=0;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-1-COLUMNS]) similar_cells++;
    if(mat.submatrix[mat.scounts_scatter[myrank]-1]==mat.submatrix[mat.scounts_scatter[myrank]-2-COLUMNS]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[mat.scounts_scatter[myrank]-1];
        freeslots[nfreeslots]=mat.scounts_scatter[myrank]-1;
        nfreeslots++;
        mat.submatrix[mat.scounts_scatter[myrank]-1]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_Ledge(Info_submatrix mat, int index, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[index];
        freeslots[nfreeslots]=index;
        nfreeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_Redge(Info_submatrix mat, int index, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[index];
        freeslots[nfreeslots]=index;
        nfreeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_Aedge(Info_submatrix mat, int index, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[index];
        freeslots[nfreeslots]=index;
        nfreeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_Bedge(Info_submatrix mat, int index, char *unsatisfied, int *freeslots, int nfreeslots){
    float similarity = PERC_SIM*5;
    int similar_cells=0;
    if(mat.submatrix[index]==mat.submatrix[index-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS-1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index-COLUMNS+1]) similar_cells++;
    if(mat.submatrix[index]==mat.submatrix[index+1]) similar_cells++;
    if(similar_cells<similarity){
        unsatisfied[nfreeslots]=mat.submatrix[index];
        freeslots[nfreeslots]=index;
        nfreeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}

int calc_center(Info_submatrix mat, int index, char *unsatisfied, int *freeslots, int nfreeslots){
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
        unsatisfied[nfreeslots]=mat.submatrix[index];
        freeslots[nfreeslots]=index;
        nfreeslots++;
        mat.submatrix[index]=' ';
    }
    return (similar_cells>=similarity);
}


void create_matrix(char* mat){

    int i=0, val_rand=0, E = (ROWS*COLUMNS*PERC_E);
    int O = (ROWS*COLUMNS-E)*PERC_O, X = (ROWS*COLUMNS-E)*PERC_X;
    if(ROWS*COLUMNS != (O+E+X)) E+=(ROWS*COLUMNS)-(O+E+X);

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
                t_mat.inizio_submatrix = i==0 ? 0 : t_mat.displ_gather[i];
                t_mat.fine_submatrix = (t_mat.displ_gather[i]+t_mat.scounts_gather[i])/COLUMNS;
                
            }
        }else{
            t_mat.scounts_scatter[0] = ROWS * COLUMNS;
            t_mat.displ_scatter[0] = 0;
        }

        t_mat.nfreeslots = 0;
        t_mat.unsatisfied = 


}