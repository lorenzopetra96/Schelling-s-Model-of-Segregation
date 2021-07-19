#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include <time.h>
#include <unistd.h>

#define COLONNE 10
#define RIGHE 10
#define N_ROUND_MAX 10
#define PERC_O 0.5
#define PERC_X (1-PERC_O)
#define PERC_E 0.3
#define PERC_SIM 0.3

typedef struct{
    //Salvo numero riga inizio e numero riga fine
    int inizio_sottomat;
    int fine_sottomat;

    char *sottomat;

    

    int *scounts_scatter;
    int *displ_scatter;

    int *scounts_gather;
    int *displ_gather;


} Info_sottomat;


void costruisci_matrice(char* mat);
void stampa_matrice(char *i_mat);
void distribuisci_matrice(Info_sottomat t_mat, int numproc);

int main(int argc, char *argv[]){
	
	int numproc, myrank, source, *arr, i, valore, round=N_ROUND_MAX, satisfied=0;
    double start, end;
    Info_sottomat t_mat;
    
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

    if(numproc<=RIGHE){
        if(myrank == 0){
            int val=1;
            mat = (char *)malloc(RIGHE * COLONNE * sizeof(char));
            costruisci_matrice(mat);
            printf("MATRICE INIZIALE\n"); val++;
            stampa_matrice(mat);
            printf("\n\n");
        }

        distribuisci_matrice(t_mat, numproc);
        t_mat.sottomat = (char *)malloc(sizeof(char)*t_mat.scounts_scatter[myrank]);

        
        
        MPI_Scatterv(mat, t_mat.scounts_scatter, t_mat.displ_scatter, MPI_CHAR, t_mat.sottomat, t_mat.scounts_scatter[myrank], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while(satisfied!=100 || round!=0){





            round--;
        }


        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();


	}else{
        //aggiungere controllo se numproc>RIGHE
        printf("numproc>RIGHE... Ripetere");
    }
	
	 
	MPI_Finalize();



	return 0;
}


void costruisci_matrice(char* mat){

    int i=0, val_rand=0, E = (RIGHE*COLONNE*PERC_E);
    int O = (RIGHE*COLONNE-E)*PERC_O, X = (RIGHE*COLONNE-E)*PERC_X;
    if(RIGHE*COLONNE != (O+E+X)) E+=(RIGHE*COLONNE)-(O+E+X);

    printf("\nPRIMA\nO:%d - X:%d - E:%d\n", O, X, E);

    for(; i<RIGHE*COLONNE; i++){
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

void stampa_matrice(char *mat) {

    printf("\n");
    for (int i = 0; i < (RIGHE * COLONNE); i++) {
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
        
        if (((i + 1) % (COLONNE) == 0) && (i != 0))
            printf("|\n");
        if ((RIGHE * COLONNE) == 1)
            printf("|\n");
        if (COLONNE == 1 && i == 0)
            printf("|\n");
    }
    printf("\n");
}

void distribuisci_matrice(Info_sottomat t_mat, int numproc){

        if(numproc>1){

            int taglia_t_intero = RIGHE / numproc;
            int taglia_t_resto = RIGHE % numproc;

            for(int i=0;i<numproc;i++){

                if(i<taglia_t_resto) t_mat.scounts_scatter[i] = t_mat.scounts_gather[i] = (taglia_t_intero+1) * COLONNE;
                else t_mat.scounts_scatter[i] = taglia_t_intero * COLONNE;
                

                if(i==0){

                    t_mat.scounts_scatter[i] += COLONNE;
                    t_mat.displ_scatter[i] = 0;

                }else if(i==(numproc-1)){

                    t_mat.scounts_scatter[i] += COLONNE;
                    t_mat.displ_scatter[i] = RIGHE*COLONNE - t_mat.scounts_scatter[i];

                }else{

                    t_mat.scounts_scatter[i] += COLONNE*2;
                    t_mat.displ_scatter[i] = t_mat.displ_scatter[i-1] + t_mat.scounts_scatter[i-1] - COLONNE*2;

                }


                t_mat.displ_gather[i] = i==0 ? 0 : t_mat.displ_gather[i-1] + t_mat.scounts_gather[i-1];
                t_mat.inizio_sottomat = i==0 ? 0 : t_mat.displ_gather[i];
                t_mat.fine_sottomat = (t_mat.displ_gather[i]+t_mat.scounts_gather[i])/COLONNE;
                
            }
        }else{
            t_mat.scounts_scatter[0] = RIGHE * COLONNE;
            t_mat.displ_scatter[0] = 0;
        }




}