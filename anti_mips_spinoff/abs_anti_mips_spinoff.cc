// 
//  To compile: c++ swarming1D_GN_V1_correct.cc cokus3.c  -o swarming1D -lgsl -lgslcblas -lm

using namespace std;

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>

#include <iostream> 
#include <stdlib.h> 

#include <math.h> 

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define PI 3.14159265358979323846

#define TS 10000	/* Conf of particles is saved every TS time steps */
#define T_OP 10		/* Order parameter is saved every T_OP time steps */
#define T_SNAP 10	/* A snapshot is saved every T_SNAP time steps */


extern void seedMT2();
extern void seedMT(unsigned long int);
extern double ranMT(void);

int main (int argc, const char * argv[]) {

	cout.precision(5);
	
	//For the random number generator
	unsigned int uin = time(NULL);
	
	srand(uin); /* This is to set the dirty random generator */
	seedMT(uin);	/* random numbers seeding FROM COKUS3 */             

	const  gsl_rng_type * T; 
	gsl_rng * rg;



 	time_t seed; seed = time(NULL);	
  	//seed = seed + uin;	
 	seed = uin;	
 	T  = gsl_rng_mt19937;		
 	rg = gsl_rng_alloc (T);		// 
 	gsl_rng_set   (rg,seed);	// 
	srand(uin);

	FILE *fpSaveConf;


	int N;
	
	double TopeX;
	double TopeY;
	
	int TempsInitial;
	int TempsTotal;
	double Dt; 
	int tSave=0;	






	double atu_dist;




	double dist;
	
	
	if (argc>2) {
	  char * pEnd;
// OUTPUTS VARIABLES NAME ENTER   

	//   N=atoi(argv[1]);
	//   TopeX=atoi(argv[2]);
	//   TopeY=TopeX;
	//   v_ext = strtod(argv[3], &pEnd);
	//   v_cont = strtod(argv[4], &pEnd);
	//   rate_cont = strtod(argv[5], &pEnd);
    //   rate_ext = strtod(argv[6], &pEnd);
    //   k_hook = strtod(argv[7], &pEnd);
      
	
	//   TempsTotal=atoi(argv[8]);
	//   TempsInitial=0.7*TempsTotal;
	//   Dt=strtod(argv[9], &pEnd);
	  
	//   int_strength=strtod(argv[10], &pEnd);

	}


	double SqrtDeltaT=sqrt(Dt);
	double distancePP=0;
	double distance=0;
	double Dx=0;
    double TheFactor;
    double Neighbors=0;
    double R=1;
	double R2=R*R; 
	double Rrep;
	double kResort;



    //===== Naming files!
	char ending[5];
	char filename[100];
	//=================


if (argc>1) {
	//   strcpy(filename,"N_"); strcat(filename,argv[1]);
	//   strcat(filename,"_L_"); strcat(filename,argv[2]);
	//   strcat(filename,"_vext_"); strcat(filename,argv[3]);
	//   strcat(filename,"_vcont_"); strcat(filename,argv[4]);
	//   strcat(filename,"_rext_"); strcat(filename,argv[5]);
	//   strcat(filename,"_rcont_"); strcat(filename,argv[6]);
	//   strcat(filename,"_hook_"); strcat(filename,argv[7]);
    //   strcat(filename,"_T_"); strcat(filename,argv[8]);
    //   strcat(filename,"_Dt_"); strcat(filename,argv[9]);
	//   strcat(filename,".dat");

    // String to name the files 

	}


    char filename_MSD[100];
 
	strcpy(filename_MSD,"MSD_"); strcat(filename_MSD,filename);

	ofstream toMSD(filename_MSD);          			 // open file 
	toMSD.precision(5);
	toMSD.setf(ios::scientific,ios::floatfield);
	
	//===================

	char filename_XY[100];
 
	strcpy(filename_XY,"XY_"); strcat(filename_XY,filename);

	ofstream toXY(filename_XY);          			 // open file 
	toXY.precision(5);
	toXY.setf(ios::scientific,ios::floatfield);
//==============================================
	char filename_spring[100];
 
	strcpy(filename_spring,"spring_"); strcat(filename_spring,filename);

	ofstream toSPRING(filename_spring);          			 // open file 
	toSPRING.precision(5);
	toSPRING.setf(ios::scientific,ios::floatfield);

	
	//===================	

	// double *RodsSpace = (double*) malloc(N * sizeof(double)); //double Rods[N][1];
	double *newRodsSpaceX = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
    double *newRodsSpaceY = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
	double *aux_unwrap = (double*) malloc(N * sizeof(double));
	double *ini_unwrap = (double*) malloc(N * sizeof(double));
    double *RodsAngle = (double*) malloc(N * sizeof(double));
	double *newRodsSpace = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
    double *newRodsAngle = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
    double *RodsSpaceX = (double*) malloc(N * sizeof(double)); //double Rods[N][1];
    double *RodsSpaceY = (double*) malloc(N * sizeof(double)); //double Rods[N][1];

	int i,j,ii,k,t;
    int part_now;
    double Dy;

    for (j=0; j<N; j++) 
	{
		//for (i=0; i<2; i++) 
		//{
			RodsSpaceX[j]=0;
            RodsSpaceY[j]=0;
			RodsAngle[j]=0;
			newRodsSpaceX[j]=0;	
            newRodsSpaceY[j]=0;
			newRodsAngle[j]=0;	
		//}
	}
	
		
	for(i=0;i<N;i++)
	{
		RodsSpaceX[i]=(rand()/(double)RAND_MAX)*TopeX;
        RodsSpaceY[i]=(rand()/(double)RAND_MAX)*TopeY;
        
		RodsAngle[i]=(rand()/(double)RAND_MAX)*2*PI; 
	}



	for(t=0;t<TempsTotal;t++)
	{

        for(i=0;i<N;i++){
            Neighbors=0;
            part_now = i;
            

            for (j = 0; i < N; j++)
            {
                if (j!= part_now)
                {
                    Dx = RodsSpaceX[j] - RodsSpaceX[j];
                    Dy = RodsSpaceY[j] - RodsSpaceY[j];

                    if ((Dx>0) && (Dx>(TopeX/2)))  Dx=Dx-TopeX;
					        if ((Dx<0) && (Dx<-(TopeX/2))) Dx=Dx+TopeX;	

                    if ((Dy>0) && (Dy>(TopeY/2)))  Dy=Dy-TopeX;
					        if ((Dy<0) && (Dy<-(TopeY/2))) Dy=Dy+TopeY;	

                distance = sqrt(Dx*Dx + Dy*Dy);
                if ((distance<R) && (distance>0))
                {
                    Neighbors++; 
                }
                

                }
                
            }
            
// UPDATE ESSA MERDA //




        }

        

    }



    return 0;

    }