// 
//  To compile: c++ abs_anti_mips_spinoff.cc cokus3.c  -o antimips -lgsl -lgslcblas -lm

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
#define T_SNAP 1	/* A snapshot is saved every T_SNAP time steps */


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
    double Vo;
    double eta;
    double R;
    double disp;



	double atu_dist;




	double dist;
	
	
	if (argc>2) {
	  char * pEnd;
// OUTPUTS VARIABLES NAME ENTER   

	  N=atoi(argv[1]);
	  TopeX=atoi(argv[2]);
	  TopeY=atoi(argv[3]);
	  R = strtod(argv[4], &pEnd);
	  eta = strtod(argv[5], &pEnd);
	    TempsTotal=atoi(argv[6]);
	  TempsInitial=0.7*TempsTotal;
      Dt = strtod(argv[7], &pEnd);
     Vo = strtod(argv[8], &pEnd);
      
	
	//   TempsTotal=atoi(argv[8]);
	//   TempsInitial=0.7*TempsTotal;
	//   Dt=strtod(argv[9], &pEnd);
	  
	//   int_strength=strtod(argv[10], &pEnd);

	}


	double SqrtDeltaT=sqrt(Dt);
	double distancePP=0;
	double distance=0;
	double Dx=0;
     double Dy;
    double TheFactor;
    double Neighbors=0;
    double Alignment;

	double R2=R*R; 
	double Rrep;
	double kResort;
    double rand_nbr;
    double dist0=sqrt(TopeX*TopeX + TopeY*TopeY)/100;



    //===== Naming files!
	char ending[5];
	char filename[100];
	//=================


if (argc>1) {
	  strcpy(filename,"N_"); strcat(filename,argv[1]);
	  strcat(filename,"_Lx_"); strcat(filename,argv[2]);
       strcat(filename,"_Ly_"); strcat(filename,argv[3]);
	  strcat(filename,"_R_"); strcat(filename,argv[4]);
      strcat(filename,"_eta_"); strcat(filename,argv[5]);
      strcat(filename,"_T_"); strcat(filename,argv[6]);
      strcat(filename,"_Dt_"); strcat(filename,argv[7]);
      strcat(filename,"_Vo_"); strcat(filename,argv[8]);
	  strcat(filename,".dat");

    // String to name the files 

	}


    // char filename_MSD[100];
 
	// strcpy(filename_MSD,"MSD_"); strcat(filename_MSD,filename);

	// ofstream toMSD(filename_MSD);          			 // open file 
	// toMSD.precision(5);
	// toMSD.setf(ios::scientific,ios::floatfield);
	
	//===================

	char filename_XY[100];
 
	strcpy(filename_XY,"XY_"); strcat(filename_XY,filename);

	ofstream toXY(filename_XY);          			 // open file 
	toXY.precision(5);
	toXY.setf(ios::scientific,ios::floatfield);
//==============================================
	char Polar[100];
 
	strcpy(Polar,"Polar_Order_"); strcat(Polar,filename);

	ofstream toPolar(Polar);          			 // open file 
	toPolar.precision(5);
	toPolar.setf(ios::scientific,ios::floatfield);

	
	//===================	

	// double *RodsSpace = (double*) malloc(N * sizeof(double)); //double Rods[N][1];
	double *newRodsSpaceX = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
    double *newRodsSpaceY = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
	// double *aux_unwrap = (double*) malloc(N * sizeof(double));
	// double *ini_unwrap = (double*) malloc(N * sizeof(double));
    double *RodsAngle = (double*) malloc(N * sizeof(double));
    double *newRodsAngle = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
    double *RodsSpaceX = (double*) malloc(N * sizeof(double)); //double Rods[N][1];
    double *RodsSpaceY = (double*) malloc(N * sizeof(double)); //double Rods[N][1];

	int i,j,ii,k,t;
    int part_now;
   

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
        // cout << "time=" << t;

        for(i=0;i<N;i++){
            Neighbors=0;
            part_now = i;
            Alignment=0;
            

            for (j = 0; j < N; j++)
            {
                if (j!= part_now)
                {
                    Dx = RodsSpaceX[j] - RodsSpaceX[i];
                    Dy = RodsSpaceY[j] - RodsSpaceY[i];

                    if ((Dx>0) && (Dx>(TopeX/2)))  Dx=Dx-TopeX;
					        if ((Dx<0) && (Dx< -TopeX/2)) Dx=Dx+TopeX;	

                    if ((Dy>0) && (Dy>(TopeY/2)))  Dy=Dy-TopeX;
					        if ((Dy<0) && (Dy<-(TopeY/2))) Dy=Dy+TopeY;	

                distance = sqrt(Dx*Dx + Dy*Dy);
                // if (i==1) cout << "part" << j <<", dis=" << distance << "\n";
                if ((distance<R) && (distance>0))
                {
                    Neighbors++; 
                    Alignment = Alignment + sin(RodsAngle[j]- RodsAngle[i]);
                }
                

                }
                
            }

if (Neighbors ==0) 
{
    rand_nbr = (rand()/(double)RAND_MAX);
    disp = -log(rand_nbr/dist0);
     
}
else {disp = Neighbors;}

            
// UPDATE ESSA MERDA //

newRodsSpaceX[i] = RodsSpaceX[i] + cos(RodsAngle[i])*Vo*(Neighbors+0.1)*Dt ;
newRodsSpaceY[i] = RodsSpaceY[i] + sin(RodsAngle[i])*Vo*(Neighbors+0.1)*Dt;
newRodsAngle[i] = RodsAngle[i] + eta*gsl_ran_gaussian(rg,SqrtDeltaT) + Alignment*Dt/(1+Neighbors);

// if (i==1) cout  << "neighbors=" << Neighbors << "\n" << "--------" << "\n";
 

// BOundary conditions 

            if (newRodsSpaceX[i]<=0)		newRodsSpaceX[i]=TopeX+newRodsSpaceX[i];
		    if (newRodsSpaceX[i]>=TopeX) 	newRodsSpaceX[i]=newRodsSpaceX[i]-TopeX;
            if (newRodsSpaceY[i]<=0)		newRodsSpaceY[i]=TopeY+newRodsSpaceY[i];
		    if (newRodsSpaceY[i]>=TopeY) 	newRodsSpaceY[i]=newRodsSpaceY[i]-TopeY;
			if (newRodsAngle[i]>2*PI) 	newRodsAngle[i]=newRodsAngle[i]-2*PI;
			if (newRodsAngle[i]<0) 		newRodsAngle[i]=newRodsAngle[i]+2*PI;
        }

    for(k=0;k<N;k++)
		 {
		   RodsSpaceX[k]=newRodsSpaceX[k];
           RodsSpaceY[k]=newRodsSpaceY[k];
		   RodsAngle[k]=newRodsAngle[k];

		   if (t%T_SNAP==0) toXY << RodsSpaceX[k] << " " << RodsSpaceY[k] << " " << RodsAngle[k] << "\n";
		 }


    }



    return 0;

    }