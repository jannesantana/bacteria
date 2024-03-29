// c++ 2nd_model_bact.cc cokus3.c -o 2_bact_model -lgsl -lgslcblas -lm
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

#define TS 100000	/* Conf of particles is saved every TS time steps */
#define T_OP 1		/* Order parameter is saved every T_OP time steps */
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
	
	double TopeX=100;
	double TopeY=100;
	
	int TempsInitial;
	int TempsTotal;
	double Dt; 
	int tSave=0;	


	double Rrep;


	double random_shit;
	double atu_dist;


    double v_ext;
    double v_cont;
	double v_now;
    double rate_ext;
    double rate_cont;
    double k_hook;
	double time_moving;
	double time_not_moving;
    double coin_pili_move;

	double dist;
	

	// WCA
	double ThePow;
	double TheFactor;
	double epsilon=1;
	double sigma=0.5;
	Rrep=0;
	//pow(2,0.1666)*sigma;
	
	if (argc>2) {
	  char * pEnd;

	  N=atoi(argv[1]);
	  TopeX=atoi(argv[2]);
	  TopeY=TopeX;
	  v_ext = strtod(argv[3], &pEnd);
	  v_cont = strtod(argv[4], &pEnd);
	  rate_cont = strtod(argv[5], &pEnd);
      rate_ext = strtod(argv[6], &pEnd);
      k_hook = strtod(argv[7], &pEnd);
      
	
	  TempsTotal=atoi(argv[8]);
	  TempsInitial=0.7*TempsTotal;
	  Dt=strtod(argv[9], &pEnd);
	  
	//   int_strength=strtod(argv[10], &pEnd);

	}
	else {
	  cout << "NbrParticules="; cin >> N;
	  cout << "TopeX="; cin >> TopeX;

	  cout << "TempsTotal="; cin >> TempsTotal;
	  cout << "Dt="; cin >> Dt;
	  TopeY=TopeX;

	  TempsInitial=0.7*TempsTotal;
	}




	double SqrtDeltaT=sqrt(Dt);
	double distancePP=0;
	double distance=0;
	double Dx=0;
	double Dx_unwrap;
    double coin_hook;

	
	

	
	
	//===== Naming files!
	char ending[5];
	char filename[100];
	//=================


if (argc>1) {
	  strcpy(filename,"N_"); strcat(filename,argv[1]);
	  strcat(filename,"_L_"); strcat(filename,argv[2]);
	  strcat(filename,"_vext_"); strcat(filename,argv[3]);
	  strcat(filename,"_vcont_"); strcat(filename,argv[4]);
	  strcat(filename,"_rext_"); strcat(filename,argv[5]);
	  strcat(filename,"_rcont_"); strcat(filename,argv[6]);
	  strcat(filename,"_hook_"); strcat(filename,argv[7]);
      strcat(filename,"_T_"); strcat(filename,argv[8]);
      strcat(filename,"_Dt_"); strcat(filename,argv[9]);
	  strcat(filename,".dat");
	} else {
	  strcpy(filename,"TEST_"); 	  
	  strcat(filename,".dat");
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

//=====================================================
char filename_time_moving[100];
 
	strcpy(filename_time_moving,"time_moving_"); strcat(filename_time_moving,filename);

	ofstream toMoving(filename_time_moving);          			 // open file 
	toMoving.precision(0);
	toMoving.setf(ios::fixed,ios::floatfield);
	
	//===================	
	char filename_time_not_moving[100];
 
	strcpy(filename_time_not_moving,"time_not_moving_"); strcat(filename_time_not_moving,filename);

	ofstream toNOTMoving(filename_time_not_moving);          			 // open file 
	toNOTMoving.precision(0);
	toNOTMoving.setf(ios::fixed,ios::floatfield);
	
	//===================	

	double *RodsSpace = (double*) malloc(N * sizeof(double)); //double Rods[N][1];
	double *newRodsSpace = (double*) malloc(N * sizeof(double)); //double newRods[N][1];
	double *aux_unwrap = (double*) malloc(N * sizeof(double));
	double *ini_unwrap = (double*) malloc(N * sizeof(double));
	double *spring_size = (double*) malloc(N * sizeof(double));
    double *new_spring_size = (double*) malloc(N * sizeof(double));
    double *v_now_all = (double*) malloc(N * sizeof(double));
    

	int i,j,ii,k,t;
    

//INITIAL CONDITIONS 
v_now=v_ext;
for(j=0;j<N;j++) {
	RodsSpace[j] = 1;
    v_now_all[j] = v_now;
    spring_size[j] = v_now*Dt;
}

time_moving = 0;
time_not_moving=0;
for(t=0;t<TempsTotal;t++){

    for(i=0;i<N;i++){

    coin_pili_move = rand()/(double)RAND_MAX;

    if (v_now_all[i] == v_ext) {
        if (coin_pili_move < rate_ext*Dt & spring_size[i] >= 0){
            v_now_all[i] = -v_cont;
        }
		else {
			v_now_all[i] = v_ext;
		}
    }

    if (v_now_all[i]==-v_cont){
        if (coin_pili_move < rate_cont*Dt) {
            v_now_all[i] = v_ext;
        }
		else if (coin_pili_move > rate_cont*Dt & spring_size[i] >=0) {v_now_all[i] = -v_cont;}
    }


    new_spring_size[i] = spring_size[i] + v_now_all[i]*Dt;

    coin_hook = rand()/(double)RAND_MAX;
        if (v_now_all[i] >= 0){
            newRodsSpace[i] = RodsSpace[i];
			aux_unwrap[i]=aux_unwrap[i];
			time_not_moving += 1;
            
        }
		else {
			if (time_not_moving>0) {
			toNOTMoving << time_not_moving << "\n";}
			time_not_moving =0;}
        if (coin_hook < k_hook*Dt & v_now_all[i] < 0) {
            newRodsSpace[i] = RodsSpace[i] + spring;
			aux_unwrap[i]=aux_unwrap[i]+v_now_all[i]*Dt;
			time_moving +=1;
    }
	else {
		if (time_moving >0) {
		toMoving << time_moving << "\n";}
		time_moving=0;
	}


			if (newRodsSpace[i]<=-TopeX/2)		newRodsSpace[i]=TopeX+newRodsSpace[i];
			if (newRodsSpace[i]>=TopeX/2) 	newRodsSpace[i]=newRodsSpace[i]-TopeX;

    }
    
atu_dist=0;
	for(k=0;k<N;k++){
		
		RodsSpace[k]=newRodsSpace[k];
		spring_size[k]=new_spring_size[k];
		

		Dx_unwrap = aux_unwrap[k]-ini_unwrap[k];
		dist=Dx_unwrap*Dx_unwrap;

			atu_dist = atu_dist + dist;
			

			

		   if (t%T_SNAP==0) toXY << (t*Dt) << " " << RodsSpace[k] << "\n";
		   if (t%T_SNAP==0) toSPRING << (t*Dt) << " " << spring_size[k] << "\n";
		   }

		   if (t%T_OP==0) toMSD << (t*Dt) << " " << atu_dist <<"\n";


}


}