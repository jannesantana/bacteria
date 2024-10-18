// 
//  To compile: c++ abs_anti_mips_spinoff_ll.cc cokus3.c  -o antimips -lgsl -lgslcblas -lm
// If in mac : c++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include/ -L/opt/homebrew/Cellar/gsl/2.7.1/lib/ abs_anti_mips_spinoff_ll.cc cokus3.c -o antimips -lgsl -lgslcblas -lm
// using namespace std;

// #include <stdio.h>
// #include <string.h>
// #include <fstream>
// #include <iomanip>

// #include <iostream> 
// #include <stdlib.h> 

// #include <math.h> 

// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>

// #define PI 3.14159265358979323846

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#define PI 3.14159265358979323846

// Constants and definitions 
const double box = 10; // Size of the 2D box (assuming a square box)
const double CUTOFF =4.0;    // Interaction cutoff distance
const int N_particles = 200;  // Number of particles
const int Cells = box / CUTOFF; // Number of cells per dimension
int T = 1000;
double Dt = 0.01;
double Vo=1.0;
double eta=0.5;
double R=1.0;
double disp;
// int N;
	
double TopeX;
double TopeY;
	
int TempsInitial;
int TempsTotal;
// double Dt; 
int tSave=0;	




double atu_dist;

double dist;


 // Duration of each time step


#define TS 10000	/* Conf of particles is saved every TS time steps */
#define T_OP 10		/* Order parameter is saved every T_OP time steps */
#define T_SNAP 1	/* A snapshot is saved every T_SNAP time steps */


// extern void seedMT2();
// extern void seedMT(unsigned long int);
// extern double ranMT(void);

struct Particle {
    double x, y; // Position
    double theta; //angle 
    double fx, fy; // colisional force 
    
    int cellIndex; // Cell index in the linked list
    double sd_t;

    double ini_posx, ini_posy;
    double aux_posx, aux_posy;
};

double applyPBC(double coord) {
    if (coord >= box) {coord = coord - box;}
    if (coord < 0) {coord = coord + box;}
     return coord;
}

 
// angular boundary conditions  
double applyPBCangle(double angle) {
    if (angle >= 2*PI) return angle - 2*PI;
    if (angle < 0) return angle + 2*PI;
    return angle;
}

// associate position to cell 
int getCellIndex(double x, double y) {

    int cellX =  static_cast<int>(x / CUTOFF) % Cells; 
    int cellY  = static_cast<int>(y / CUTOFF) % Cells;

    if (cellX < 0) cellX += Cells; // negative correction
    if (cellY < 0) cellY += Cells;
    if (cellX < 0 || cellX >= Cells || cellY < 0 || cellY >= Cells) { //  in case shit goes wrong 
        
        std::cerr << "Cell index out of bounds: (" << cellX << ", " << cellY << ")" <<  std::endl;
        exit(1);
    }
    return cellY * Cells + cellX;
    
}


    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);// Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis2(0.25, 0.75);// Standard mersenne_twister_engine seeded with rd()

void initializeParticles(Particle* particles) {
    
    for (int i = 0; i < N_particles; ++i) {

        particles[i].x = dis(gen) * box; // random initial consitions 
        particles[i].y = dis(gen) * box;

        particles[i].theta = dis(gen) * 2 * PI;
		particles[i].fx = 0.0;
        particles[i].fy = 0.0;

		 particles[i].ini_posx = 0.0;
        particles[i].ini_posy = 0.0;
        particles[i].aux_posx=0.0;
        particles[i].aux_posy=0.0;
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // get cell index 

	}

	}

void savePositions(const Particle* particles) { 
    std::ofstream file;
   
    file.open("particle_positions.dat", std::ios::app); // open the file 

    // file << "Timestep " << timestep << "\n";
    for (int i = 0; i < N_particles; ++i) {
        file << particles[i].x << " " << particles[i].y << " " << particles[i].xp << " " << particles[i].yp << "\n";
    }

    

    file.close();
}


void buildLinkedList(Particle* particles, int* head, int* linkedList) {
    // Initialize head and linked list
    for (int i = 0; i < Cells * Cells; ++i) {
        // std::cout << "Cell " << i << ": ";
        head[i] = -1;
    }
    for (int i = 0; i < N_particles; ++i) {
        linkedList[i] = -1;
    }

    for (int i = 0; i < N_particles; ++i) { // initialize the complete linked list 
        int cellIndex = particles[i].cellIndex;
        linkedList[i] = head[cellIndex];
        head[cellIndex] = i;
    }

    for (int i = 0; i < Cells * Cells; ++i) { // neighbouring linked list 
        // std::cout << "Cell " << i << ": ";
        int j = head[i];
        
    }

   

}

void simulateInteractions(Particle* particles, int* head,  int* linkedList) {
    for (int i = 0; i < N_particles; ++i) {
        
        particles[i].fx = 0.0; // Reset forces
        particles[i].fy = 0.0;
        // particles[i].force_pili_x = 0.0;
        // particles[i].force_pili_y = 0.0;
       
    }

 for (int i = 0; i < N_particles; ++i) {

	
 }
