#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#define PI 3.14159265358979323846

// Constants
const double BOX_SIZE = 5; // Size of the 2D box (assuming a square box)
const double CUTOFF =1;    // Interaction cutoff distance
const int N_PARTICLES = 1;  // Number of particles
const int N_CELLS = BOX_SIZE / CUTOFF; // Number of cells per dimension
int T = 100;
double TIME_STEP = 0.01; // Duration of each time step
const double vo = 0.8;
double new_coord;

// Particle structure
struct Particle {
    double x, y; // Position
    double theta; //angle 
    double fx, fy; // force
    double l;// pili
    int cellIndex; // Cell index in the linked list
};




// periodic boundary conditions
double applyPBC(double coord) {
    if (coord >= BOX_SIZE) {coord = coord - BOX_SIZE;}
    if (coord < 0) {coord = coord + BOX_SIZE;}
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

    int cellX =  static_cast<int>(x / CUTOFF) % N_CELLS; 
    int cellY  = static_cast<int>(y / CUTOFF) % N_CELLS;

    if (cellX < 0) cellX += N_CELLS; // negative correction
    if (cellY < 0) cellY += N_CELLS;
    if (cellX < 0 || cellX >= N_CELLS || cellY < 0 || cellY >= N_CELLS) { //  in case shit goes wrong 
        
        std::cerr << "Cell index out of bounds: (" << cellX << ", " << cellY << ")" <<  std::endl;
        exit(1);
    }
    return cellY * N_CELLS + cellX;
    
}


    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);// Standard mersenne_twister_engine seeded with rd()
// Function to initialize particles with random positions
void initializeParticles(Particle* particles) {
    
    for (int i = 0; i < N_PARTICLES; ++i) {

        particles[i].x = dis(gen) * BOX_SIZE; // random initial consitions 
        particles[i].y = dis(gen) * BOX_SIZE;
        particles[i].theta = dis(gen) * 2 * PI;
        particles[i].l = 1.0;
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // get cell index 
    }
}

// save poisitons in file 
void savePositions(const Particle* particles) { 
    std::ofstream file;
    file.open("particle_positions.dat", std::ios::app); // open the file 

    // file << "Timestep " << timestep << "\n";
    for (int i = 0; i < N_PARTICLES; ++i) {
        file << particles[i].x << " " << particles[i].y << "\n";
    }

    file.close();
}

// build the linked list
void buildLinkedList(Particle* particles, int* head, int* linkedList) {
    // Initialize head and linked list
    for (int i = 0; i < N_CELLS * N_CELLS; ++i) {
        // std::cout << "Cell " << i << ": ";
        head[i] = -1;
    }
    for (int i = 0; i < N_PARTICLES; ++i) {
        linkedList[i] = -1;
    }

    for (int i = 0; i < N_PARTICLES; ++i) { // initialize the complete linked list 
        int cellIndex = particles[i].cellIndex;
        linkedList[i] = head[cellIndex];
        head[cellIndex] = i;
    }

    for (int i = 0; i < N_CELLS * N_CELLS; ++i) { // neighbouring linked list 
        // std::cout << "Cell " << i << ": ";
        int j = head[i];
        
    }

   

}

// Function to perform particle interactions
void simulateInteractions(Particle* particles, int* head,  int* linkedList) {
    for (int i = 0; i < N_PARTICLES; ++i) {
        particles[i].fx = 0.0; // Reset forces
        particles[i].fy = 0.0;
    }

    const double epsilon = 0.5; //potential well
    const double sigma = 3.0;   // potential zero cutoff


    for (int i = 0; i < N_PARTICLES; ++i) {
        Particle& pi = particles[i];
        int cellX = static_cast<int>(pi.x / CUTOFF)% N_CELLS;
        int cellY = static_cast<int>(pi.y / CUTOFF)% N_CELLS;
        // std::cout <<"cell index for particle at the moment" << std::endl;

        // Loop over neighboring cells
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighborX = (cellX + dx + N_CELLS) % N_CELLS;
                int neighborY = (cellY + dy + N_CELLS) % N_CELLS;
                int neighborCellIndex = neighborY * N_CELLS + neighborX; // pick up neighboring cell index 

                //  std::cout <<" searching for neighbors" << std::endl;

                int j = head[neighborCellIndex];
                while (j != -1) {
                    if (i != j) {
                        Particle& pj = particles[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;

                        // PBC in the distance between the particles!! 
                       if ((dx>0) && (dx>(BOX_SIZE/2)))  dx=dx-BOX_SIZE;
					        if ((dx<0) && (dx<-(BOX_SIZE/2))) dx=dx+BOX_SIZE;	
                        if ((dy>0) && (dy>(BOX_SIZE/2)))  dy=dy-BOX_SIZE;
					        if ((dy<0) && (dy<-(BOX_SIZE/2))) dx=dy+BOX_SIZE;	
                    //     dx -= BOX_SIZE * round(dx / BOX_SIZE);
                    //     dy -= BOX_SIZE * round(dy / BOX_SIZE);

                       

                        double distanceSquared = dx * dx + dy * dy;
                        if (distanceSquared < CUTOFF * CUTOFF) { 

                            // ---- FORCE ---- // 
                            double force;

                            double distance = sqrt(distanceSquared); // LJ potential 
                            double sigma_over_r = sigma / distance;
                            double sigma_over_r_12 = pow(sigma_over_r, 12);
    
                            double forceMagnitude = 4 * epsilon * (sigma_over_r_12);     
                            if (distance < pow(2,1/6)*sigma) {force = forceMagnitude + epsilon;}
                            else {force = 0.0;}
                            

                            double fx = force * dx / distance;
                            double fy = force * dy / distance;
                            //std::cout << "force= " <<forceMagnitude << "\n";
                            // std::cout << "distance= " <<distance << "\n"; 
                            // std::cout << "sigma/r= " << sigma_over_r << "\n";
                            pi.fx += fx;
                            pi.fy += fy;
                            pj.fx -= fx; 
                            pj.fy -= fy;
                            // std::cout << "Particle " << i << " interacts with Particle " << j << std::endl;
                        }
                    }
                    j = linkedList[j];
                }
            }
        }
    }
}

// UPDATE PARTICLES POSITIONS AND ANGLES 
void updatePositions(Particle* particles) {
    for (int i = 0; i < N_PARTICLES; ++i) {
        // std::cout <<"particle "<< i <<"-> pos (" <<  particles[i].x << ", " <<  particles[i].y << ")" <<  std::endl;


        
        particles[i].x += vo*cos(particles[i].theta) * TIME_STEP + particles[i].fx * TIME_STEP;
        particles[i].y += vo*sin(particles[i].theta) * TIME_STEP + particles[i].fy*TIME_STEP;

        // Apply periodic boundary conditions
        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta = applyPBCangle(particles[i].theta);

        // Update cell index
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // update cell index 

    }
}

int main() {
    Particle* particles = (Particle*)malloc(N_PARTICLES * sizeof(Particle)); // use malloc to handle large arrays 
    int* head = (int*)malloc(N_CELLS * N_CELLS * sizeof(int));
    int* linkedList = (int*)malloc(N_PARTICLES * sizeof(int));


    if (!particles || !head || !linkedList) {
        std::cerr << "Memory allocation failed" << std::endl;
        return -1;
    }
  
    initializeParticles(particles);
        std::ofstream file("particle_positions.dat"); // Create a new file or overwrite if it already exists
        file.close();
    for (int t = 0; t < T; ++t) {
        // std::cout << "time= " << t << std::endl << "---------------------" << "\n";
        
  
        buildLinkedList(particles, head, linkedList);


      
   
        simulateInteractions(particles, head, linkedList);



        updatePositions(particles);


        savePositions(particles);
    }

    free(particles);
    free(head);
    free(linkedList);

    return 0;
}