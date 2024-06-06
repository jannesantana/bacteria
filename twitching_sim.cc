#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#define PI 3.14159265358979323846

// Constants
double BOX_SIZE = 10; // Size of the 2D box (assuming a square box)
double CUTOFF = 2.5;    // Interaction cutoff distance
int N_PARTICLES = 9;  // Number of particles
int N_CELLS = BOX_SIZE / CUTOFF; // Number of cells per dimension
int T = 100;
double TIME_STEP = 0.01; // Duration of each time step
double vo = 2;
// Particle structure
struct Particle {
    double x, y; // Position
    double theta; //angle 
    double fx, fy; // force
    double l;// pili
    int cellIndex; // Cell index in the linked list
};

// get the cell index of a particle


// periodic boundary conditions
double applyPBC(double coord) {
    if (coord >= BOX_SIZE) return coord - BOX_SIZE;
    if (coord < 0) return coord + BOX_SIZE;
    return coord;
}

double applyPBCangle(double angle) {
    if (angle >= 2*PI) return angle - 2*PI;
    if (angle < 0) return angle + 2*PI;
    return angle;
}

int getCellIndex(double x, double y) {
    int cellX = static_cast<int>(applyPBC(x) / CUTOFF);
    int cellY = static_cast<int>(applyPBC(y) / CUTOFF);
    if (cellX < 0 || cellX >= N_CELLS || cellY < 0 || cellY >= N_CELLS) {
        std::cerr << "Cell index out of bounds: (" << cellX << ", " << cellY << ")" << std::endl;
        exit(1);
    }
    return cellY * N_CELLS + cellX;
    
}


// Function to initialize particles with random positions
void initializeParticles(Particle* particles) {
    for (int i = 0; i < N_PARTICLES; ++i) {
        particles[i].x = static_cast<double>(rand()) / RAND_MAX * BOX_SIZE;
        particles[i].y = static_cast<double>(rand()) / RAND_MAX * BOX_SIZE;
        particles[i].theta = static_cast<double>(rand()) / RAND_MAX * 2*PI;
        particles[i].l = 0.0;
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y);
    }
}

void savePositions(const Particle* particles) {
    std::ofstream file;
    file.open("particle_positions.dat", std::ios::app); 

    // file << "Timestep " << timestep << "\n";
    for (int i = 0; i < N_PARTICLES; ++i) {
        file << particles[i].x << " " << particles[i].y << "\n";
    }

    file.close();
}

// Function to build the linked list
void buildLinkedList(Particle* particles, int* head, int* linkedList) {
    // Initialize head and linked list
    for (int i = 0; i < N_CELLS * N_CELLS; ++i) {
        // std::cout << "Cell " << i << ": ";
        head[i] = -1;
    }
    for (int i = 0; i < N_PARTICLES; ++i) {
        linkedList[i] = -1;
    }

    for (int i = 0; i < N_PARTICLES; ++i) {
        int cellIndex = particles[i].cellIndex;
        linkedList[i] = head[cellIndex];
        head[cellIndex] = i;
    }

    for (int i = 0; i < N_CELLS * N_CELLS; ++i) {
        std::cout << "Cell " << i << ": ";
        int j = head[i];
        while (j != -1) {
            std::cout << j << " -> ";
            j = linkedList[j];
        }
        std::cout << "NULL" << std::endl;
    }

    // // Populate linked list
    // for (int i = 0; i < N_PARTICLES; ++i) {
    //     int cellIndex = particles[i].cellIndex;
    //     linkedList[i] = head[cellIndex];
    //     head[cellIndex] = i;
    // }

}

// Function to perform particle interactions
void simulateInteractions(Particle* particles, const int* head, const int* linkedList) {
    for (int i = 0; i < N_PARTICLES; ++i) {
        particles[i].fx = 0.0; // Reset forces
        particles[i].fy = 0.0;
    }

    const double epsilon = 1.0; // Depth of the potential well
    const double sigma = 1.0;   // Distance at which the potential is zero


    for (int i = 0; i < N_PARTICLES; ++i) {
        Particle& pi = particles[i];
        int cellX = static_cast<int>(pi.x / CUTOFF);
        int cellY = static_cast<int>(pi.y / CUTOFF);

        // Loop over neighboring cells
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighborX = (cellX + dx + N_CELLS) % N_CELLS;
                int neighborY = (cellY + dy + N_CELLS) % N_CELLS;
                int neighborCellIndex = neighborY * N_CELLS + neighborX;

                int j = head[neighborCellIndex];
                while (j != -1) {
                    if (i != j) {
                        Particle& pj = particles[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;

                        //PBC
                        dx = applyPBC(dx);
                        dy = applyPBC(dy);

                        double distanceSquared = dx * dx + dy * dy;
                        if (distanceSquared < CUTOFF * CUTOFF && distanceSquared >= 0.0) {
                            double distance = sqrt(distanceSquared);
                            double sigma_over_r = sigma / distance;
                            double sigma_over_r_6 = pow(sigma_over_r, 6);
                            double sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
                            double forceMagnitude = 24 * epsilon * (2 * sigma_over_r_12 - sigma_over_r_6) / distanceSquared;    

                            double fx = forceMagnitude * dx;
                            double fy = forceMagnitude * dy;

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

void updatePositions(Particle* particles) {
    for (int i = 0; i < N_PARTICLES; ++i) {
        particles[i].x += vo*cos(particles[i].theta) * TIME_STEP + particles[i].fx * TIME_STEP;
        particles[i].y += vo*sin(particles[i].theta) * TIME_STEP + particles[i].fy*TIME_STEP;

        // Apply periodic boundary conditions
        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta = applyPBCangle(particles[i].theta);

        // Update cell index
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y);
    }
}

int main() {
    Particle* particles = (Particle*)malloc(N_PARTICLES * sizeof(Particle));
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
        std::cout << "time= " << t << std::endl;
        for (int i = 0; i < N_PARTICLES; ++i) {
            std::cout << "Particle " << i << ": (" << particles[i].x << ", " << particles[i].y << ") "
                      << "Cell Index: " << particles[i].cellIndex << " Linked List: " << linkedList[i] << std::endl;
        }
  
        buildLinkedList(particles, head, linkedList);

        // std::cout << " linkedlist" << std::endl;
      
   
        simulateInteractions(particles, head, linkedList);
        // std::cout << "interactions" << std::endl;  


        updatePositions(particles);
        // std::cout << "update" << std::endl;  

        savePositions(particles);
    }

    free(particles);
    free(head);
    free(linkedList);

    return 0;
}