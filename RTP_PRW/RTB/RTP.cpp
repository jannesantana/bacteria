#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream> 
#define PI 3.14159265358979323846



const double box = 10;// Size of the 2D box (assuming a square box)
const int N_particles = 2000;  // Number of particles
int T = 200;
double Dt = 0.1;
const double rate=40.0;
const double v0=2;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); 
std::uniform_real_distribution<> dis(0.0, 1.0);// Standard mersenne_twister_engine seeded with rd()

struct Particle {
    double x, y; // Position
    double theta; //angle 
    double ini_posx, ini_posy;
    double aux_posx, aux_posy;
};

double applyPBC(double coord) {
    if (coord >= box) {coord = coord - box;}
    if (coord < 0) {coord = coord + box;}
     return coord;
}
double applyPBCangle(double angle) {
    if (angle >= 2*PI) return angle - 2*PI;
    if (angle < 0) return angle + 2*PI;
    return angle;
}

void initializeParticles(Particle* particles) {
    
    for (int i = 0; i < N_particles; ++i) {

        particles[i].x = 0.2*box; // in the middle
        particles[i].y = 0.2*box;

        particles[i].theta = dis(gen) * 2 * PI;
        
        particles[i].ini_posx = 0.0;
        particles[i].ini_posy = 0.0;
        particles[i].aux_posx=0.0;
        particles[i].aux_posy=0.0;
    }
}


void savePositions(const Particle* particles, const std::string& filename) { 
    std::ofstream file;
   
    file.open(filename, std::ios::app); // open the file 

    // file << "Timestep " << timestep << "\n";
    for (int i = 0; i < N_particles; ++i) {
        file << particles[i].x << " " << particles[i].y << " " << particles[i].theta <<"\n";
    }

    

    file.close();
}


double updatePositions(Particle* particles) {
   
    for (int i = 0; i < N_particles; ++i) {

        double random_nbr;
        random_nbr = dis(gen);
        // std::cout << "Rate " << rate*Dt <<"\n";
        if (rate*Dt > random_nbr) {
            // std::cout << "random number " << random_nbr <<"\n";
            // std::cout << "tumble" << "\n";
            particles[i].theta = dis(gen)*2*PI;
        }
        else {
            // std::cout << "no tumble" << "\n";
            particles[i].theta = particles[i].theta;
        }

        particles[i].x += v0*cos(particles[i].theta)*Dt;
        particles[i].y += v0*sin(particles[i].theta)*Dt;
        particles[i].aux_posx += v0*cos(particles[i].theta)*Dt;
        particles[i].aux_posy += v0*sin(particles[i].theta)*Dt;

        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta = applyPBCangle(particles[i].theta);
        
        }

        double dist=0.0;
        double aux_dist_x=0.0;
        double aux_dist_y=0.0;
        for (int f=0; f <N_particles; f++) {
        aux_dist_x = particles[f].aux_posx - particles[f].ini_posx;
        aux_dist_y = particles[f].aux_posy - particles[f].ini_posy;
        aux_dist_x = aux_dist_x*aux_dist_x;
        aux_dist_y = aux_dist_y*aux_dist_y;
        dist += aux_dist_x + aux_dist_y;
       
    }

return dist;
        }


        int main() {
    Particle* particles = (Particle*)malloc(N_particles * sizeof(Particle)); // use malloc to handle large arrays 
    // int* head = (int*)malloc(Cells * Cells * sizeof(int));
    // int* linkedList = (int*)malloc(N_particles * sizeof(int));


    // if (!particles || !head || !linkedList) {
    //     std::cerr << "Memory allocation failed" << std::endl;
    //     return -1;
    // }
  
    initializeParticles(particles);
   std::stringstream particle_filename;
    particle_filename << "particle_positions_N_particles_" << N_particles
                      << "_box_" << box << "_T_" << T << "_Dt_" << Dt
                      << "_rate_" << rate << "_v0_" << v0 << ".dat";

    std::stringstream sd_filename;
    sd_filename << "squared_disp_N_particles_" << N_particles
                << "_box_" << box << "_T_" << T << "_Dt_" << Dt
                << "_rate_" << rate << "_v0_" << v0 << ".dat";

    std::ofstream particle_file(particle_filename.str());
    particle_file.close();

    std::ofstream sd(sd_filename.str());

        for (int t = 0; t < T; ++t) {

             double store_msd = updatePositions(particles);

            

        // std:: cout << store_msd << "\n";
            sd << (t*Dt) << " " << store_msd/N_particles << "\n";
            savePositions(particles,particle_filename.str());
        }
  free(particles);
    // free(head);
    // free(linkedList);
    sd.close();
    return 0;
        }

