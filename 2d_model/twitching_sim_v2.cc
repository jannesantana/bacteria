// c++ twitching_sim_v2.cc -o twitching_sim_2 
// To compile Mac: c++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include/ -L/opt/homebrew/Cellar/gsl/2.7.1/lib/ twitching_sim_v2.cc cokus3.c -o twitching_sim_2 -lgsl -lgslcblas -lm

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <fstream>
#include <random>
#define PI 3.14159265358979323846

// Constants and definitions 
double box; // Size of the 2D box (assuming a square box)
double CUTOFF;    // Interaction cutoff distance
int N_particles;  // Number of particles
int Cells = box / CUTOFF; // Number of cells per dimension
int T ;
double Dt ; // Duration of each time step
double new_coord;
double k_spring; // spring constand
double rate; // extending rate 
double lo;
double epsilon; //potential well
double sigma;   // potential zero cutoff
// const double sigma_radius=0.8;

std::map<std::string, double> readParams(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, double> params;
    std::string line, param_name;
    double param_value;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> param_name >> param_value)) { break; } // Error in line format
            params[param_name] = param_value;
        }
        file.close();
    } else {
        std::cerr << "Unable to open parameter file: " << filename << std::endl;
        exit(1);
    }
    return params;
}

std::string createFileName(const std::map<std::string, double>& params, const std::string& baseName) {
    std::ostringstream oss;
    oss << baseName;
    for (const auto& param : params) {
        oss << "_" << param.first << "_" << param.second;
    }
    oss << ".dat";
    return oss.str();
}

// Particle structure
struct Particle {
    double x, y; // Position
    double theta; //angle 
    double fx, fy; // colisional force 
    double force_pili_x, force_pili_y;
    double xp,yp;// pili
    double A; //prob of extending 
    int cellIndex; // Cell index in the linked list
    double sd_t;
    double prev_A;  // previous time step probability of extending 
    double ini_posx, ini_posy;
    double aux_posx, aux_posy;
};

    



// periodic boundary conditions
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
    std::normal_distribution<double> dis(0.0, 1.0);// Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis2(0.0,1.0);// Standard mersenne_twister_engine seeded with rd()

// Function to initialize particles with random positions
void initializeParticles(Particle* particles) {
    
    for (int i = 0; i < N_particles; ++i) {

        particles[i].x = dis2(gen) * box; // random initial consitions 
        particles[i].y = dis2(gen) * box;

        particles[i].theta = dis(gen) * 2 * PI;
        particles[i].xp = 0.0;
        particles[i].yp  =0.0;
        particles[i].A = 0.0;
        particles[i].prev_A = 0;
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].sd_t = 0.0;
        particles[i].ini_posx = 0.0;
        particles[i].ini_posy = 0.0;
        particles[i].aux_posx=0.0;
        particles[i].aux_posy=0.0;
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // get cell index 
    }
}



// save poisitons in file 
void savePositions(const Particle* particles, const std::string& file_name) {
    std::ofstream file;
    file.open(file_name, std::ios::app);
    for (int i = 0; i < N_particles; ++i) {
        file << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    file.close();
}

// build the linked list
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

// Function to perform particle interactions
void simulateInteractions(Particle* particles, int* head,  int* linkedList) {
    for (int i = 0; i < N_particles; ++i) {
        
        particles[i].fx = 0.0; // Reset forces
        particles[i].fy = 0.0;
        particles[i].force_pili_x = 0.0;
        particles[i].force_pili_y = 0.0;
       
    }




    for (int i = 0; i < N_particles; ++i) {

        particles[i].prev_A = particles[i].A;


        particles[i].A = dis(gen);
        // double force_pili_x;
        // double force_pili_y;
        double l = sigma + lo;


        // std::cout << "prob of flipping = "<< particles[i].A << "\n";
        // std::cout << "rate = "<< rate*Dt << "\n";
        if ( (particles[i].A <= rate*Dt)  & (particles[i].prev_A <= rate*Dt)) {
            
            // std::cout << "MOVE + KEEP THETA" << "\n";
            particles[i].theta = particles[i].theta;
            particles[i].xp = particles[i].x + l*cos(particles[i].theta);
            particles[i].yp = particles[i].y + l*sin(particles[i].theta);
           
          
            //std::cout << "angle = " << particles[i].theta << "\n";
            particles[i].force_pili_x = -k_spring*(particles[i].x - particles[i].xp);
            particles[i].force_pili_y = -k_spring*(particles[i].y - particles[i].yp);
            //std::cout << "force_x = "<< particles[i].force_pili_x << "\n";
            //std::cout << "force_y = "<< particles[i].force_pili_y << "\n";
            // force_pili_x = -k()

            
        }
        else if ( (particles[i].A <= rate*Dt) & (particles[i].prev_A > rate*Dt) ) {

            // std::cout << "MOVE + NOT KEEP THETA" << "\n";
            particles[i].theta = dis(gen)*2*PI;
            particles[i].xp = particles[i].x + l*cos(particles[i].theta);
            particles[i].yp = particles[i].y + l*sin(particles[i].theta);

            particles[i].force_pili_x = -k_spring*(particles[i].x - particles[i].xp);
            particles[i].force_pili_y = -k_spring*(particles[i].y - particles[i].yp);

        }
        else if (particles[i].A > rate*Dt) {

            // std::cout << "NOT MOVE" << "\n";
            l=0.0;
            // particles[i].prob_A == 0;
            particles[i].force_pili_x = 0.0;
            particles[i].force_pili_y= 0.0;
        
        }


        Particle& pi = particles[i];
        int cellX = static_cast<int>(pi.x / CUTOFF)% Cells;
        int cellY = static_cast<int>(pi.y / CUTOFF)% Cells;
        // std::cout <<"cell index for particle at the moment" << std::endl;

        // Loop over neighboring cells
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighborX = (cellX + dx + Cells) % Cells;
                int neighborY = (cellY + dy + Cells) % Cells;
                int neighborCellIndex = neighborY * Cells + neighborX; // pick up neighboring cell index 

                //  std::cout <<" searching for neighbors" << std::endl;

                int j = head[neighborCellIndex];
                while (j != -1) {
                    if (i != j) {
                        Particle& pj = particles[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double dx_pili = pi.xp - pj.x;
                        double dy_pili = pi.yp - pj.y;

                        // PBC in the distance between the particles!! 
                       if ((dx>0) && (dx>(box/2)))  dx=dx-box;
					        if ((dx<0) && (dx< -(box/2))) dx=dx+box;	
                        if ((dy>0) && (dy>(box/2)))  dy=dy-box;
					        if ((dy<0) && (dy< -(box/2))) dy=dy+box;	
                        
                        if ((dx_pili>0) && (dx_pili>(box/2)))  dx_pili=dx-box;
					        if ((dx_pili<0) && (dx_pili< -(box/2))) dx_pili=dx_pili+box;	
                        if ((dy_pili>0) && (dy_pili>(box/2)))  dy_pili=dy_pili-box;
					        if ((dy_pili<0) && (dy_pili< -(box/2))) dy_pili=dy+box;	
                    //     dx -= box * round(dx / box);
                    //     dy -= box * round(dy / box);

                       

                        double distanceSquared = dx * dx + dy * dy;
                        if (distanceSquared < CUTOFF * CUTOFF) { 

                            // ---- FORCE ---- // 
                            double force=0.0;

                            double distance = sqrt(distanceSquared); // LJ potential 
                            double sigma_12 = pow(sigma,12);

                            double sigma12_over_r_11 = sigma_12* pow(distance, 11);
                             //std::cout << "sigma/r= " << sigma_over_r_12 << "\n";
    
                            double forceMagnitude = 48 * epsilon * (sigma12_over_r_11);     
                            if (distance < pow(2,1/6)*sigma) {force = forceMagnitude + epsilon;}
                            else {force = 0.0;}
                            

                            double fx = force * dx / distance;
                            double fy = force * dy / distance;



                            // double fx = 0.0;
                            // double fy = 0.0;
                            // std::cout << "repulsion force= " <<force << "\n";
                            
                            //std::cout << "distance= " <<distance << "\n"; 
                            
                            pi.fx += fx;
                            pi.fy += fy;
                            pj.fx -= fx; 
                            pj.fy -= fy;
                            // std::cout << "Particle " << i << " interacts with Particle " << j << std::endl;
                        }

                    double distance_pili = sqrt(dx_pili*dx_pili + dy_pili*dy_pili);
                    if (distance_pili <= sigma) {
                        // std::cout << "---------ATTACHED!!!------" << "\n";
                        pi.force_pili_x = 2*pi.force_pili_x;
                        pi.force_pili_y = 2*pi.force_pili_y;
                        // std::cout << "force_x = "<< particles[i].force_pili_x << "\n";
                        // std::cout << "force_y = "<< particles[i].force_pili_y << "\n";
                    }

                    }
                    j = linkedList[j];
                }
            }
        }
    }
}

// UPDATE PARTICLES POSITIONS AND ANGLES 
// double store_msd;
double updatePositions(Particle* particles) {
   
    for (int i = 0; i < N_particles; ++i) {
        // particles[i].A = dis(gen);
        // double force_pili_x;
        // double force_pili_y;
        // double l = sigma_radius + lo;


        // std::cout << "prob of flipping = "<< particles[i].A << "\n";
        // std::cout << "rate = "<< rate*Dt << "\n";
        // if (particles[i].A <= rate*Dt) {
            
        //     particles[i].xp = particles[i].x + l;
        //     particles[i].yp = particles[i].y + l;
        //     particles[i].theta = dis(gen) * PI * 0.5;
        //     force_pili_x = -k_spring*(particles[i].x - particles[i].xp)*cos(particles[i].theta);
        //     force_pili_y = -k_spring*(particles[i].y - particles[i].yp)*sin(particles[i].theta);
        //     std::cout << "force_x = "<< force_pili_x << "\n";
        //     std::cout << "force_y = "<< force_pili_y << "\n";
        //     // force_pili_x = -k()
        // }
        // else {
        //     l=0.0;
        //     force_pili_x = 0.0;
        //     force_pili_y= 0.0;
        // }
       // std::cout <<" BEFORE particle "<< i <<"-> pos (" <<  particles[i].x << ", " <<  particles[i].y << ")" <<  std::endl;


        
        particles[i].x +=particles[i].force_pili_x * Dt + particles[i].fx * Dt;
        particles[i].y += particles[i].force_pili_y * Dt + particles[i].fy * Dt;
        particles[i].aux_posx += particles[i].force_pili_x * Dt + particles[i].fx * Dt;
        particles[i].aux_posy += particles[i].force_pili_y * Dt + particles[i].fy * Dt;
      
        // std:: cout << dist << "\n";
         //std::cout <<" AFTER particle "<< i <<"-> pos (" <<  particles[i].x << ", " <<  particles[i].y << ")" <<  std::endl;
        // Apply periodic boundary conditions
        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta = applyPBCangle(particles[i].theta);
         //std::cout <<" AFTER PBC particle "<< i <<"-> pos (" <<  particles[i].x << ", " <<  particles[i].y << ")" <<  std::endl;
       

        //std::cout <<"particle "<< i <<"-> pos (" <<  particles[i].x << ", " <<  particles[i].y << ")" <<  std::endl;

        // Update cell index
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // update cell index 

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

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string param_file = argv[1];
    // gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
    // gsl_rng_set(r, time(NULL));

    std::map<std::string, double> params = readParams(param_file);


    box = params["box"];
    CUTOFF = params["cutoff"];
    N_particles = static_cast<int>(params["Nparticles"]);
    T = static_cast<int>(params["T"]);
    Dt = params["Dt"];
    epsilon = params["epsilon"];
    k_spring= params["kspring"];
    rate = params["rate"];
    lo = params["lo"];
    epsilon = params["epsilon"];
    sigma = params["sigma"];

    Cells = box / CUTOFF;

    std::string positions_file_name = createFileName(params, "particle_positions");
    std::string squared_disp_file_name = createFileName(params, "squared_disp");




    Particle* particles = (Particle*)malloc(N_particles * sizeof(Particle)); // use malloc to handle large arrays 
    int* head = (int*)malloc(Cells * Cells * sizeof(int));
    int* linkedList = (int*)malloc(N_particles * sizeof(int));


    if (!particles || !head || !linkedList) {
        std::cerr << "Memory allocation failed" << std::endl;
        return -1;
    }
  
    initializeParticles(particles);
    std::ofstream file(positions_file_name);  // Dynamic output file name for particle positions
    file.close();   
    

    std::ofstream sd;
    sd.open(squared_disp_file_name); 

    // sd.close();
        
    for (int t = 0; t < T; ++t) {
       // std::cout << "time= " << t << std::endl << "---------------------" << "\n";
        
  
        buildLinkedList(particles, head, linkedList);


      
   
        simulateInteractions(particles, head, linkedList);



        double store_msd = updatePositions(particles);

            

        // std:: cout << store_msd << "\n";
        sd << (t*Dt) << " " << store_msd/N_particles << "\n";
        

        savePositions(particles, positions_file_name);
        
    }

    free(particles);
    free(head);
    free(linkedList);
    sd.close();
    return 0;
}