// c++ rods_twitching.cc -o twitching_sim_2 
// To compile Mac: clang++ -std=c++17 -stdlib=libc++ rods_twitching.cc -o rods_twitch_sim

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
double L,R;  // rod length and half width
double k_hard; // strength of hardcore repulsion
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

struct Vector2 {
    double x, y;
    Vector2(): x(0), y(0) {}
    Vector2(double _x, double _y): x(_x), y(_y) {}
    Vector2 operator+(const Vector2& o) const { return Vector2(x+o.x, y+o.y); }
    Vector2 operator-(const Vector2& o) const { return Vector2(x-o.x, y-o.y); }
    Vector2 operator*(double s       ) const { return Vector2(x*s,   y*s  ); }
};

inline double dot(const Vector2& a, const Vector2& b) {
    return a.x*b.x + a.y*b.y;
}

inline double clamp(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

void shortestSegmentSegment(
    const Vector2& P1, const Vector2& Q1,
    const Vector2& P2, const Vector2& Q2,
    Vector2& outP, Vector2& outQ, double& outDist)
{
    Vector2 d1 = Q1 - P1;
    Vector2 d2 = Q2 - P2;
    Vector2 r  = P1 - P2;
    double a = dot(d1,d1);
    double e = dot(d2,d2);
    double f = dot(d2,r);

    const double EPS = 1e-12;
    double s = 0, t = 0;

    if (a <= EPS && e <= EPS) {
        // both segments degenerate to points
        outP = P1;
        outQ = P2;
    }
    else if (a <= EPS) {
        // first segment is a point
        s = 0;
        t = clamp(f/e, 0.0, 1.0);
        outP = P1;
        outQ = P2 + d2 * t;
    }
    else if (e <= EPS) {
        // second segment is a point
        t = 0;
        s = clamp(-dot(d1,r)/a, 0.0, 1.0);
        outP = P1 + d1 * s;
        outQ = P2;
    }
    else {
        double b = dot(d1,d2);
        double c = dot(d1,r);
        double denom = a*e - b*b;

        // find s on [0,1]
        if (denom != 0.0)
            s = clamp((b*f - c*e)/denom, 0.0, 1.0);
        else
            s = 0.0;

        // find t on [0,1]
        t = (b*s + f) / e;
        if (t < 0.0)       { t = 0.0; s = clamp(-c/a,        0.0, 1.0); }
        else if (t > 1.0)  { t = 1.0; s = clamp((b-c)/a,     0.0, 1.0); }

        outP = P1 + d1 * s;
        outQ = P2 + d2 * t;
    }

    Vector2 diff = outP - outQ;
    outDist = std::sqrt(dot(diff,diff));
}





// Particle structure
struct Particle {
    double x, y; // Position
    double theta_b; // axis angle of the body
    double theta_p; // pili angle
    double L,R; //rod length and half width
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
    std::uniform_real_distribution<> ang(PI*0.25,3*PI*0.25);// Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> ang2(5*PI*0.25,7*PI*0.25);// Standard mersenne_twister_engine seeded with rd()

// Function to initialize particles with random positions
void initializeParticles(Particle* particles) {
    
    // for (int i = 0; i < N_particles; ++i) {

    //     particles[i].x = ; // random initial conditions 
    //     particles[i].y = dis2(gen) * box;

    //     particles[i].theta_b = dis(gen) * 2 * PI;
    //     particles[i].theta_p = dis(gen) * 2 * PI;

    //     particles[i].xp = 0.0;
    //     particles[i].yp  =0.0;
    //     particles[i].A = 0.0;
    //     particles[i].prev_A = 0;
    //     particles[i].fx = 0.0;
    //     particles[i].fy = 0.0;
    //     particles[i].sd_t = 0.0;
    //     particles[i].ini_posx = 0.0;
    //     particles[i].ini_posy = 0.0;
    //     particles[i].aux_posx=0.0;
    //     particles[i].aux_posy=0.0;
    //     particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // get cell index 
    // }

    int ncol = std::ceil(std::sqrt(N_particles));
    double dx = box / ncol;
    int idx = 0;
    for (int i = 0; i < ncol && idx < N_particles; ++i) {
        

      for (int j = 0; j < ncol && idx < N_particles; ++j) {
        double x = (i + L/2)*dx;
        double y = (j + L/2)*dx;
        // double theta = dis(gen)*2*PI;
        particles[idx].x = applyPBC(x);
        particles[idx].y = applyPBC(y);
        particles[idx].theta_p = dis(gen) * 2 * PI;
        particles[idx].xp = 0.0;
        particles[idx].yp  =0.0;
        particles[idx].A = 0.0;
        particles[idx].prev_A = 0;
        particles[idx].fx = 0.0;
        particles[idx].fy = 0.0;
        particles[idx].sd_t = 0.0;
        particles[idx].ini_posx = 0.0;
        particles[idx].ini_posy = 0.0;
        particles[idx].aux_posx=0.0;
        particles[idx].aux_posy=0.0;
        particles[idx].theta_b = dis(gen) * 2 * PI;
        
        particles[idx].cellIndex = getCellIndex(particles[idx].x, particles[idx].y);
        ++idx;
      }

      
    }
}



// save poisitons in file 
void savePositions(const Particle* particles, const std::string& file_name) {
    std::ofstream file;
    file.open(file_name, std::ios::app);
    for (int i = 0; i < N_particles; ++i) {
        file << particles[i].x << " " << particles[i].y << " " << particles[i].theta_b << " " << particles[i].theta_p << "\n";
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


void resolveOverlaps( Particle* particles,int* head, int* linkedList )
{
// A) Recompute each particle’s cellIndex in O(N):
    for(int i=0; i<N_particles; ++i) {
            particles[i].cellIndex = getCellIndex(particles[i].x,
                        particles[i].y);
        }

// B) Rebuild the cell list in O(N):
buildLinkedList(particles, head, linkedList);

// C) Precompute segment endpoints A[i], B[i] in O(N):
std::vector<Vector2> A(N_particles), B(N_particles);
    for(int i=0; i<N_particles; ++i) {
        double ux = std::cos(particles[i].theta_b),
        uy = std::sin(particles[i].theta_b);
        A[i] = { particles[i].x + 0.5*L*ux,
        particles[i].y + 0.5*L*uy };
        B[i] = { particles[i].x - 0.5*L*ux,
        particles[i].y - 0.5*L*uy };
        }

// D) Loop only over neighboring cells—still ~O(N) or O(N·neighbors):
    for(int cell=0; cell<Cells*Cells; ++cell) {
// for each of the 9 neighboring cells
        for(int dx=-1; dx<=1; ++dx) for(int dy=-1; dy<=1; ++dy) {
            int nbx = (cell % Cells + dx + Cells) % Cells;
            int nby = (cell / Cells + dy + Cells) % Cells;
            int nbrCell = nby*Cells + nbx;

// walk the linked list in cell ‘cell’ and ‘nbrCell’
            for(int i = head[cell]; i!=-1; i = linkedList[i]) {
                double nematic_alignment = 0;
                double neighbors = 1;
                for(int j = head[nbrCell]; j!=-1; j = linkedList[j]) {
                    if (j <= i) continue;    // avoid double‐checking & self

// compute closest points & distance
                Vector2 P, Q; double dist;
                shortestSegmentSegment(A[i], B[i],A[j], B[j],P, Q, dist);
                
                
                double overlap = 2*R - dist;
                if (overlap > 0) {
                    neighbors +=1 ;

                    nematic_alignment += sin(2*(particles[j].theta_b - particles[i].theta_b));

                    // std::cout << "OVERLAP CORRECTION ACTIVATED" << std::endl;
                // unit push direction
                    Vector2 diff = P - Q;
                    Vector2 n;
                    const double EPS = 1e-8;

                    if (dist > EPS) {
                        // normal case
                        n = diff * (1.0 / dist);
                    } else {
                        // fallback: try center-to-center vector
                        double cx = particles[i].x - particles[j].x;
                        double cy = particles[i].y - particles[j].y;
                        double cdist = std::sqrt(cx*cx + cy*cy);

                        if (cdist > EPS) {
                            n = Vector2{cx/cdist, cy/cdist};
                        } else {
                            // completely coincident: pick a random direction
                            double angle = dis(gen) * 2 * PI;  
                            n = Vector2{std::cos(angle), std::sin(angle)};
                        }
                    }

                    // 2) Project apart by half the penetration
                    particles[i].x +=  0.5 * overlap * n.x;
                    particles[i].y +=  0.5 * overlap * n.y;
                    particles[j].x -=  0.5 * overlap * n.x;
                    particles[j].y -=  0.5 * overlap * n.y;

                    // 3) Re‐apply PBC immediately
                    particles[i].x = applyPBC(particles[i].x);
                    particles[i].y = applyPBC(particles[i].y);
                    particles[j].x = applyPBC(particles[j].x);
                    particles[j].y = applyPBC(particles[j].y);

                    }
                }
                if (particles[i].A <= rate*Dt) {
                    particles[i].theta_b += 5*nematic_alignment*Dt/neighbors + 2*sin(particles[i].theta_p - particles[i].theta_b)*Dt;
                }
                else {
                    particles[i].theta_b += 5*nematic_alignment*Dt/neighbors;
                }
                

                particles[i].theta_b = applyPBCangle(particles[i].theta_b);
            }
        }
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


        particles[i].A = dis2(gen);
        // double force_pili_x;
        // double force_pili_y;
        
        // std::cout << "prob of moving before = " << particles[i].prev_A << "\n";
        // std::cout << "prob pof moving after = " << particles[i].A << "\n";
        if ( (particles[i].A <= rate*Dt)  & (particles[i].prev_A <= rate*Dt)) {
            // move and keep angle
            
            particles[i].theta_p = particles[i].theta_p;
            // particles[i].theta_b += sin(particles[i].theta_p - particles[i].theta_b)*Dt;
            // particles[i].theta_b = particles[i].theta_p;
            particles[i].xp = particles[i].x + 0.5*L*cos(particles[i].theta_b) + lo*cos(particles[i].theta_p);
            particles[i].yp = particles[i].y + 0.5*L*sin(particles[i].theta_b) + lo*sin(particles[i].theta_p);
           
        
            particles[i].force_pili_x = -k_spring*(particles[i].x - particles[i].xp);
            particles[i].force_pili_y = -k_spring*(particles[i].y - particles[i].yp);
       

            
        }
        else if ( (particles[i].A <= rate*Dt) & (particles[i].prev_A > rate*Dt) ) {

            std::cout << "MOVE + NOT KEEP THETA" << "\n";
            
            if ( dis2(gen) < 0.6) {
                std::cout << "FORWARDS" << "\n";
                particles[i].theta_p = ang(gen);
            }
            else {
                std::cout << "REVERSAL" << "\n";
                particles[i].theta_p = ang2(gen);
            }
            
            // particles[i].theta_b += sin(particles[i].theta_p - particles[i].theta_b)*Dt;
            particles[i].xp = particles[i].x + 0.5*L*cos(particles[i].theta_b) + lo*cos(particles[i].theta_p);
            particles[i].yp = particles[i].y + 0.5*L*sin(particles[i].theta_b) + lo*sin(particles[i].theta_p);
           
        
            particles[i].force_pili_x = -k_spring*(particles[i].x - particles[i].xp);
            particles[i].force_pili_y = -k_spring*(particles[i].y - particles[i].yp);

        }
        else {
        // if (particles[i].A >= rate*Dt) {
        // if (i==5){
        //     std::cout << "NOT MOVING"<< "\n";
        //     std::cout << "theta_b = "<< particles[i].theta_b <<"\n";
        //     std::cout << "theta_p = "<< particles[i].theta_p <<"\n";}
            particles[i].theta_p = particles[i].theta_p;
            
            particles[i].force_pili_x = 0.0;
            particles[i].force_pili_y= 0.0;
        
        }


        // Particle& pi = particles[i];
        // // Unit vector along rod
        // double ux = cos(pi.theta_b);
        // double uy = sin(pi.theta_b);

        // Vector2 Ai = { particles[i].x + 0.5*particles[i].L * ux,
        //     particles[i].y + 0.5*particles[i].L * uy };
        // Vector2 Bi = { particles[i].x - 0.5*particles[i].L * ux,    
        //     particles[i].y - 0.5*particles[i].L * uy };

        // int cellX = static_cast<int>(pi.x / CUTOFF)% Cells;
        // int cellY = static_cast<int>(pi.y / CUTOFF)% Cells;
        // // std::cout <<"cell index for particle at the moment" << std::endl;

        // // Loop over neighboring cells
        // for (int dx = -1; dx <= 1; ++dx) {
        //     for (int dy = -1; dy <= 1; ++dy) {
        //         int neighborX = (cellX + dx + Cells) % Cells;
        //         int neighborY = (cellY + dy + Cells) % Cells;
        //         int neighborCellIndex = neighborY * Cells + neighborX; // pick up neighboring cell index 

        //         //  std::cout <<" searching for neighbors" << std::endl;

        //         int j = head[neighborCellIndex];
        //         while (j != -1) {
        //             if (i != j) {
        //                 Particle& pj = particles[j];
        //                 double uxj = cos(pj.theta_b);
        //                 double uyj = sin(pj.theta_b);
        //                 Vector2 Aj = { pj.x + 0.5*pj.L * uxj,
        //                     pj.y + 0.5*pj.L * uyj };
        //                 Vector2 Bj = { pj.x - 0.5*pj.L * uxj,    
        //                     pj.y - 0.5*pj.L * uyj };

        //                     Vector2 P, Q;
        //                     double dist;
        //                     shortestSegmentSegment(Ai, Bi, Aj, Bj, P, Q, dist);
        //                     // double Alignment = 0;
        //                     if (dist < 2*R) {
        //                         // std::cout << "COLLISION" << std::endl;
        //                         double overlap = 2*R - dist;
        //                         // Alignment += sin(2*(pi.theta_b - pj.theta_b));
        //                         // normalized push‐apart direction:
        //                         Vector2 n = (P - Q) * (1.0/dist);
        //                         double Fmag = 0.0 * overlap;
        //                         // std::cout << "Fmag = " << Fmag <<std::endl;
                                
        //                         // apply to centers:
        //                         pi.fx +=  Fmag * n.x;
        //                         pi.fy +=  Fmag * n.y;
        //                         pj.fx += -Fmag * n.x;
        //                         pj.fy += -Fmag * n.y;
        //                     }


        //                 // double dx_pili = pi.xp - pj.x;
        //                 // double dy_pili = pi.yp - pj.y;

                     
                        
        //                 // if ((dx_pili>0) && (dx_pili>(box/2)))  dx_pili=dx-box;
		// 			    //     if ((dx_pili<0) && (dx_pili< -(box/2))) dx_pili=dx_pili+box;	
        //                 // if ((dy_pili>0) && (dy_pili>(box/2)))  dy_pili=dy_pili-box;
		// 			    //     if ((dy_pili<0) && (dy_pili< -(box/2))) dy_pili=dy+box;	
                   

        //             }
        //             j = linkedList[j];
        //         }
        //     }
        // }
    }
}

// UPDATE PARTICLES POSITIONS AND ANGLES 
// double store_msd;
double updatePositions(Particle* particles) {
   
    for (int i = 0; i < N_particles; ++i) {
    

        // std::cout << "fx = "<< particles[i].fx * Dt << ", " << "fy = " << particles[i].fy * Dt << std::endl;
        // std::cout << "fx pili = " << particles[i].force_pili_x * Dt << ", " << "fy pili = " << particles[i].force_pili_y * Dt << std::endl;
        
        particles[i].x +=particles[i].force_pili_x * Dt + particles[i].fx * Dt;
        particles[i].y += particles[i].force_pili_y * Dt + particles[i].fy * Dt;
        particles[i].aux_posx += particles[i].force_pili_x * Dt + particles[i].fx * Dt;
        particles[i].aux_posy += particles[i].force_pili_y * Dt + particles[i].fy * Dt;
      
        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta_b = applyPBCangle(particles[i].theta_b);
        particles[i].theta_p = applyPBCangle(particles[i].theta_p);


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
    L = params["rod_length"];
    R = params["rod_radius"];
    k_hard = params["khardcore"];
    k_spring= params["kspring"];
    rate = params["rate"];
    lo = params["lo"];
    
    

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

        resolveOverlaps(particles, head, linkedList);

            

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