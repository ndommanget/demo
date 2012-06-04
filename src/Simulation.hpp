// Simulation.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__


#include "Scene.hpp"
#include "Shaders.hpp"

// Mersenne Twister random generation
#include <../api/random/MersenneTwister.h>
#include <../api/eigen-eigen-3.0.3/Eigen/Dense>



// A container for objects
class Simulation
{
     private:

        // General data
        Scene & scene;                      // Scene in which to draw the simulation
        const Shaders & shaders;            // Shaders to use for rendering
        const unsigned int defaultShader;   // Shader used for test rendering
        const unsigned int spriteShader;    // Shader used for particles
        MTRand * rand;                      // Mersenne-Twister random generator

        // Geometric parameters
        const bool surface;                 // True if 2D simulation
        const float size;                   // Width of the grid
        float h;                            // Cell side length (square or cubes)
        float offset[3];                    // Start point left/bottom/back of the grid
        float offsetInCell;                 // Offset from left/bottom/back of grid cell
        const bool solidWalls;              // True if grid boundary is expected to be solid
        float vectorScale;                  // Scale coefficient to display vector lengths
        
        // Simulation parameters
        const float density;                // Expected density for the fluid
        const float viscosity;              // Expected viscosity for the fluid
        float g[3];                         // Vector of acceleration due to gravity

        // About the data stored on centers of the MAC grid
        const unsigned int nbSamplesX;      // Samples number on X axis
        const unsigned int nbSamplesY;      // Samples number on Y axis
        unsigned int nbSamplesZ;            // Samples number on Z axis
        unsigned int nbSamples;             // Total samples number
        float * samples;                    // Samples positions
        float * pressures;                  // Samples pressures
        float * colors;                     // Samples colors
        float * forces;                     // Samples forces
        unsigned int  * types;              // Types (0 : cell is fluid, 
                                            // 1 : cell is empty, 
                                            // 2 : cell is solid)
           
        // About the data stored on borders of the MAC grid
        unsigned int nbBordersX;            // Number of left/right cell borders
        unsigned int nbBordersY;            // Number of bottom/top cell borders
        unsigned int nbBordersZ;            // Number of back/front cell borders
        float * velocitiesX;                // X velocity components, sampled on left/right borders
        float * velocitiesY;                // Y velocity components, sampled on bottom/top borders
        float * velocitiesZ;                // Z velocity components, sampled on back/front borders
  
        // About the particles used for visualization (advected on the whole grid)
        const unsigned int nbParticlesCoef; // Coefficient for density of particles per cell
        unsigned int nbParticles;           // Number of rendered particles          
                                            // (=nbSamples*[4,8]^nbParticlesCoef)
        float * particles;                  // Particle positions
        float * particleColors;             // Grid interpolated particles colors
        float * particleVelocities;         // Grid interpolated particles velocities
        
        // Objects for OpenGL checking and visualization
        Object * objectSamples;             // To render samples as positions and colors
        Object * objectPressures;           // To render pressures as color
        Object * objectForces;              // To render forces as vectors and colors
        Object * objectTypes;               // To render types as colors
        Object * objectGrid;                // To render the grid as lines    
        Object * objectVelocitiesCenters;   // To render velocities interpolated on samples as vectors and colors
        Object * objectVelocitiesBorders;   // To render velocities components as colors
        Object * objectParticles;           // To render particles as positions and colors 
        Object * objectParticleVelocities;  // To render particles interpolated velocities as vectors and colors

        
        
        // Initialization
        void initSimulation();
        void initVisualization();
        void initSamples();
        void initVelocitiesBorders();
        void initParticles();

        // Elements building
        void buildSamples();
        void buildPressures();
        void buildForces();
        void buildTypes();
        void buildGrid();
        void buildVelocitiesBorders();
        void enforceVelocitiesBorders();
        void buildVelocitiesCenters();
        void buildParticles();
        void buildParticleVelocities();
        
        // Interpolation
        void interpolateFromCenters(const float * const data, const float * const position, float * const result);
        void interpolateFromBorders(const float * const dataX, const float * const dataY, const float * const dataZ, const float * const position, float * const result);
        void forwardEuler1stOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue);
        void rungeKutta2ndOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue);
        void rungeKutta3rdOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue);
        void integrate(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue);
        unsigned int type(int iX, int iY, int iZ);
        void getCell(const float * const pos, int * const iX, int * const iY, int * const iZ);
        
        // System resolution
        void MICPreconditioner(const double * const Adiag, const double * const Aright, const double * const Atop, Eigen::VectorXd  & precon);
        void applyPreconditioner(const double * const Aright, const double * const Atop, const Eigen::VectorXd & precon, const Eigen::VectorXd & r, Eigen::VectorXd & z);
        void multiplySparseMatrix(const double * const Adiag, const double * const Aright, const double * const Atop, const Eigen::VectorXd & v, Eigen::VectorXd & result);
        void conjugateGradient(const double * const Adiag, const double * const Aright, const double * const Atop, const double * const b);
        void setDivergences(double * const divergences);

        // Main steps
        void advectColors(float dt);
        void updateVelocitiesFromBorders();
        void updateVelocitiesFromCenters(const float * const centeredVel);
        void advectVelocities(float dt);
        void applyForces(float dt);
        void project(float dt);
        void advectParticles(float dt);

    public:

        // Constructors / destructor
        Simulation(Scene & scene,
                   const Shaders & shaders,
                   const unsigned int defaultShader=0, 
                   const unsigned int spriteShader=2,
                   const bool surface=true, 
                   const float size=2.0, 
                   const bool solidWalls=true,
                   const float density=1000.0, 
                   const float viscosity=0.0, 
                   const unsigned int nbSamplesX=10, 
                   const unsigned int nbSamplesY=10, 
                   unsigned int nbSamplesZ=1, 
                   const unsigned int nbParticlesCoef=1);
        ~Simulation();

        // Drawing set-up
        unsigned int drawSamples();
        unsigned int drawPressures();
        unsigned int drawForces();
        unsigned int drawTypes();
        unsigned int drawGrid();
        unsigned int drawVelocitiesBorders();
        unsigned int drawVelocitiesCenters();
        unsigned int drawParticles();
        unsigned int drawParticlesVelocities();

        // Main actions
        void update();
        void render();
};


#endif // __SIMULATION_HPP__
