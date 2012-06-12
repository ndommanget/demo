// main.c
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Application.hpp"



// Entry point in the program
int main(int argc, char **argv)
{
    //__________________________________________________________________________ 
    // Application creation

    Application * application=new Application();

    //__________________________________________________________________________ 
    // Application content definition

    application->setTitle("Fluid simulation");
    application->setSimpleFloor();
    application->setSimulation();

    //__________________________________________________________________________  
    // Loop

    application->loop();

    //__________________________________________________________________________
    // Cleaning and finishing

    delete application;
	return 0;
}
