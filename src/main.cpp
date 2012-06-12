// main.c
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Application.hpp"



// Entry point in the program
int main(int argc, char *argv[])
{
    //__________________________________________________________________________
    // Main parameters setting

    unsigned int cellsNbWidth=20;
    unsigned int particlesDensity=2;
    bool recording=false;

    for (int iParam=1 ; iParam<argc ; ++iParam)
    {
        if ((iParam<(argc-1)) && (strcmp(argv[iParam], "-c")==0))
        {
            cellsNbWidth=(unsigned int)atoi(argv[iParam+1]);
            iParam++;
        }
        else if ((iParam<(argc-1)) && (strcmp(argv[iParam], "-p")==0))
        {
            particlesDensity=(unsigned int)atoi(argv[iParam+1]);
            iParam++;
        }
        else if (strcmp(argv[iParam], "-r")==0)
        {
            recording=true;
        }
        else
        {
            std::cout<<"Invalid arguments, usage is :"<<std::endl<<"[-c <Cells number on width>]"<<std::endl<<"[-p <particles density [1-3]>]"<<std::endl<<"[-r]>, r optional parameter to record."<<std::endl;
            return 1;
        }
    }
    std::cout<<"Cells number on width : "<<cellsNbWidth<<std::endl;
    std::cout<<"Particles density : "<<particlesDensity<<std::endl;
    if (recording)
        std::cout<<"Recording an image each frame."<<std::endl<<"To make an animation (in ../videos/), run 'sh video [name]'."<<std::endl;
    //__________________________________________________________________________ 
    // Application creation

    Application * application=new Application(recording);

    //__________________________________________________________________________ 
    // Application content definition

    application->setTitle("Fluid simulation");
    application->setSimpleFloor();
    application->setSimulation(cellsNbWidth, particlesDensity);

    //__________________________________________________________________________  
    // Loop

    application->loop();

    //__________________________________________________________________________
    // Cleaning and finishing

    delete application;
	return 0;
}
