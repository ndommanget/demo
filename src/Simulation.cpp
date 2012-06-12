// Simulation.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Simulation.hpp"
#include "Tools.hpp"
#include "Scene.hpp"



//______________________________________________________________________________
// Initialization


// Inits the data sampled on centers and borders of the MAC grid
void Simulation::initSimulation()
{
    // Inits the samples positions and colors
    this->initSamples();
    
    // Inits the cell borders velocity components
    this->initVelocitiesBorders();
}


// Inits the particles spread over the grid for visualization
void Simulation::initVisualization()
{
    this->initParticles();
}


// Inits the samples positions and marks 4 zones with different colors
void Simulation::initSamples()
{
	for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
	{
	    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
	    {
	        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
	        {
	            unsigned int index=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
	            samples[index*4+0]=this->offset[0]+this->offsetInCell+iX*this->h;
	            samples[index*4+1]=this->offset[1]+this->offsetInCell+iY*this->h;
	            samples[index*4+2]=this->offset[2]+this->offsetInCell+iZ*this->h;
	            samples[index*4+3]=1.0;
	            
	            pressures[index]=1.0;
	            
	            types[index]=1;
	            if ( (iX>=(nbSamplesX/2))
                  && (iY>=(nbSamplesY/2))
                  && (iZ<=(nbSamplesZ/2)) )
	                types[index]=0;

	            if ( (iX>=(1*nbSamplesX/2)) && (iX<=(2*nbSamplesX/3)) 
	              && (iY>=(1*nbSamplesY/2)) && (iY<=(2*nbSamplesY/3))
                  && (iZ>=(1*nbSamplesZ/2)) && (iZ<=(2*nbSamplesZ/3)) )
                    types[index]=2;
	              
	            forces[index*4+0]=g[0];
	            //forces[index*4+1]=0.0;
	            forces[index*4+1]=g[1];
	            forces[index*4+2]=g[2];
	            forces[index*4+3]=0.0;
	            //if (iX==(nbSamplesX/4))   forces[index*4+1]=10.0;
	            //if (iX==(3*nbSamplesX/4)) forces[index*4+1]=-10.0;	 
	            if ((iX>=(nbSamplesX/2)) && (iX<=(nbSamplesX/2+2))) forces[index*4+1]+=12.0;	           
	            colors[index*4+0]=0.0;
	            colors[index*4+1]=0.0;
	            colors[index*4+2]=0.0;
	            colors[index*4+3]=1.0;
	            // Right : red
	            if (samples[index*4+0]>0.0) colors[index*4+0]=1.0;
	            // Top : Green
	            if (samples[index*4+1]>0.0) colors[index*4+1]=1.0;
	            // Bottom-left : blue / Top-right : white
	            if ((int(colors[index*4+0]+colors[index*4+1])%2)==0)
	                colors[index*4+2]=1.0;
	        }
	    }
	}
}


// Inits velocity directional components (sampled on cell borders)
void Simulation::initVelocitiesBorders()
{
    float coefVelocity=1.0/1.0;

	unsigned int iBordersX, iBordersY, iBordersZ;
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    	{
    	    for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
    	    {   
    	        iBordersX=iZ*(nbSamplesY*(nbSamplesX+1))+iY*(nbSamplesX+1)+iX;
    	        
    	        // Inits circular velocity field (velocityX=-borderPosition[1])
                //velocitiesX[iBordersX]=(-(offset[1]+offsetInCell+iY*h)/size)*coefVelocity;
                
                // Inits diagonally growing velocity field (velocityX=k*borderPosition[0])
                //velocitiesX[iBordersX]=(((float)iX/(float)nbSamplesX)/size)*coefVelocity;
                
                velocitiesX[iBordersX]=0.0*coefVelocity;
            }
        }
    }
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
    	{
    	    for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
    	    {   
    	        iBordersY=iZ*((nbSamplesY+1)*nbSamplesX)+iY*nbSamplesX+iX;
    	        
    	        // Inits circular velocity field (velocityY=borderPosition[0])
                //velocitiesY[iBordersY]=((offset[0]+offsetInCell+iX*h)/size)*coefVelocity;
                
                // Inits diagonally growing velocity field (velocityY=k*borderPosition[1])
                //velocitiesY[iBordersY]=(((float)iY/(float)nbSamplesY)/size)*coefVelocity;
                
                velocitiesY[iBordersY]=0.0*coefVelocity;
            }
        }
    }

    for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {   
                iBordersZ=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                    
                // Inits circular velocity field (velocityY=borderPosition[0])
                //velocitiesZ[iBordersZ]=((offset[0]+offsetInCell+iX*h)/size)*coefVelocity;
                    
                // Inits diagonally growing velocity field (velocityY=k*borderPosition[1])
                //velocitiesZ[iBordersZ]=(((float)iZ/(float)nbSamplesZ)/size)*coefVelocity;
                    
                velocitiesZ[iBordersZ]=0.0*coefVelocity;
            }
        }
    }
    enforceVelocitiesBorders();
}


// Init the particles positions, and interpolates colors and velocities from grid
void Simulation::initParticles()
{
    unsigned int nbAxis=pow(2.0, nbParticlesCoef);  
    
    float halfCell=h/2.0;
    float step=h/float(nbAxis);
    float sideOffsetX=h/(4.0*nbAxis);
    
    unsigned int iParticles=0;
    float * particlesTmp=new float[this->nbParticles*4]; //new
        
    // Inits 8 particles in 3D cells (optimal distribution)
    // if (nbParticlesCoef==1) : corners of a losange in the square cell
    // if (nbParticlesCoef==2) : twice as more ...
	for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
	{
	    if (types[iSamples]==0)
	    {
            for (unsigned int iZ=0 ; iZ<nbAxis ; ++iZ)
            {
                for (unsigned int iY=0 ; iY<nbAxis ; ++iY)
                {
                    float offsetY=-halfCell + step/2.0 + iY*step;
                    for (unsigned int iX=0 ; iX<nbAxis ; ++iX)
                    {
                        float offsetZ=-(-halfCell  + sideOffsetX + iZ*step);
                        float offsetX=-halfCell + sideOffsetX + iX*step;
                        if (iY%2==1)
                        {
                            offsetZ-=step/2.0;
                            offsetX+=step/2.0;
                        }
                        particlesTmp[iParticles*4+0]=samples[iSamples*4+0]+offsetX;
                        particlesTmp[iParticles*4+1]=samples[iSamples*4+1]+offsetY;
                        particlesTmp[iParticles*4+2]=samples[iSamples*4+2]+offsetZ;
                        particlesTmp[iParticles*4+3]=samples[iSamples*4+3];
                        iParticles++;
                    }
                }
            }
	    }
	}

	nbParticles=iParticles;
	delete [] particles;
    delete [] particleColors;
    delete [] particleVelocities;
    this->particles=new float[this->nbParticles*4];
    this->particleColors=new float[this->nbParticles*4];
    this->particleVelocities=new float[this->nbParticles*4*2];
    
    for (unsigned int iParticlesCoord=0 ; iParticlesCoord<4*this->nbParticles ; ++iParticlesCoord)
        this->particles[iParticlesCoord]=particlesTmp[iParticlesCoord];
    delete [] particlesTmp;
    
    // Indices to draw depth sorted particles
    this->particleIndices.resize(nbParticles, 0);
    for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
        particleIndices[iParticles]=iParticles; 
     

	for (unsigned int iParticles=0 ; iParticles<this->nbParticles ; ++iParticles)
	{
	    // Moves particles randomly and at close range around initial position
        float radius=getRand(*this->rand, true, 0.0, 1.0);
        float t=getRand(*this->rand, true, 0.0, 2.0*M_PI);
        float z=getRand(*this->rand, true, -1.0, 1.0);

        radius=2.0*sideOffsetX*pow(radius, 1.0/3.0);
        float r=sqrt(1.0-z*z);
        particles[iParticles*4+0]+=radius*r*cos(t);
        particles[iParticles*4+1]+=radius*r*sin(t);
        particles[iParticles*4+2]+=-radius*z;

        // Interpolates color and velocity bilinearly for each particle
        for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        {
            particleColors[iParticles*4+iCoord]=0.0;
            if (iCoord==3) particleColors[iParticles*4+iCoord]=1.0;
            particleVelocities[iParticles*4+iCoord]=0.0;
        }

        //interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0]));
        /*particleColors[iParticles*4+0]=1.0;
        particleColors[iParticles*4+3]=1.0;
        if (iParticles<2*nbParticles/3) particleColors[iParticles*4+1]=1.0;
        if (iParticles<1*nbParticles/3) particleColors[iParticles*4+2]=1.0;*/

        particleColors[iParticles*4+0]=2.0*(particles[iParticles*4+0]/(h*(float)nbSamplesX));
        particleColors[iParticles*4+1]=2.0*(particles[iParticles*4+1]/(h*(float)nbSamplesY));
        particleColors[iParticles*4+2]=2.0*(particles[iParticles*4+2]/(h*(float)nbSamplesZ));
        particleColors[iParticles*4+3]=1.0; 

        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), &(particleVelocities[iParticles*4+0]));  
    }
}


//______________________________________________________________________________
// Elements building


// Builds an object to visualize the samples (positions and colors)
// The object will then be accessed and updated each frame (advectColors(...))
void Simulation::buildSamples()
{
    //std::cout<<"Building simulation : samples positions."<<std::endl;
    
    objectSamples->setNbVertices(nbSamples);
    
    // No indices are necessary since we draw GL_POINTS
	unsigned int * indices=NULL;
	
    // Sends the data into buffers on the GPU
    objectSamples->sendPrimitives(samples, indices);
    objectSamples->sendColors(colors, true);
}


// Builds an object to visualize the pressures
// The object will then be accessed and updated each frame (project(...))
void Simulation::buildPressures()
{
    //std::cout<<"Building simulation : pressures."<<std::endl;
    
    objectPressures->setNbVertices(nbSamples);
    // No indices are necessary since we draw GL_POINTS
	unsigned int * indices=NULL;
	
    float colors[objectPressures->getNbVertices()*4];
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
    {
            colors[iSamples*4+0]=pressures[iSamples];
            colors[iSamples*4+1]=pressures[iSamples];
            colors[iSamples*4+2]=pressures[iSamples];
            colors[iSamples*4+3]=1.0;
    }
    // Sends the data into buffers on the GPU
    objectPressures->sendPrimitives(samples, indices);
    objectPressures->sendColors(colors, true);
}


// Builds an object to visualize the forces
void Simulation::buildForces()
{
    //std::cout<<"Building simulation : cell forces."<<std::endl;
    
    objectForces->setNbVertices(nbSamples*2);
    objectForces->setNbIndices(objectForces->getNbVertices());
    
    float vertices[objectForces->getNbVertices()*4];
    float colors[objectForces->getNbVertices()*4];
    unsigned int indices[objectForces->getNbIndices()];
    
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                unsigned int index=iSamples*2+0;
                // Position of vectors beginning
                vertices[index*4+0]=samples[iSamples*4+0];
                vertices[index*4+1]=samples[iSamples*4+1];
                vertices[index*4+2]=samples[iSamples*4+2];
                vertices[index*4+3]=1.0;
                colors[index*4+0]=1.0;
                colors[index*4+1]=1.0;
                colors[index*4+2]=1.0;
                colors[index*4+3]=1.0;
                indices[index]=index;
                        
                index=iSamples*2+1;

                // The vector is printed from white to red
                vertices[index*4+0]=samples[iSamples*4+0]+forces[iSamples*4+0]*0.01;
                vertices[index*4+1]=samples[iSamples*4+1]+forces[iSamples*4+1]*0.01;
                vertices[index*4+2]=samples[iSamples*4+2]+forces[iSamples*4+2]*0.01;
                vertices[index*4+3]=1.0;
                colors[index*4+0]=1.0;
                colors[index*4+1]=0.0;
                colors[index*4+2]=0.0;
                colors[index*4+3]=1.0;
                indices[index]=index;
            }  
        }
    }
    // Sends the data into buffers on the GPU
    objectForces->sendPrimitives(vertices, indices);
    objectForces->sendColors(colors);
}


// Builds an object to visualize the types
void Simulation::buildTypes()
{
    //std::cout<<"Building simulation : samples types."<<std::endl;
    
    objectTypes->setNbVertices(nbSamples);
    // No indices are necessary since we draw GL_POINTS
	unsigned int * indices=NULL;
	
    float colors[objectTypes->getNbVertices()*4];
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
    {
            colors[iSamples*4+0]=0.0;
            colors[iSamples*4+1]=0.0;
            colors[iSamples*4+2]=0.0;
            colors[iSamples*4+3]=1.0;
            colors[iSamples*4+types[iSamples]]=1.0;
    }
    // Sends the data into buffers on the GPU
    objectTypes->sendPrimitives(samples, indices);
    objectTypes->sendColors(colors, true);
}


// Builds an object to visualize the grid (borders of cells as lines)
void Simulation::buildGrid()
{
    //std::cout<<"Building simulation : grid borders."<<std::endl;

    objectGrid->setNbVertices((nbSamplesX+1)*(nbSamplesY+1)*(nbSamplesZ+1));

    unsigned int nbLinesX=this->nbSamplesX*(this->nbSamplesY+1)*(this->nbSamplesZ+1);
    unsigned int nbLinesY=(this->nbSamplesX+1)*this->nbSamplesY*(this->nbSamplesZ+1);
    unsigned int nbLinesZ=(this->nbSamplesX+1)*(this->nbSamplesY+1)*this->nbSamplesZ;   
    objectGrid->setNbIndices((nbLinesX+nbLinesY+nbLinesZ)*2);

    float vertices[objectGrid->getNbVertices()*4];
    unsigned int indices[objectGrid->getNbIndices()];
    
    unsigned int iIndices=0;
    for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
        {
            for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
            {
                unsigned int index=iZ*((nbSamplesY+1)*(nbSamplesX+1))
                                  +iY*                (nbSamplesX+1)
                                  +iX;
                // Builds a vertex at each cell corner
                vertices[index*4+0]=offset[0]+iX*h;
                vertices[index*4+1]=offset[1]+iY*h;
                vertices[index*4+2]=offset[2]+iZ*h;
                vertices[index*4+3]=1.0;
                
                // Builds x oriented lines
                if (iX<nbSamplesX) // x oriented segments
                {
                    indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
                                       + iY   *                (nbSamplesX+1)
                                       + iX;
                    indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
                                       + iY   *                (nbSamplesX+1)
                                       +(iX+1);
                }
                // Builds y oriented lines
                if (iY<nbSamplesY) // y oriented segments
                {
                    indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
                                       + iY   *                (nbSamplesX+1)
                                       + iX;
                    indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
                                       +(iY+1)*                (nbSamplesX+1)
                                       + iX;
                }
                // Builds z oriented lines
                if (iZ<nbSamplesZ) // z oriented segments
                {
                    indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
                                       + iY   *                (nbSamplesX+1)
                                       + iX;
                    indices[iIndices++]=(iZ+1)*((nbSamplesY+1)*(nbSamplesX+1))
                                       + iY*                   (nbSamplesX+1)
                                       + iX;
                }
            }
        }
    }
    // Sends the data into buffers on the GPU
    objectGrid->sendPrimitives(vertices, indices);
}


// Builds an object to visualize the velocity components as colors (at borders)
// The object will then be accessed and updated each frame (advectVelocities(...))
void Simulation::buildVelocitiesBorders()
{
    //std::cout<<"Building simulation : cell borders velocities."<<std::endl;
    
    objectVelocitiesBorders->setNbVertices(nbBordersX+nbBordersY+nbBordersZ);
	unsigned int * indices=NULL;
	
    float vertices[objectVelocitiesBorders->getNbVertices()*4];
    float colors[objectVelocitiesBorders->getNbVertices()*4];

	unsigned int iBorders=0;
	unsigned int iBordersX, iBordersY, iBordersZ;
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    	{
    	    for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
    	    {   
    	        iBordersX=iZ*(nbSamplesY*(nbSamplesX+1))+iY*(nbSamplesX+1)+iX;
                // Velocity x components are sampled on the middle 
                // of the yz oriented cell borders
                vertices[iBorders*4+0]=offset[0]+iX*h;
                vertices[iBorders*4+1]=offset[1]+offsetInCell+iY*h;
                vertices[iBorders*4+2]=offset[2]+offsetInCell+iZ*h;
                vertices[iBorders*4+3]=1.0;
                
                // Velocity x components are black to red
                float normVelocityX=(velocitiesX[iBordersX]);
                if (normVelocityX<0.0) normVelocityX=-normVelocityX;
                colors[iBorders*4+0]=normVelocityX;
                colors[iBorders*4+1]=0.0;
                colors[iBorders*4+2]=0.0;
                colors[iBorders*4+3]=1.0;
                
                iBorders++;
            }
        }
    }
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
    	{
    	    for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
    	    {   
                iBordersY=iZ*((nbSamplesY+1)*nbSamplesX)+iY*nbSamplesX+iX;
                // Velocity y components are sampled on the middle 
                // of the xz oriented cell borders
                vertices[iBorders*4+0]=offset[0]+offsetInCell+iX*h;
                vertices[iBorders*4+1]=offset[1]+iY*h;
                vertices[iBorders*4+2]=offset[2]+offsetInCell+iZ*h;
                vertices[iBorders*4+3]=1.0;
                
                // Velocity y components are black to green
                float normVelocityY=(velocitiesY[iBordersY]);
                if (normVelocityY<0.0) normVelocityY=-normVelocityY;
                colors[iBorders*4+0]=0.0;
                colors[iBorders*4+1]=normVelocityY;
                colors[iBorders*4+2]=0.0;
                colors[iBorders*4+3]=1.0;

                iBorders++;
            }
        }
    }
    for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {   
                iBordersZ=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                // Velocity y components are sampled on the middle 
                // of the xz oriented cell borders
                vertices[iBorders*4+0]=offset[0]+offsetInCell+iX*h;
                vertices[iBorders*4+1]=offset[1]+offsetInCell+iY*h;
                vertices[iBorders*4+2]=offset[2]+iZ*h;
                vertices[iBorders*4+3]=1.0;
                    
                // Velocity y components are black to green
                float normVelocityZ=(velocitiesZ[iBordersZ]);
                if (normVelocityZ<0.0) normVelocityZ=-normVelocityZ;
                colors[iBorders*4+0]=0.0;
                colors[iBorders*4+1]=0.0;
                colors[iBorders*4+2]=normVelocityZ;
                colors[iBorders*4+3]=1.0;

                iBorders++;
            }
        }
    }

    // Sends the data into buffers on the GPU
    objectVelocitiesBorders->sendPrimitives(vertices, indices);
    objectVelocitiesBorders->sendColors(colors, true);
}


// Cancels normal velocities on solid boundaries
void Simulation::enforceVelocitiesBorders()
{
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                // If solid cell
                if (types[iSamples]==2)   
                { 
                    //  normal velocity component must be null
                    velocitiesX[indexVelocitiesXLeft]  =0.0;
                    velocitiesX[indexVelocitiesXRight] =0.0;
                    velocitiesY[indexVelocitiesYBottom]=0.0;
                    velocitiesY[indexVelocitiesYTop]   =0.0;
                    velocitiesZ[indexVelocitiesZBack]  =0.0;
                    velocitiesZ[indexVelocitiesZFront] =0.0;
                }
                
                if (solidWalls)
                {
                    if (iX==0)              velocitiesX[indexVelocitiesXLeft]  =0.0;
                    if (iX==(nbSamplesX-1)) velocitiesX[indexVelocitiesXRight] =0.0;
                    if (iY==0)              velocitiesY[indexVelocitiesYBottom]=0.0;
                    if (iY==(nbSamplesY-1)) velocitiesY[indexVelocitiesYTop]   =0.0;
                    if (iZ==0)              velocitiesZ[indexVelocitiesZBack]  =0.0;
                    if (iZ==(nbSamplesZ-1)) velocitiesZ[indexVelocitiesZFront] =0.0; 
                }
            }
        }
    }
}


// Builds an object to visualize the velocity vectors interpolated at cell centers
// The object will then be accessed and updated each frame (advectVelocities(...))
void Simulation::buildVelocitiesCenters()
{
    //std::cout<<"Building simulation : cell centers velocities."<<std::endl;
    
    objectVelocitiesCenters->setNbVertices(nbSamples*2);
    objectVelocitiesCenters->setNbIndices(objectVelocitiesCenters->getNbVertices());
    
    float vertices[objectVelocitiesCenters->getNbVertices()*4];
    float colors[objectVelocitiesCenters->getNbVertices()*4];
    unsigned int indices[objectVelocitiesCenters->getNbIndices()];
    
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                unsigned int index=iSamples*2+0;
                // Position of vectors beginning
                vertices[index*4+0]=samples[iSamples*4+0];
                vertices[index*4+1]=samples[iSamples*4+1];
                vertices[index*4+2]=samples[iSamples*4+2];
                vertices[index*4+3]=1.0;
                colors[index*4+0]=1.0;
                colors[index*4+1]=1.0;
                colors[index*4+2]=0.0;
                colors[index*4+3]=1.0;
                indices[index]=index;
                        
                index=iSamples*2+1;
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;

                float velocity[]={0.0, 0.0, 0.0, 0.0};
                // X component is interpolated from left and right
                velocity[0]=(velocitiesX[indexVelocitiesXLeft]
                            +velocitiesX[indexVelocitiesXRight])/2.0;
                // Y component is interpolated from bottom and top           
                velocity[1]=(velocitiesY[indexVelocitiesYBottom]
                            +velocitiesY[indexVelocitiesYTop])/2.0;
                // Z component is interpolated from back and front           
                velocity[2]=(velocitiesZ[indexVelocitiesZBack]
                            +velocitiesZ[indexVelocitiesZFront])/2.0;                
                            
                // The vector is printed inversed and transparent at the end (for effect)
                vertices[index*4+0]=samples[iSamples*4+0]-velocity[0]*vectorScale;
                vertices[index*4+1]=samples[iSamples*4+1]-velocity[1]*vectorScale;
                vertices[index*4+2]=samples[iSamples*4+2]-velocity[2]*vectorScale;
                vertices[index*4+3]=1.0;
                colors[index*4+0]=0.0;
                colors[index*4+1]=0.0;
                colors[index*4+2]=0.0;
                colors[index*4+3]=0.0;
                indices[index]=index;
            }  
        }
    }
    // Sends the data into buffers on the GPU
    objectVelocitiesCenters->sendPrimitives(vertices, indices, true);
    objectVelocitiesCenters->sendColors(colors);
}


// Builds an object to visualize the particles 
// The object will then be accessed and updated each frame (advectParticles(...))
void Simulation::buildParticles()
{
    //std::cout<<"Building simulation : particles."<<std::endl;
    
    objectParticles->setNbVertices(nbParticles);
    objectParticles->setNbIndices(nbParticles);

    // Sends the data into buffers on the GPU
    objectParticles->sendPrimitives(particles, particleIndices.data(), true, true);
    objectParticles->sendColors(particleColors, true);
}


// Builds an object to visualize the particles interpolated velocities
// The object will then be accessed and updated each frame (advectParticles(...))
void Simulation::buildParticleVelocities()
{
    //std::cout<<"Building simulation : particles velocities."<<std::endl;
    
    objectParticleVelocities->setNbVertices(nbParticles*2);
    objectParticleVelocities->setNbIndices(objectParticleVelocities->getNbVertices());
    
	float colors[objectParticleVelocities->getNbVertices()*2*4];
	unsigned int indices[objectParticleVelocities->getNbIndices()];
	
	for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
	{
	    indices[iParticles*2+0]=iParticles*2+0;
	    indices[iParticles*2+1]=iParticles*2+1;

        // Trilinear inteprolation of velocity
        float particleVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), particleVelocity);

	    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        {
	        // Vector start on particle location
	        particleVelocities[(iParticles*2+0)*4+iCoord]=particles[iParticles*4+iCoord];
	        // Vector starts with the particle interpolated color
	        colors[(iParticles*2+0)*4+iCoord]=particleColors[iParticles*4+iCoord];
            if (iCoord==3) colors[(iParticles*2+0)*4+iCoord]=0.1;

            // The vector is printed inversed and transparent at the end (for effect)
	        particleVelocities[(iParticles*2+1)*4+iCoord]=particles[iParticles*4+iCoord]-(particleVelocity[iCoord]*vectorScale);
	        colors[(iParticles*2+1)*4+iCoord]=0.0;
	    }
	}
    // Sends the data into buffers on the GPU
    objectParticleVelocities->sendPrimitives(particleVelocities, indices, true, true);
    objectParticleVelocities->sendColors(colors, true);
}


// Builds an object to visualize the solids borders of cells as faces
void Simulation::buildSolids()
{
    //std::cout<<"Building simulation : solids."<<std::endl;

    unsigned int nbVerticesMax=(nbSamplesX+1)*(nbSamplesY+1)*(nbSamplesZ+1)*4;
    unsigned int nbIndicesMax=(nbBordersX+nbBordersY+nbBordersZ)*3*2;

    float verticesTmp[nbVerticesMax*4];
    float normalsTmp[nbVerticesMax*3];
    unsigned int indicesTmp[nbIndicesMax];
    
    unsigned int iVertices=0;
    unsigned int iIndices=0;
    for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
        {
            for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
            {  
                int iXint=(int)iX;
                int iYint=(int)iY;
                int iZint=(int)iZ;
                unsigned int leftType  =type(iXint-1, iYint,   iZint  );
                unsigned int rightType =type(iXint,   iYint,   iZint  );
                unsigned int bottomType=type(iXint,   iYint-1, iZint  );
                unsigned int topType   =type(iXint,   iYint,   iZint  );
                unsigned int backType  =type(iXint,   iYint,   iZint-1);
                unsigned int frontType =type(iXint,   iYint,   iZint  );                

                // Builds x oriented faces
                if ((iY<nbSamplesY) && (iZ<nbSamplesZ))
                {
                    bool solidBorderRight=(leftType!=2) && (rightType==2);
                    bool solidBorderLeft =(leftType==2) && (rightType!=2);

                    if ((solidBorderRight) || (solidBorderLeft))
                    {
                        verticesTmp[(iVertices+0)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+0)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+0)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+0)*4+3]=1.0;

                        verticesTmp[(iVertices+1)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+1)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+1)*4+2]=offset[2]+(iZ+1)*h;
                        verticesTmp[(iVertices+1)*4+3]=1.0;                      

                        verticesTmp[(iVertices+2)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+2)*4+1]=offset[1]+(iY+1)*h;
                        verticesTmp[(iVertices+2)*4+2]=offset[2]+(iZ+1)*h;
                        verticesTmp[(iVertices+2)*4+3]=1.0;

                        verticesTmp[(iVertices+3)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+3)*4+1]=offset[1]+(iY+1)*h;
                        verticesTmp[(iVertices+3)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+3)*4+3]=1.0; 

                        for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                            for (unsigned int iCoord=0 ; iCoord<3 ; ++iCoord)
                                normalsTmp[(iVertices+iCorners)*3+iCoord]=0.0;

                        if (solidBorderRight)
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+1;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+3;
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+0]=-1.0;
                        }
                        else
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+3;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+1;
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+0]=1.0; 
                        }
                        iVertices+=4;
                    }                  
                }
                // Builds y oriented faces
                if ((iX<nbSamplesX) && (iZ<nbSamplesZ))
                {
                    bool solidBorderTop   =(bottomType!=2) && (topType==2);
                    bool solidBorderBottom=(bottomType==2) && (topType!=2);

                    if ( solidBorderTop || solidBorderBottom )
                    {
                        verticesTmp[(iVertices+0)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+0)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+0)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+0)*4+3]=1.0;

                        verticesTmp[(iVertices+1)*4+0]=offset[0]+(iX+1)*h;
                        verticesTmp[(iVertices+1)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+1)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+1)*4+3]=1.0;                      

                        verticesTmp[(iVertices+2)*4+0]=offset[0]+(iX+1)*h;
                        verticesTmp[(iVertices+2)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+2)*4+2]=offset[2]+(iZ+1)*h;
                        verticesTmp[(iVertices+2)*4+3]=1.0;

                        verticesTmp[(iVertices+3)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+3)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+3)*4+2]=offset[2]+(iZ+1)*h;
                        verticesTmp[(iVertices+3)*4+3]=1.0;

                        for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                            for (unsigned int iCoord=0 ; iCoord<3 ; ++iCoord)
                                normalsTmp[(iVertices+iCorners)*3+iCoord]=0.0;

                        if (solidBorderTop)
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+1;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+3;
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+1]=-1.0;  
                        }
                        else
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+3;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+1; 
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+1]=1.0; 
                        }
                        iVertices+=4;
                    }
                }
                // Builds z oriented faces
                if ((iX<nbSamplesX) && (iY<nbSamplesY))
                {
                    bool solidBorderFront=(backType!=2) && (frontType==2);
                    bool solidBorderBack =(backType==2) && (frontType!=2);

                    if ( solidBorderFront || solidBorderBack )
                    {
                        verticesTmp[(iVertices+0)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+0)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+0)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+0)*4+3]=1.0;

                        verticesTmp[(iVertices+1)*4+0]=offset[0]+(iX+1)*h;
                        verticesTmp[(iVertices+1)*4+1]=offset[1]+(iY+0)*h;
                        verticesTmp[(iVertices+1)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+1)*4+3]=1.0;                      

                        verticesTmp[(iVertices+2)*4+0]=offset[0]+(iX+1)*h;
                        verticesTmp[(iVertices+2)*4+1]=offset[1]+(iY+1)*h;
                        verticesTmp[(iVertices+2)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+2)*4+3]=1.0;

                        verticesTmp[(iVertices+3)*4+0]=offset[0]+(iX+0)*h;
                        verticesTmp[(iVertices+3)*4+1]=offset[1]+(iY+1)*h;
                        verticesTmp[(iVertices+3)*4+2]=offset[2]+(iZ+0)*h;
                        verticesTmp[(iVertices+3)*4+3]=1.0; 

                        for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                            for (unsigned int iCoord=0 ; iCoord<3 ; ++iCoord)
                                normalsTmp[(iVertices+iCorners)*3+iCoord]=0.0;

                        if (solidBorderBack)
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+1;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+3;
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+2]=1.0;  
                        }
                        else
                        {
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+3;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+0;
                            indicesTmp[iIndices++]=iVertices+2;
                            indicesTmp[iIndices++]=iVertices+1;
                            for (unsigned int iCorners=0 ; iCorners<4 ; ++iCorners)
                                normalsTmp[(iVertices+iCorners)*3+2]=-1.0;   
                        }
                        iVertices+=4;
                    } 
                }
            }
        }
    }
    objectSolids->setNbVertices(iVertices);
    objectSolids->setNbIndices(iIndices);

    float vertices[objectSolids->getNbVertices()*4];
    float normals[objectSolids->getNbVertices()*3];
    unsigned int indices[objectSolids->getNbIndices()];

    for (iVertices=0 ; iVertices<objectSolids->getNbVertices() ; ++iVertices)
        for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        {
            vertices[iVertices*4+iCoord]=verticesTmp[iVertices*4+iCoord];
            if (iCoord<3) 
                normals[iVertices*3+iCoord]=normalsTmp[iVertices*3+iCoord];
        }
    
    for (iIndices=0 ; iIndices<objectSolids->getNbIndices() ; ++iIndices)
        indices[iIndices]=indicesTmp[iIndices];

    // Sends the data into buffers on the GPU
    objectSolids->sendPrimitives(vertices, indices);
    objectSolids->sendNormals(normals);
}


//______________________________________________________________________________
// Interpolation


// Interpolates trilinearly a data field sampled on cell centers
void Simulation::interpolateFromCenters(const float * const data, const float * const position, float * const result)
{
    float x=(position[0]-offset[0])/h; // x position scaled between 0 and nbSamplesX
    float y=(position[1]-offset[1])/h; // y position scaled between 0 and nbSamplesY
    float z=(position[2]-offset[2])/h; // z position scaled between 0 and nbSamplesZ
    if ((x<0.0) || (x>(float)nbSamplesX) 
     || (y<0.0) || (y>(float)nbSamplesY)
     || (z<0.0) || (z>(float)nbSamplesZ))
        return;
        
    int nbX=(int)nbSamplesX;
    int nbY=(int)nbSamplesY;
    int nbZ=(int)nbSamplesZ;

    // Left-bottom-back closest sample
    float xLeftSample  =x-0.5;
    float yBottomSample=y-0.5;
    float zBackSample  =z-0.5;   
    
    // Distances to left-bottom-back closest sample
    float toLeftSample  =xLeftSample  -floor(xLeftSample  );
    float toBottomSample=yBottomSample-floor(yBottomSample);
    float toBackSample  =zBackSample  -floor(zBackSample  );   
    
    // Corresponding indices on x/y/z axis
    int iLeftSample=int(floor(xLeftSample));
    int iBottomSample=int(floor(yBottomSample));
    int iBackSample=int(floor(zBackSample));
    
    // Height neighbours indices
    int x0y0z0= iBackSample   *(nbY*nbX)+ iBottomSample   *nbX+iLeftSample  ;
    int x1y0z0= iBackSample   *(nbY*nbX)+ iBottomSample   *nbX+iLeftSample+1;
    int x0y1z0= iBackSample   *(nbY*nbX)+(iBottomSample+1)*nbX+iLeftSample  ;
    int x1y1z0= iBackSample   *(nbY*nbX)+(iBottomSample+1)*nbX+iLeftSample+1;
    int x0y0z1=(iBackSample+1)*(nbY*nbX)+ iBottomSample   *nbX+iLeftSample  ;
    int x1y0z1=(iBackSample+1)*(nbY*nbX)+ iBottomSample   *nbX+iLeftSample+1;
    int x0y1z1=(iBackSample+1)*(nbY*nbX)+(iBottomSample+1)*nbX+iLeftSample  ;
    int x1y1z1=(iBackSample+1)*(nbY*nbX)+(iBottomSample+1)*nbX+iLeftSample+1;

    if (iLeftSample  < 0)   {x0y0z0=x1y0z0; x0y1z0=x1y1z0; x0y0z1=x1y0z1; x0y1z1=x1y1z1;}
    if (iLeftSample  >=nbX) {x1y0z0=x0y0z0; x1y1z0=x0y1z0; x1y0z1=x0y0z1; x1y1z1=x0y1z1;}
    if (iBottomSample< 0)   {x0y0z0=x0y1z0; x1y0z0=x1y1z0; x0y0z1=x0y1z1; x1y0z1=x1y1z1;}
    if (iBottomSample>=nbY) {x0y1z0=x0y0z0; x1y1z0=x1y0z0; x0y1z1=x0y0z1; x1y1z1=x1y0z1;}
    if (iBackSample  < 0)   {x0y0z0=x0y1z1; x1y1z0=x1y1z1; x0y0z0=x0y1z1; x1y1z0=x1y1z1;}
    if (iBackSample  >=nbZ) {x0y0z1=x0y0z0; x1y1z1=x1y0z0; x0y0z1=x0y0z0; x1y1z1=x1y0z0;}

    // Trilinear interpolation (4 components) from the height neighbours
    triLinearInterpolation(4, &(data[x0y0z0*4]), &(data[x1y0z0*4]), 
                              &(data[x0y1z0*4]), &(data[x1y1z0*4]), 
                              &(data[x0y0z1*4]), &(data[x1y0z1*4]), 
                              &(data[x0y1z1*4]), &(data[x1y1z1*4]), 
                           toLeftSample, toBottomSample, toBackSample, result);
}


// Interpolates trilinearly a data field sampled on cell borders per components
void Simulation::interpolateFromBorders(const float * const dataX, const float * const dataY, const float * const dataZ, const float * const position, float * const result)
{
    float x=(position[0]-offset[0])/h; // x position scaled between 0 and nbSamplesX
    float y=(position[1]-offset[1])/h; // y position scaled between 0 and nbSamplesY
    float z=(position[2]-offset[2])/h; // z position scaled between 0 and nbSamplesZ
    if ((x<0.0) || (x>(float)nbSamplesX) 
     || (y<0.0) || (y>(float)nbSamplesY)
     || (z<0.0) || (z>(float)nbSamplesZ))
        return;

    int nbX=(int)nbSamplesX;
    int nbY=(int)nbSamplesY;
    int nbZ=(int)nbSamplesZ;
    float xLeftBorder  =x;
    float yBottomBorder=y-0.5;
    float zBackBorder  =z-0.5;

    // Distances to left-bottom-back closest sample
    float toLeftBorder  =xLeftBorder  -floor(xLeftBorder  );
    float toBottomBorder=yBottomBorder-floor(yBottomBorder);
    float toBackBorder  =zBackBorder  -floor(zBackBorder  );
    
    // Corresponding indices on x/y/z axis
    int iLeftBorder  =int(floor(xLeftBorder  ));
    int iBottomBorder=int(floor(yBottomBorder));
    int iBackBorder  =int(floor(zBackBorder  ));     

    // Corresponding global index in dataX
    int iDataXLeft=iBackBorder*(nbY*(nbX+1))+iBottomBorder*(nbX+1)+iLeftBorder;
    
    // Height dataX neighbours indices
    int x0y0z0=iDataXLeft                        ;
    int x1y0z0=iDataXLeft                      +1;
    int x0y1z0=iDataXLeft              +(nbX+1)  ;
    int x1y1z0=iDataXLeft              +(nbX+1)+1;
    int x0y0z1=iDataXLeft+(nbY*(nbX+1))          ;
    int x1y0z1=iDataXLeft+(nbY*(nbX+1))        +1;
    int x0y1z1=iDataXLeft+(nbY*(nbX+1))+(nbX+1)  ;
    int x1y1z1=iDataXLeft+(nbY*(nbX+1))+(nbX+1)+1;

    if (iBottomBorder<0)        {x0y0z0=x0y1z0; x1y0z0=x1y1z0; x0y0z1=x0y1z1; x1y0z1=x1y1z1;}
    if (iBottomBorder>=(nbY-1)) {x0y1z0=x0y0z0; x1y1z0=x1y0z0; x0y1z1=x0y0z1; x1y1z1=x1y0z1;}
    if (iBackBorder  < 0)       {x0y0z0=x0y0z1; x1y0z0=x1y0z1; x0y1z0=x0y1z1; x1y1z0=x1y1z1;}
    if (iBackBorder  >=(nbZ-1)) {x0y0z1=x0y0z0; x1y0z1=x1y0z0; x0y1z1=x0y1z0; x1y1z1=x1y1z0;}

    triLinearInterpolation(1, &(dataX[x0y0z0]), &(dataX[x1y0z0]), 
                              &(dataX[x0y1z0]), &(dataX[x1y1z0]), 
                              &(dataX[x0y0z1]), &(dataX[x1y0z1]), 
                              &(dataX[x0y1z1]), &(dataX[x1y1z1]), 
                           toLeftBorder, toBottomBorder, toBackBorder, &(result[0]));

    xLeftBorder  =x-0.5;
    yBottomBorder=y;
    zBackBorder  =z-0.5;

    // Distances to left-bottom-back closest sample 
    toLeftBorder  =xLeftBorder  -floor(xLeftBorder  );
    toBottomBorder=yBottomBorder-floor(yBottomBorder);
    toBackBorder  =zBackBorder  -floor(zBackBorder  );

    // Corresponding indices on x/y/z axis
    iLeftBorder  =int(floor(xLeftBorder  )); 
    iBottomBorder=int(floor(yBottomBorder));
    iBackBorder  =int(floor(zBackBorder  ));
    
    // Corresponding global index in dataY
    int iDataYBottom=iBackBorder*((nbY+1)*nbX)+iBottomBorder*nbX+iLeftBorder;

    // Height dataY neighbours indices
    x0y0z0=iDataYBottom                    ;
    x1y0z0=iDataYBottom                  +1;
    x0y1z0=iDataYBottom              +nbX  ;
    x1y1z0=iDataYBottom              +nbX+1;
    x0y0z1=iDataYBottom+((nbY+1)*nbX)      ;
    x1y0z1=iDataYBottom+((nbY+1)*nbX)    +1;
    x0y1z1=iDataYBottom+((nbY+1)*nbX)+nbX  ;
    x1y1z1=iDataYBottom+((nbY+1)*nbX)+nbX+1;

    if (xLeftBorder< 0)       {x0y0z0=x1y0z0; x0y1z0=x1y1z0; x0y0z1=x1y0z1; x0y1z1=x1y1z1;}
    if (xLeftBorder>=(nbX-1)) {x1y0z0=x0y0z0; x1y1z0=x0y1z0; x1y0z1=x0y0z1; x1y1z1=x0y1z1;}
    if (zBackBorder< 0)       {x0y0z0=x0y0z1; x1y0z0=x1y0z1; x0y1z0=x0y1z1; x1y1z0=x1y1z1;}
    if (zBackBorder>=(nbZ-1)) {x0y0z1=x0y0z0; x1y0z1=x1y0z0; x0y1z1=x0y1z0; x1y1z1=x1y1z0;}

    triLinearInterpolation(1, &(dataY[x0y0z0]), &(dataY[x1y0z0]), 
                              &(dataY[x0y1z0]), &(dataY[x1y1z0]), 
                              &(dataY[x0y0z1]), &(dataY[x1y0z1]), 
                              &(dataY[x0y1z1]), &(dataY[x1y1z1]), 
                           toLeftBorder, toBottomBorder, toBackBorder, &(result[1]));

    xLeftBorder  =x-0.5;
    yBottomBorder=y-0.5;
    zBackBorder  =z;

    // Distances to left-bottom-back closest sample 
    toLeftBorder  =xLeftBorder  -floor(xLeftBorder  );
    toBottomBorder=yBottomBorder-floor(yBottomBorder);
    toBackBorder  =zBackBorder  -floor(zBackBorder  );

    // Corresponding indices on x/y/z axis
    iLeftBorder  =int(floor(xLeftBorder  )); 
    iBottomBorder=int(floor(yBottomBorder));
    iBackBorder  =int(floor(zBackBorder  ));
    
    // Corresponding global index in dataY
    int iDataZBack=iBackBorder*(nbY*nbX)+iBottomBorder*nbX+iLeftBorder;

    // Height dataZ neighbours indices
    x0y0z0=iDataZBack                ;
    x1y0z0=iDataZBack              +1;
    x0y1z0=iDataZBack          +nbX  ;
    x1y1z0=iDataZBack          +nbX+1;
    x0y0z1=iDataZBack+(nbY*nbX)      ;
    x1y0z1=iDataZBack+(nbY*nbX)    +1;
    x0y1z1=iDataZBack+(nbY*nbX)+nbX  ;
    x1y1z1=iDataZBack+(nbY*nbX)+nbX+1;

    if (xLeftBorder  < 0)       {x0y0z0=x1y0z0; x0y1z0=x1y1z0; x0y0z1=x1y0z1; x0y1z1=x1y1z1;}
    if (xLeftBorder  >=(nbX-1)) {x1y0z0=x0y0z0; x1y1z0=x0y1z0; x1y0z1=x0y0z1; x1y1z1=x0y1z1;}
    if (yBottomBorder< 0)       {x0y0z0=x0y1z0; x1y0z0=x1y1z0; x0y0z1=x0y1z1; x1y0z1=x1y1z1;}
    if (yBottomBorder>=(nbY-1)) {x0y1z0=x0y0z0; x1y1z0=x1y0z0; x0y1z1=x0y0z1; x1y1z1=x1y0z1;}

    triLinearInterpolation(1, &(dataZ[x0y0z0]), &(dataZ[x1y0z0]), 
                              &(dataZ[x0y1z0]), &(dataZ[x1y1z0]), 
                              &(dataZ[x0y0z1]), &(dataZ[x1y0z1]), 
                              &(dataZ[x0y1z1]), &(dataZ[x1y1z1]), 
                           toLeftBorder, toBottomBorder, toBackBorder, &(result[2]));
}


// Simple and intuitive (but inaccurate) time integration method (Forward Euler order 1)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::forwardEuler1stOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue)
{
    // New value = old value + dt * speed of variation
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        resultValue[iCoord]=knownValue[iCoord]+direction*dt*knownVariation[iCoord];
}


// Sufficienty accurate time integration method (Runge-Kutta order 2)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::rungeKutta2ndOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue)
{
    // New value = old value + dt * speed of variation
    // Only this time "variation" is approximated at what it would be at t+dt/2
    float intermediateValue[4];
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        intermediateValue[iCoord]=knownValue[iCoord]+direction*(1.0/2.0)*dt*knownVariation[iCoord];

    float intermediateVariation[]={0.0, 0.0, 0.0, 0.0};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue, intermediateVariation);

    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        resultValue[iCoord]=knownValue[iCoord]+direction*dt*intermediateVariation[iCoord];
}


// Sophisticated (best trade-of accuracy/cost) time integration method (Runge-Kutta order 3)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::rungeKutta3rdOdrTimeIntegration(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue)
{
    // Variation is approximated at several instants between t and t+dt
    // Final value uses these successive variations instead of only the starting one as in Forward Euler
    float intermediateValue1[4];
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        intermediateValue1[iCoord]=knownValue[iCoord]+direction*(1.0/2.0)*dt*knownVariation[iCoord];

    float intermediateVariation1[]={0.0, 0.0, 0.0, 0.0};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue1, intermediateVariation1);

    float intermediateValue2[4];
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        intermediateValue2[iCoord]=knownValue[iCoord]+direction*(3.0/4.0)*dt*intermediateVariation1[iCoord];        
    float intermediateVariation2[]={0.0, 0.0, 0.0, 0.0};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue2, intermediateVariation2);
    
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        resultValue[iCoord]=knownValue[iCoord]
                           +direction*(2.0/9.0)*dt*knownVariation[iCoord]
                           +direction*(3.0/9.0)*dt*intermediateVariation1[iCoord]
                           +direction*(4.0/9.0)*dt*intermediateVariation2[iCoord];
}


// Interpolation scheme choice for the whole simulation
// if direction is -1 : evaluates value is value at t-dt
void Simulation::integrate(float direction, const float * const knownValue, const float * const knownVariation, float dt, float * const resultValue)
{
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        resultValue[iCoord]=knownValue[iCoord];

    //forwardEuler1stOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
    //rungeKutta2ndOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
    rungeKutta3rdOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
}


// Returns type flag for neighbours with provided coordinates (works for boundaries as well)
unsigned int Simulation::type(int iX, int iY, int iZ)
{
    if ( (iX<0) || (iX>=(int)nbSamplesX)
      || (iY<0) || (iY>=(int)nbSamplesY) 
      || (iZ<0) || (iZ>=(int)nbSamplesZ))
    {
        if (solidWalls) return 2;
        else return 1;
    }
    return types[iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX];
}


// Sets integer indices of the cell position pos overlaps
void Simulation::getCell(const float * const pos, int * const iX, int * const iY, int * const iZ)
{
    (*iX)=int(floor((pos[0]-offset[0])/h));
    (*iY)=int(floor((pos[1]-offset[1])/h));
    (*iZ)=int(floor((pos[2]-offset[2])/h));
    if ((*iX)<0) (*iX)=0;
    if ((*iY)<0) (*iY)=0;
    if ((*iZ)<0) (*iZ)=0;  
}


//______________________________________________________________________________
// System resolution


// Evaluates Modified Incomplete Choleski (level 0) preconditioner
void Simulation::MICPreconditioner(const double * const Adiag, const double * const Aright, const double * const Atop, const double * const Afront, Eigen::VectorXd  & precon)
{
    double tau=0.97;
    double safetyConstant=0.25;

    for (int iZ=0 ; iZ<(int)nbSamplesZ ; ++iZ)
    {
        for (int iY=0 ; iY<(int)nbSamplesY ; ++iY)
        {
            for (int iX=0 ; iX<(int)nbSamplesX ; ++iX)
            {
                int iS=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                int iSLeft  = iZ   *(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+(iX-1);
                int iSBottom= iZ   *(nbSamplesY*nbSamplesX)+(iY-1)*nbSamplesX+ iX;
                int iSBack  =(iZ-1)*(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+ iX;

                double e=Adiag[iS];
                if (iX>0)
                {
                    e-=pow(Aright[iSLeft]*(precon)[iSLeft], 2)
                      +tau*(Aright[iSLeft]*Atop[iSLeft]*Afront[iSLeft]
                      *pow((precon)[iSLeft], 2));
                }
                if (iY>0)
                {                  
                    e-=pow(Atop[iSBottom]*(precon)[iSBottom], 2)
                      +tau*(Atop[iSBottom]*Aright[iSBottom]*Afront[iSBottom]
                      *pow((precon)[iSBottom], 2));
                }
                if (iZ>0)
                {                  
                    e-=pow(Afront[iSBack]*(precon)[iSBack], 2)
                      +tau*(Atop[iSBack]*Aright[iSBack]*Afront[iSBack]
                      *pow((precon)[iSBack], 2));
                }                
                if (e<(safetyConstant*Adiag[iS])) e=Adiag[iS];
                if (e!=0.0) precon[iS]=1.0/sqrt(e);
            }
        }
    }
}
    

// Multiplies preconditionneur to r and stores it in z
void Simulation::applyPreconditioner(const double * const Aright, const double * const Atop, const double * const Afront, const Eigen::VectorXd & precon, const Eigen::VectorXd & r, Eigen::VectorXd & z)
{
    double t=0.0;
    Eigen::VectorXd q=Eigen::VectorXd::Zero(nbSamples);

    for (int iZ=0 ; iZ<(int)nbSamplesZ ; ++iZ)
    {
        for (int iY=0 ; iY<(int)nbSamplesY ; ++iY)
        {
            for (int iX=0 ; iX<(int)nbSamplesX ; ++iX)
            {
                int iS=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;

                int iSLeft  = iZ   *(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+(iX-1);
                int iSBottom= iZ   *(nbSamplesY*nbSamplesX)+(iY-1)*nbSamplesX+ iX;
                int iSBack  =(iZ-1)*(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+ iX;

                t=r[iS];
                if (iX>0) t-=Aright[iSLeft  ]*precon[iSLeft]  *q[iSLeft  ];
                if (iY>0) t-=Atop  [iSBottom]*precon[iSBottom]*q[iSBottom];
                if (iZ>0) t-=Afront[iSBack  ]*precon[iSBack]  *q[iSBack  ];
                q[iS]=t*precon[iS];
            }
        }
    }
    for (int iZ=((int)nbSamplesZ-1) ; iZ>=0; --iZ)
    {
        for (int iY=((int)nbSamplesY-1) ; iY>=0; --iY)
        {
            for (int iX=((int)nbSamplesX-1) ; iX>=0 ; --iX)
            {
                int iS=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;

                int iSRight= iZ   *(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+(iX+1);
                int iSTop  = iZ   *(nbSamplesY*nbSamplesX)+(iY+1)*nbSamplesX+ iX;
                int iSFront=(iZ+1)*(nbSamplesY*nbSamplesX)+ iY   *nbSamplesX+ iX;

                t=q[iS];
                if (iX<((int)nbSamplesX-1)) t-=Aright[iS]*precon[iS]*z[iSRight];
                if (iY<((int)nbSamplesY-1)) t-=Atop  [iS]*precon[iS]*z[iSTop  ];
                if (iZ<((int)nbSamplesZ-1)) t-=Afront[iS]*precon[iS]*z[iSFront];
                z[iS]=t*precon[iS];
            }
        }
    }
}


// Multiplies sparse matrix to v and stores to result
void Simulation::multiplySparseMatrix(const double * const Adiag, const double * const Aright, const double * const Atop, const double * const Afront, const Eigen::VectorXd & v, Eigen::VectorXd & result)
{
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                result[iSamples]=Adiag[iSamples]*v[iSamples];
                
                if (iX<(nbSamplesX-1)) // Right
                    result[iSamples]+=Aright[iSamples]*v[iSamples+1         ];
                if (iY<(nbSamplesY-1)) // Top
                    result[iSamples]+=Atop  [iSamples]*v[iSamples+nbSamplesX];
                if (iZ<(nbSamplesZ-1)) // Front
                    result[iSamples]+=Afront[iSamples]*v[iSamples+(nbSamplesY*nbSamplesX)];


                if (iX>0) // Left
                    result[iSamples]+=Aright[iSamples-1]         *v[iSamples-1         ];
                if (iY>0) // Bottom
                    result[iSamples]+=Atop  [iSamples-nbSamplesX]*v[iSamples-nbSamplesX];
                if (iZ>0) // Back
                    result[iSamples]+=Afront[iSamples-(nbSamplesY*nbSamplesX)]*v[iSamples-(nbSamplesY*nbSamplesX)];                  
            }
        }
    }
}


// Fills pressures from sparse matrix and rhs=b
void Simulation::conjugateGradient(const double * const Adiag, const double * const Aright, const double * const Atop, const double * const Afront, const double * const b)
{
    unsigned int maxIterations=100;
    double tol=1e-6;
    
    // Initializes pressures to null
    Eigen::VectorXd p=Eigen::VectorXd::Zero(nbSamples);
    // Copy r from b and return if r null
    Eigen::VectorXd r(nbSamples);

    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        r[iSamples]=b[iSamples];
        
    // rInfinityNorm is the maximum of absolute values of r
    double rInfinityNorm=r.lpNorm<Eigen::Infinity>(); // maximum of absolute values
    if (rInfinityNorm<=tol)
    {
        for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
            pressures[iSamples]=(float)p[iSamples];
        std::cout<<"Stability reached."<<std::endl;
        return;
    }

    // Builds Preconditionner
    Eigen::VectorXd precon=Eigen::VectorXd::Zero(nbSamples);
    MICPreconditioner(Adiag, Aright, Atop, Afront, precon);

    // Apply preconditionner from r to z
    Eigen::VectorXd z=Eigen::VectorXd::Zero(nbSamples);
    applyPreconditioner(Aright, Atop, Afront, precon, r, z);
    // Copy z to s
    Eigen::VectorXd s(z);
    double sigma=z.dot(r);
    unsigned int iLast=maxIterations;
    for (unsigned int i=0 ; i<maxIterations ; ++i)
    {
        // Multiply matrix A to s to get z;
        multiplySparseMatrix(Adiag, Aright, Atop, Afront, s, z);
        double alpha=0.0;
        double div=z.dot(s);
        if (div!=0.0) alpha=sigma/div;
        p+=alpha*s;
        r-=alpha*z;
        rInfinityNorm=r.lpNorm<Eigen::Infinity>();
        if (rInfinityNorm<=tol)
        {
            iLast=i;
            i=maxIterations;
        }
        else
        {
            applyPreconditioner(Aright, Atop, Afront, precon, r, z);
            double beta=0.0;
            double newSigma=z.dot(r);
            if (sigma!=0.0) beta=newSigma/sigma;
            s=z+beta*s;
            sigma=newSigma;
        }
    }
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        pressures[iSamples]=(float)p[iSamples];

    //std::cout<<"Resolution stopped after "<<iLast<<" iterations";
    //std::cout<<"(max div = "<<rInfinityNorm<<")"<<std::endl;
}


// Fills a table of every cell divergence
void Simulation::setDivergences(double * const divergences)
{
    double scale=1.0/h;
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {  
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;

                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;

                // Sets Right Hand Side        
                divergences[iSamples]=scale*( velocitiesX[indexVelocitiesXRight]
                                             -velocitiesX[indexVelocitiesXLeft]
                                             +velocitiesY[indexVelocitiesYTop]
                                             -velocitiesY[indexVelocitiesYBottom] 
                                             +velocitiesZ[indexVelocitiesZFront]
                                             -velocitiesZ[indexVelocitiesZBack] );
            }
        }
    }
}


//______________________________________________________________________________
// Main steps


// Modifies the color samples so that color derivative is null (particles moving through the grid keep their color)
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectColors()
{
    float newColors[nbSamples*4];
    
    float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
    float virtualParticlePos[]={0.0, 0.0, 0.0, 1.0};
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                
                // Interpolates the velocity at cell center (sampleVelocity)
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                                  +velocitiesX[indexVelocitiesXRight])/2.0;
                sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                                  +velocitiesY[indexVelocitiesYTop])/2.0;
                sampleVelocity[2]=(velocitiesZ[indexVelocitiesZBack]
                                  +velocitiesZ[indexVelocitiesZFront])/2.0;                
                
                // At t-dt, where was the material now arriving at sample location ?
                integrate(-1.0,                        // direction
                        &(samples[iSamples*4+0]),      // startValue
                          sampleVelocity,              // startVariation 
                          dt,                          // dt
                          virtualParticlePos);         // resultValue

                // What is the color of the material now arriving on sample location ?
                newColors[iSamples*4+0]=colors[iSamples*4+0];
                newColors[iSamples*4+1]=colors[iSamples*4+1];
                newColors[iSamples*4+2]=colors[iSamples*4+2];
                newColors[iSamples*4+3]=colors[iSamples*4+3];
                
                int iXvirtual, iYvirtual, iZvirtual;
                getCell(virtualParticlePos, &iXvirtual, &iYvirtual, &iZvirtual);

                if (type(iXvirtual, iYvirtual, iZvirtual)==0)             
                    interpolateFromCenters(colors, virtualParticlePos, &(newColors[iSamples*4+0]));
            }  
        }
    }
    // Updates the "colors" buffer
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
            colors[iSamples*4+iCoord]=newColors[iSamples*4+iCoord];
    
    // Updates the corresponding data if stored on GPU
    if (objectSamples!=NULL) objectSamples->updateColors(colors, true);
}


// Updates the velocities of centers from known borders
void Simulation::updateVelocitiesFromBorders()
{
    if (objectVelocitiesCenters!=NULL)
    {
        float vertices[nbSamples*4*2];
        for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
            {
                for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
                {
                    unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                    
                    float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
                    
                    // Interpolates the velocity at cell center (sampleVelocity)
                    unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                    unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                    unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                    unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                    unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                    unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                    sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                                      +velocitiesX[indexVelocitiesXRight])/2.0;
                    sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                                      +velocitiesY[indexVelocitiesYTop])/2.0;
                    sampleVelocity[2]=(velocitiesZ[indexVelocitiesZBack]
                                      +velocitiesZ[indexVelocitiesZFront])/2.0;   
                                      
                    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                    {
                        vertices[(iSamples*2+0)*4+iCoord]=samples[iSamples*4+iCoord];
                        vertices[(iSamples*2+1)*4+iCoord]=samples[iSamples*4+iCoord]-sampleVelocity[iCoord]*vectorScale;
                    }
                }
            }
        }
        objectVelocitiesCenters->updateVertices(vertices, true);
    }
    
    if (objectVelocitiesBorders!=NULL)
    {
        float colors[objectVelocitiesBorders->getNbVertices()*4];
        unsigned int iBorders=0;
        for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    	    {
    	        for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
    	        {   
                    float normVelocityX=velocitiesX[iZ*(nbSamplesY*(nbSamplesX+1))+iY*(nbSamplesX+1)+iX];
                    if (normVelocityX<0.0) normVelocityX=-normVelocityX;
                    colors[iBorders*4+0]=normVelocityX*10.0;
                    colors[iBorders*4+1]=0.0; 
                    colors[iBorders*4+2]=0.0; 
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
    	    {
    	        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
    	        {   
                    float normVelocityY=velocitiesY[iZ*((nbSamplesY+1)*nbSamplesX)+iY*nbSamplesX+iX];
                    if (normVelocityY<0.0) normVelocityY=-normVelocityY;
                    colors[iBorders*4+0]=0.0;
                    colors[iBorders*4+1]=normVelocityY*10.0;
                    colors[iBorders*4+2]=0.0; 
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
            {
                for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
                {   
                    float normVelocityZ=velocitiesZ[iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX];
                    if (normVelocityZ<0.0) normVelocityZ=-normVelocityZ;
                    colors[iBorders*4+0]=0.0;
                    colors[iBorders*4+1]=0.0;
                    colors[iBorders*4+2]=normVelocityZ*10.0;
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        objectVelocitiesBorders->updateColors(colors, true);
    }
}


// Updates the velocities of borders from known centers
void Simulation::updateVelocitiesFromCenters(const float * const centeredVel)
{
    // Resets velocities components to 0
    for (unsigned int iBordersX=0 ; iBordersX<nbBordersX ; ++iBordersX)
        velocitiesX[iBordersX]=0.0;
    for (unsigned int iBordersY=0 ; iBordersY<nbBordersY ; ++iBordersY)
        velocitiesY[iBordersY]=0.0;
    for (unsigned int iBordersZ=0 ; iBordersZ<nbBordersZ ; ++iBordersZ)
        velocitiesZ[iBordersZ]=0.0;    
    
    // Uses the new centers velocity to update the separated velocity components
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                
                float coefLeft=0.5;
                float coefRight=0.5;
                float coefBottom=0.5;
                float coefTop=0.5;
                float coefBack=0.5;
                float coefFront=0.5;

                if (iX==0)              coefLeft  =1.0;
                if (iX==(nbSamplesX-1)) coefRight =1.0;
                if (iY==0)              coefBottom=1.0;
                if (iY==(nbSamplesY-1)) coefTop   =1.0;  
                if (iZ==0)              coefBack  =1.0;
                if (iZ==(nbSamplesZ-1)) coefFront =1.0;                            
                
                velocitiesX[indexVelocitiesXLeft]  +=coefLeft  *centeredVel[iSamples*4+0];
                velocitiesX[indexVelocitiesXRight] +=coefRight *centeredVel[iSamples*4+0];
                velocitiesY[indexVelocitiesYBottom]+=coefBottom*centeredVel[iSamples*4+1];
                velocitiesY[indexVelocitiesYTop]   +=coefTop   *centeredVel[iSamples*4+1];
                velocitiesZ[indexVelocitiesZBack]  +=coefBack  *centeredVel[iSamples*4+2];
                velocitiesZ[indexVelocitiesZFront] +=coefFront *centeredVel[iSamples*4+2];
            }  
        }
    }
    enforceVelocitiesBorders();
    
    // Updates the corresponding data if stored on GPU
    if (objectVelocitiesCenters!=NULL)
    {
        float vertices[nbSamples*4*2];
        for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        {
            for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
            {
                vertices[(iSamples*2+0)*4+iCoord]=samples[iSamples*4+iCoord];
                vertices[(iSamples*2+1)*4+iCoord]=samples[iSamples*4+iCoord]-centeredVel[iSamples*4+iCoord]*vectorScale;
            }
        }
        objectVelocitiesCenters->updateVertices(vertices, true);
    }
    
    if (objectVelocitiesBorders!=NULL)
    {
        float colors[objectVelocitiesBorders->getNbVertices()*4];
        unsigned int iBorders=0;
        for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
            {
                for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
                {   
                    float normVelocityX=velocitiesX[iZ*(nbSamplesY*(nbSamplesX+1))+iY*(nbSamplesX+1)+iX];
                    if (normVelocityX<0.0) normVelocityX=-normVelocityX;
                    colors[iBorders*4+0]=normVelocityX*10.0;
                    colors[iBorders*4+1]=0.0; 
                    colors[iBorders*4+2]=0.0; 
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
            {
                for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
                {   
                    float normVelocityY=velocitiesY[iZ*((nbSamplesY+1)*nbSamplesX)+iY*nbSamplesX+iX];
                    if (normVelocityY<0.0) normVelocityY=-normVelocityY;
                    colors[iBorders*4+0]=0.0;
                    colors[iBorders*4+1]=normVelocityY*10.0;
                    colors[iBorders*4+2]=0.0; 
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        for (unsigned int iZ=0 ; iZ<(nbSamplesZ+1) ; ++iZ)
        {
            for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
            {
                for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
                {   
                    float normVelocityZ=velocitiesZ[iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX];
                    if (normVelocityZ<0.0) normVelocityZ=-normVelocityZ;
                    colors[iBorders*4+0]=0.0;
                    colors[iBorders*4+1]=0.0;
                    colors[iBorders*4+2]=normVelocityZ*10.0;
                    colors[iBorders*4+3]=1.0;
                    iBorders++;
                }
            }
        }
        objectVelocitiesBorders->updateColors(colors, true);
    }
}


// Modifies the velocity samples so that velocity derivative is null (particles moving through the grid keep their velocity)
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectVelocities()
{
    float newVelocities[nbSamples*4];
    
    float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
    float virtualParticlePos[]={0.0, 0.0, 0.0, 1.0};

    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {    
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                
                // Interpolates the velocity at cell center (sampleVelocity)
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                                  +velocitiesX[indexVelocitiesXRight])/2.0;
                sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                                  +velocitiesY[indexVelocitiesYTop])/2.0;
                sampleVelocity[2]=(velocitiesZ[indexVelocitiesZBack]
                                  +velocitiesZ[indexVelocitiesZFront])/2.0; 
                
                // At t-dt, where was the material now arriving at sample location ?
                integrate(-1.0,                        // direction
                        &(samples[iSamples*4+0]),      // startValue
                          sampleVelocity,              // startVariation 
                          dt,                          // dt
                          virtualParticlePos);         // resultValue
                
                // What is the velocity of the material now arriving on sample location ?
                newVelocities[iSamples*4+0]=sampleVelocity[0];
                newVelocities[iSamples*4+1]=sampleVelocity[1];
                newVelocities[iSamples*4+2]=sampleVelocity[2];
                newVelocities[iSamples*4+3]=0.0;
                
                int iXvirtual, iYvirtual, iZvirtual;
                getCell(virtualParticlePos, &iXvirtual, &iYvirtual, &iZvirtual);

                if (type(iXvirtual, iYvirtual, iZvirtual)==0)
                    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, virtualParticlePos, &(newVelocities[iSamples*4+0]));
            }  
        }
    }
    updateVelocitiesFromCenters(newVelocities);
}


// Change velocities according to forces
void Simulation::applyForces()
{
    float newVelocities[nbSamples*4];
    
    // Uses the new centers velocity to update the separated velocity components
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {   
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;  
                            
                float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
                
                // Interpolates the velocity at cell center (sampleVelocity)
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                                  +velocitiesX[indexVelocitiesXRight])/2.0;
                sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                                  +velocitiesY[indexVelocitiesYTop])/2.0;
                sampleVelocity[2]=(velocitiesZ[indexVelocitiesZBack]
                                  +velocitiesZ[indexVelocitiesZFront])/2.0; 
                                      
                for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                    newVelocities[iSamples*4+iCoord]=sampleVelocity[iCoord];
                    
                if (type(iX, iY, iZ)==0)
                {        
                    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                        newVelocities[iSamples*4+iCoord]+=dt*forces[iSamples*4+iCoord];
                }
            }
        }
    }       
    updateVelocitiesFromCenters(newVelocities);
}


// Projection : substracting the pressure gradient
// Pressures evaluated so that incompressibility and boundary conditions are enforced 
void Simulation::project()
{
    // To evaluate Right Hand Side of equation
    double scaleRhs=1.0/h;
    double rhs[nbSamples];
    
    // To evaluate A matrix
    double scaleA=dt/(density*h*h);
    double Adiag[nbSamples];
    double Aright[nbSamples];
    double Atop[nbSamples];
    double Afront[nbSamples];
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
    {
        rhs[iSamples]=0.0;
        Adiag[iSamples]=0.0;
        Aright[iSamples]=0.0;
        Atop[iSamples]=0.0;
        Afront[iSamples]=0.0;
    }
    
    double velocitySolids=0.0;
    setDivergences(rhs);

    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                
                unsigned int typeSample=type(int(iX), int(iY), int(iZ));
                unsigned int typeNeighbourLeft  =type(int(iX)-1, int(iY),   int(iZ)  );
                unsigned int typeNeighbourRight =type(int(iX)+1, int(iY),   int(iZ)  );
                unsigned int typeNeighbourBottom=type(int(iX),   int(iY)-1, int(iZ)  );
                unsigned int typeNeighbourTop   =type(int(iX),   int(iY)+1, int(iZ)  );
                unsigned int typeNeighbourBack  =type(int(iX),   int(iY),   int(iZ)-1);
                unsigned int typeNeighbourFront =type(int(iX),   int(iY),   int(iZ)+1);

                
                // If Fluid cell
                if (typeSample==0)
                {
                    // Sets Right Hand Side
                    rhs[iSamples]=-rhs[iSamples];
                
                    if (typeNeighbourLeft==2) // If solid on the left
                        rhs[iSamples]-=scaleRhs*(velocitiesX[indexVelocitiesXLeft]
                                                -velocitySolids); //velocitySolidXLeft
                    if (typeNeighbourRight==2) // If solid on the right
                        rhs[iSamples]+=scaleRhs*(velocitiesX[indexVelocitiesXRight]
                                                -velocitySolids); //velocitySolidXRight
                    if (typeNeighbourBottom==2) // If solid on the bottom
                        rhs[iSamples]-=scaleRhs*(velocitiesY[indexVelocitiesYBottom]
                                                -velocitySolids); //velocitySolidYBottom
                    if (typeNeighbourTop==2) // If solid on the top
                        rhs[iSamples]+=scaleRhs*(velocitiesY[indexVelocitiesYTop]
                                                -velocitySolids); //velocitySolidYTop
                    if (typeNeighbourBack==2) // If solid on the back
                        rhs[iSamples]-=scaleRhs*(velocitiesZ[indexVelocitiesZBack]
                                                -velocitySolids); //velocitySolidZBack
                    if (typeNeighbourFront==2) // If solid on the front
                        rhs[iSamples]+=scaleRhs*(velocitiesZ[indexVelocitiesZFront]
                                                -velocitySolids); //velocitySolidZFront


                    // Sets matrix A compact form
                    if (typeNeighbourLeft!=2)   Adiag[iSamples]+=scaleA;
                    if (typeNeighbourRight!=2)  Adiag[iSamples]+=scaleA;
                    if (typeNeighbourBottom!=2) Adiag[iSamples]+=scaleA;
                    if (typeNeighbourTop!=2)    Adiag[iSamples]+=scaleA;
                    if (typeNeighbourBack!=2)   Adiag[iSamples]+=scaleA;
                    if (typeNeighbourFront!=2)  Adiag[iSamples]+=scaleA;

                        
                    if (typeNeighbourRight==0) Aright[iSamples]=-scaleA;
                    if (typeNeighbourTop==0)   Atop  [iSamples]=-scaleA;
                    if (typeNeighbourFront==0) Afront[iSamples]=-scaleA;
                }
                else
                {
                    rhs[iSamples]=0.0;
                }
            }
        }
    }
    // Resolves system with preconditionned conjugate gradient method
    conjugateGradient(Adiag, Aright, Atop, Afront, rhs);
    
    // Substracts pressure gradient
    float scale=dt/(density*h);
    for (unsigned int iZ=0 ; iZ<nbSamplesZ ; ++iZ)
    {
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
                unsigned int indexVelocitiesXLeft  = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+ iX   ;
                unsigned int indexVelocitiesXRight = iZ  *( nbSamplesY   *(nbSamplesX+1))+ iY   *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom= iZ  *((nbSamplesY+1)*nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesYTop   = iZ  *((nbSamplesY+1)*nbSamplesX    )+(iY+1)* nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZBack  = iZ  *( nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                unsigned int indexVelocitiesZFront =(iZ+1)*(nbSamplesY   *nbSamplesX    )+ iY   * nbSamplesX   + iX   ;
                
                if (types[iSamples]==0)
                {
                    velocitiesX[indexVelocitiesXLeft]  -=scale*pressures[iSamples];
                    velocitiesX[indexVelocitiesXRight] +=scale*pressures[iSamples];
                    velocitiesY[indexVelocitiesYBottom]-=scale*pressures[iSamples];
                    velocitiesY[indexVelocitiesYTop]   +=scale*pressures[iSamples];
                    velocitiesZ[indexVelocitiesZBack]  -=scale*pressures[iSamples];
                    velocitiesZ[indexVelocitiesZFront] +=scale*pressures[iSamples];
                }
            }
        }
    }
    
    enforceVelocitiesBorders();
    updateVelocitiesFromBorders();

    if (objectPressures!=NULL)
    {
        float colors[objectPressures->getNbVertices()*4];
        for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        {
                colors[iSamples*4+0]=-pressures[iSamples]/100.0;
                colors[iSamples*4+1]=pressures[iSamples]/100.0;
                colors[iSamples*4+2]=0.0;
                colors[iSamples*4+3]=1.0;
        }
        objectPressures->updateColors(colors, true);
    }
}


// Moves the particles by interpolating velocity and integrating the new position from it, changes the color as well by interpolating on the color field
void Simulation::advectParticles()
{
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        if (types[iSamples]!=2) types[iSamples]=1;
        
    //float velocityColors[nbParticles*4*2];
    for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
	{
	    int iX, iY, iZ;
        getCell(&(particles[iParticles*4+0]), &iX, &iY, &iZ);
            if (type(iX, iY, iZ)!=2)
                types[iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX]=0;
        
        float knownVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), knownVelocity);

        // What is the new position after dt ?
        integrate(1.0,                            // direction
                &(particles[iParticles*4+0]),     // startValue
                  knownVelocity,                  // startVariation 
                  dt,                             // dt
                &(particles[iParticles*4+0]));    // resultValue 

        
        // What is the velocity at this new position ?
        float newVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), newVelocity);
        for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        {
            particleVelocities[(iParticles*2+0)*4+iCoord]=particles[iParticles*4+iCoord];
            particleVelocities[(iParticles*2+1)*4+iCoord]=particles[iParticles*4+iCoord]-(newVelocity[iCoord]*vectorScale);
        }

        // What is the color at this new position ?
        //interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0]));
        //for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        //{
        //    velocityColors[(iParticles*2+0)*4+iCoord]=particleColors[iParticles*4+iCoord];
	    //    velocityColors[(iParticles*2+1)*4+iCoord]=0.0;
	    //S}

        // Cancelling color of stuck particles
        // Particles coordinates scaled between 0 and nbSamplesX/Y/Z
        float x=(particles[iParticles*4+0]-offset[0])/h;
        float y=(particles[iParticles*4+1]-offset[1])/h;
        float z=(particles[iParticles*4+2]-offset[2])/h;
        if ((x<0.0) || (x>(float)nbSamplesX) 
         || (y<0.0) || (y>(float)nbSamplesY)
         || (z<0.0) || (z>(float)nbSamplesZ))
            for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                particleColors[iParticles*4+iCoord]=0.0;
    }

    // Particles depth sorting (only if particles or particleVelocities visualization)
    if ((objectParticles!=NULL) || (objectParticleVelocities!=NULL))
    {
        const float * const c=camera.getCenter();
        float center[]={c[0], c[1], c[2], 1.0};
        if (this->model!=NULL)
        {
            float modelInv[16];
            getInverseMatrix(this->model, modelInv); 
            multMatrixToPoint(center, modelInv);
        }
        float distances[nbParticles];
        for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
            distances[iParticles]=distance(center, &(particles[iParticles*4+0]), true);

        std::sort(particleIndices.begin(), particleIndices.end(), compareParticlesDepth(distances));
    }

    // Updates the corresponding data if stored on GPU
    if (objectParticles!=NULL)
    {
        objectParticles->updateVertices(particles, true);
        objectParticles->updateIndices(particleIndices.data(), true);
        objectParticles->updateColors(particleColors, true);
    }
    if (objectParticleVelocities!=NULL)
    {
        unsigned int velocitiesIndices[nbParticles*2];
        for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
        {
            velocitiesIndices[iParticles*2+0]=particleIndices[iParticles]*2+0;
            velocitiesIndices[iParticles*2+1]=particleIndices[iParticles]*2+1;
            //std::cout<<" "<<velocitiesIndices[iParticles*2+0];
        }
        objectParticleVelocities->updateVertices(particleVelocities, true);
        objectParticleVelocities->updateIndices(velocitiesIndices, true);
        //objectParticleVelocities->updateColors(velocityColors, true);
    }
    if (objectTypes!=NULL)
    {
        float colorsTypes[objectTypes->getNbVertices()*4];
        for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        {
            colorsTypes[iSamples*4+0]=0.0;
            colorsTypes[iSamples*4+1]=0.0;
            colorsTypes[iSamples*4+2]=0.0;
            colorsTypes[iSamples*4+3]=1.0;
            colorsTypes[iSamples*4+types[iSamples]]=1.0;
        }
        objectTypes->updateColors(colorsTypes, true);
    }
}


//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Constructor
Simulation::Simulation(const Shaders & shaders,
                       Camera & camera,
                       Scene & scene, 
                       const unsigned int defaultShader,
                       const unsigned int spriteShader,
                       const float size,
                       const bool solidWalls,
                       const float density,
                       const float viscosity,
                       const float dt,
                       const unsigned int nbSamplesX,
                       const unsigned int nbSamplesY,
                       const unsigned int nbSamplesZ,
                       const unsigned int nbParticlesCoef)
                       : shaders(shaders),
                         camera(camera),
                         scene(scene), 
                         defaultShader(defaultShader),
                         spriteShader(spriteShader),
                         size(size),
                         solidWalls(solidWalls),
                         density(density),
                         viscosity(viscosity),
                         dt(dt),
                         nbSamplesX(nbSamplesX),
                         nbSamplesY(nbSamplesY),
                         nbSamplesZ(nbSamplesZ),
                         nbParticlesCoef(nbParticlesCoef)
{
        // About the data stored on centers of the MAC grid
    // Total samples number
    this->nbSamples=this->nbSamplesX*this->nbSamplesY*this->nbSamplesZ;

    // Samples positions
    this->samples=new float[this->nbSamples*4];
    // Samples pressures
    this->pressures=new float[this->nbSamples];
    // Samples colors
    this->colors=new float[this->nbSamples*4];
    // Samples forces
    this->forces=new float[this->nbSamples*4];
    // Samples types
    this->types=new unsigned int[this->nbSamples];    
           
           
        // About the data stored on borders of the MAC grid (on faces)
    // Number of left/right cell borders
    this->nbBordersX=(this->nbSamplesX+1)*this->nbSamplesY*this->nbSamplesZ;
    // Number of bottom/top cell borders
    this->nbBordersY=this->nbSamplesX*(this->nbSamplesY+1)*this->nbSamplesZ;
    // Number of back/front cell borders
    this->nbBordersZ=this->nbSamplesX*this->nbSamplesY*(this->nbSamplesZ+1);


    // X velocity components, sampled on left/right borders
    this->velocitiesX=new float[this->nbBordersX];
    // Y velocity components, sampled on bottom/top borders
    this->velocitiesY=new float[this->nbBordersY];
    // Z velocity components, sampled on back/front borders
    this->velocitiesZ=new float[this->nbBordersZ];
    
    
        // About the particles used for visualization (advected on the whole grid)
    // Number of rendered particles (=nbSamples*8^nbParticlesCoef)
    this->nbParticles=this->nbSamples*(unsigned int)pow(8.0, this->nbParticlesCoef);
    // Particle positions
    this->particles=new float[this->nbParticles*4];
    // Grid interpolated particles colors
    this->particleColors=new float[this->nbParticles*4];
    // Grid interpolated particles velocities
    this->particleVelocities=new float[this->nbParticles*4*2];
        
        // Random generator
    this->rand = new MTRand(1995);

        // Geometric parameters
    // Cell side length (square or cubes)
    this->h=this->size/float(this->nbSamplesX);
    // Start point left/bottom/back of the grid
    this->offset[0]=-this->size/2.0;
    this->offset[1]=-this->size*(float(this->nbSamplesY)/float(this->nbSamplesX))/2.0;
    this->offset[2]=-this->size*(float(this->nbSamplesZ)/float(this->nbSamplesX))/2.0;

    // Offset from left/bottom/back of grid cell
    this->offsetInCell=this->h/2.0;
    // Scale coefficient to display vector lengths
    this->vectorScale=1.0/2.0;
    // Model matrix for the simulation
    setToIdentity(this->model);

        // Simulation parameters
    // Vector of acceleration due to gravity
    this->g[0]=0.0; this->g[1]=-9.81; this->g[2]=0.0;


        // Objects for OpenGL visualization
    // To render samples as positions and colors
    this->objectSamples=NULL;  
    // To render pressures as color
    this->objectPressures=NULL;
    // To render forces as vectors and colors
    this->objectForces=NULL;
    // To render types as colors
    this->objectTypes=NULL;
    // To render the grid as lines
    this->objectGrid=NULL;
    // To render velocities interpolated on samples as vectors and colors
    this->objectVelocitiesCenters=NULL;  
    // To render velocities components as colors
    this->objectVelocitiesBorders=NULL; 
    // To render particles as positions and colors
    this->objectParticles=NULL;   
    // To render particles interpolated velocities as vectors and colors
    this->objectParticleVelocities=NULL; 


    // Printing a few informations
    std::cout<<"Grid : ["<<this->nbSamplesX<<"x"<<this->nbSamplesY<<"x"<<this->nbSamplesZ<<"]"<<std::endl;
    std::cout<<"Cells : "<<this->nbSamples<<std::endl;
    std::cout<<"Particles : "<<this->nbParticles<<std::endl;


        // Initialisation of buffers
    // Initialisation of required positions colors and velocities for simulation
    this->initSimulation();
    // Initialisation of particles for visualization
    this->initVisualization();
}


// Cleans memory for Simulation
Simulation::~Simulation()
{
    delete [] samples;
    delete [] pressures;
    delete [] colors;
    delete [] forces;
    delete [] types;
    
    delete [] velocitiesX;
    delete [] velocitiesY;
    delete [] velocitiesZ;
    
    delete [] particles;
    delete [] particleColors;
    delete [] particleVelocities;

    delete rand;
}


//______________________________________________________________________________
// Setters / getters
void Simulation::setModel(const float * const model)
{
    for (unsigned int i=0 ; i<16 ; ++i)
        this->model[i]=model[i];
}


//______________________________________________________________________________
// Drawing set-up


// Creates, builds and add to draw list an object for samples visualization
unsigned int Simulation::drawSamples()
{
    this->objectSamples=new Object(GL_POINTS);
    this->buildSamples();
    unsigned int storedObjectSamples=scene.storeObject(objectSamples);
    unsigned int samplesID=scene.addObjectToDraw(storedObjectSamples);
    scene.setDrawnObjectShader(samplesID, defaultShader);
    return samplesID;
}


// Creates, builds and add to draw list an object for pressures visualization
unsigned int Simulation::drawPressures()
{
    this->objectPressures=new Object(GL_POINTS);
    this->buildPressures();
    unsigned int storedObjectPressures=scene.storeObject(objectPressures);
    unsigned int pressuresID=scene.addObjectToDraw(storedObjectPressures);
    scene.setDrawnObjectShader(pressuresID, defaultShader);
    return pressuresID;
}


// Creates, builds and add to draw list an object for forces visualization
unsigned int Simulation::drawForces()
{
    this->objectForces=new Object(GL_LINES);
    this->buildForces();
    unsigned int storedObjectForces=scene.storeObject(objectForces);
    unsigned int forcesID=scene.addObjectToDraw(storedObjectForces);
    scene.setDrawnObjectShader(forcesID, defaultShader);
    return forcesID;
}


// Creates, builds and add to draw list an object for types visualization
unsigned int Simulation::drawTypes()
{
    this->objectTypes=new Object(GL_POINTS);
    this->buildTypes();
    unsigned int storedObjectTypes=scene.storeObject(objectTypes);
    unsigned int typesID=scene.addObjectToDraw(storedObjectTypes);
    scene.setDrawnObjectShader(typesID, defaultShader);
    return typesID;
}


// Creates, builds and add to draw list an object for grid visualization
unsigned int Simulation::drawGrid()  
{
    this->objectGrid=new Object(GL_LINES);
    this->buildGrid();
    unsigned int storedObjectGrid=scene.storeObject(objectGrid);
    unsigned int gridID=scene.addObjectToDraw(storedObjectGrid);  
    scene.setDrawnObjectShader(gridID, defaultShader);
    float shade=0.6;
    float grey[]={shade, shade, shade, 1.0};
    scene.setDrawnObjectColor(gridID, grey);
    return gridID;
}
   
   
// Creates, builds and add to draw list an object for velocity components visualization
unsigned int Simulation::drawVelocitiesBorders() 
{
    this->objectVelocitiesBorders=new Object(GL_POINTS);
    this->buildVelocitiesBorders();
    unsigned int storedObjectVelocitiesBorders=scene.storeObject(objectVelocitiesBorders);
    unsigned int velocitiesBordersID=scene.addObjectToDraw(storedObjectVelocitiesBorders);
    scene.setDrawnObjectShader(velocitiesBordersID, defaultShader);
    return velocitiesBordersID;
}


// Creates, builds and add to draw list an object for velocity vectors visualization
unsigned int Simulation::drawVelocitiesCenters() 
{
    this->objectVelocitiesCenters=new Object(GL_LINES);
    this->buildVelocitiesCenters();
    unsigned int storedObjectVelocitiesCenters=scene.storeObject(objectVelocitiesCenters);
    unsigned int velocitiesCentersID=scene.addObjectToDraw(storedObjectVelocitiesCenters);
    scene.setDrawnObjectShader(velocitiesCentersID, defaultShader);
    return velocitiesCentersID;
}


// Creates, builds and add to draw list an object for particles visualization
unsigned int Simulation::drawParticles()
{
    this->objectParticles=new Object(GL_POINTS);
    this->buildParticles();
    unsigned int storedObjectParticles=scene.storeObject(objectParticles);
    unsigned int particlesID=scene.addObjectToDraw(storedObjectParticles);
    
    // Particles are drawn with point sprites 
    // (textures aligned on screen for each point)
    // Black parts of sprites are transparent (blending enabled, depth test disabled)
    scene.setDrawnObjectShader(particlesID, spriteShader);
    scene.setDrawnObjectTransparency(particlesID, true);
    unsigned int starTextureID=loadTexture("../textures/star1.ppm");
    scene.setDrawnObjectTextureID(particlesID, 0, starTextureID);
    return particlesID;
}


// Creates, builds and add to draw list an object for particles velocities visualization
unsigned int Simulation::drawParticlesVelocities()
{
    this->objectParticleVelocities=new Object(GL_LINES);
    this->buildParticleVelocities();
    unsigned int storedObjectParticleVelocities=scene.storeObject(objectParticleVelocities);
    unsigned int particleVelocitiesID=scene.addObjectToDraw(storedObjectParticleVelocities);
    scene.setDrawnObjectShader(particleVelocitiesID, defaultShader);
    scene.setDrawnObjectTransparency(particleVelocitiesID, true);
    return particleVelocitiesID;
}

// Creates, builds and add to draw list an object for solids visualization
unsigned int Simulation::drawSolids()  
{
    this->objectSolids=new Object(GL_TRIANGLES);
    this->buildSolids();
    unsigned int storedObjectSolids=scene.storeObject(objectSolids);
    unsigned int solidsID=scene.addObjectToDraw(storedObjectSolids);  
    scene.setDrawnObjectShader(solidsID, defaultShader);
    float shade=0.2;
    float grey[]={shade*0.8, shade*0.85, shade*1.0, 0.6};
    scene.setDrawnObjectColor(solidsID, grey);
    scene.setDrawnObjectTransparency(solidsID, true);
    scene.setDrawnObjectShader(solidsID, defaultShader);
    return solidsID;
}

//______________________________________________________________________________
// Main actions


// Simulation update
void Simulation::update()
{
    // Colors semi-lagrangian advection (and potential visualization update)
    this->advectColors();
    
    // Velocities semi-lagrangian advection (and potential visualization update)
    this->advectVelocities();
    
    // Extern forces are applied
    this->applyForces();

    // Velocities are uptated to satisfy incompressiblity and limit conditons
    this->project();
}


// Particle visualization for the simulation
void Simulation::render()
{
    // Moves particles according to current velocity field 
    // Updates potential visualization
    if ((objectParticles!=NULL) || (objectParticleVelocities!=NULL))
        this->advectParticles();
}
