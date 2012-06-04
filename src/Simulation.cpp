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
// 3D working
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
	            if (surface) samples[index*4+2]=this->offset[2]+iZ*this->h;
	            else         samples[index*4+2]=this->offset[2]+this->offsetInCell+iZ*this->h;
	            samples[index*4+3]=1.0;
	            
	            pressures[index]=1.0;
	            
	            types[index]=0;
	            if (((iY>=(nbSamplesY-1)) || (iY<(1*nbSamplesY/3))) || (iX<(nbSamplesY/2)))
	                types[index]=1;
	            if ( ((iX>=(nbSamplesX/2+1)) && (iX<=(3*nbSamplesX/4))) 
	              && ((iY>=(nbSamplesY/2)) && (iY<=(nbSamplesY/2+2)))) types[index]=2;
	              
	            forces[index*4+0]=0.0;
	            //forces[index*4+1]=0.0;
	            forces[index*4+1]=-2.0;
	            forces[index*4+2]=0.0;
	            forces[index*4+3]=0.0;
	            //if (iX==(nbSamplesX/4))   forces[index*4+1]=10.0;
	            //if (iX==(3*nbSamplesX/4)) forces[index*4+1]=-10.0;	 
	            if ((iX>=(nbSamplesX/2)) && (iX<=(nbSamplesX/2+2))) forces[index*4+1]=10.0;	 
	                             
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
// 3D not implemented
void Simulation::initVelocitiesBorders()
{
    float coefVelocity=1.0/1.0;

	unsigned int iBordersX, iBordersY;
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
	{
	    for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
	    {   
	        iBordersX=iY*(nbSamplesX+1)+iX;
	        
	        // Inits circular velocity field (velocityX=-borderPosition[1])
            //velocitiesX[iBordersX]=(-(offset[1]+offsetInCell+iY*h)/size)*coefVelocity;
            
            // Inits diagonally growing velocity field (velocityX=k*borderPosition[0])
            //velocitiesX[iBordersX]=(((float)iX/(float)nbSamplesX)/size)*coefVelocity;
            
            velocitiesX[iBordersX]=0.0*coefVelocity;
        }
    }
    for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
	{
	    for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
	    {   
	        iBordersY=iY*nbSamplesX+iX;
	        
	        // Inits circular velocity field (velocityY=borderPosition[0])
            //velocitiesY[iBordersY]=((offset[0]+offsetInCell+iX*h)/size)*coefVelocity;
            
            // Inits diagonally growing velocity field (velocityY=k*borderPosition[1])
            //velocitiesY[iBordersY]=(((float)iY/(float)nbSamplesY)/size)*coefVelocity;
            
            velocitiesY[iBordersY]=0.0*coefVelocity;
        }
    }
    enforceVelocitiesBorders();
}


// Init the particles positions, and interpolates colors and velocities from grid
// 3D not implemented
void Simulation::initParticles()
{
    unsigned int nbAxis=pow(2.0, nbParticlesCoef);  
    
    float halfCell=h/2.0;
    float step=h/float(nbAxis);
    float sideOffsetX=h/(4.0*nbAxis);
    
    unsigned int iParticles=0;
    float * particlesTmp=new float[this->nbParticles*4]; //new
        
    // Inits 4 particles in 2D cells (optimal distribution)
    // if (nbParticlesCoef==1) : corners of a losange in the square cell
    // if (nbParticlesCoef==2) : twice as more ...
	for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
	{
	    if (types[iSamples]==0)
	    {
	        for (unsigned int iY=0 ; iY<nbAxis ; ++iY)
	        {
	            float offsetY=-halfCell + step/2.0 + iY*step;
	            for (unsigned int iX=0 ; iX<nbAxis ; ++iX)
	            {
	                float offsetX=-halfCell + sideOffsetX + iX*step;
	                if (iY%2==1) offsetX+=step/2.0;

	                particlesTmp[iParticles*4+0]=samples[iSamples*4+0]+offsetX;
	                particlesTmp[iParticles*4+1]=samples[iSamples*4+1]+offsetY;
	                particlesTmp[iParticles*4+2]=samples[iSamples*4+2];
	                particlesTmp[iParticles*4+3]=samples[iSamples*4+3];
	                iParticles++;
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
    
     
	for (unsigned int iParticles=0 ; iParticles<this->nbParticles ; ++iParticles)
	{
	    // Moves particles randomly and at close range around initial position
	    float radius=getRand(*this->rand, true, 0.0, 2.0*sideOffsetX);
	    float angle=getRand(*this->rand, true, 0.0, 2.0*M_PI);
	    particles[iParticles*4+0]+=radius*cos(angle);
	    particles[iParticles*4+1]+=radius*sin(angle);

        // Interpolates color and velocity bilinearly for each particle
        for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        {
            particleColors[iParticles*4+iCoord]=0.0;
            if (iCoord==3) particleColors[iParticles*4+iCoord]=1.0;
            particleVelocities[iParticles*4+iCoord]=0.0;
        }

        //interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0]));
        particleColors[iParticles*4+0]=1.0;
        particleColors[iParticles*4+3]=1.0;
        if (iParticles<2*nbParticles/3) particleColors[iParticles*4+1]=1.0;
        if (iParticles<1*nbParticles/3) particleColors[iParticles*4+2]=1.0;
            
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), &(particleVelocities[iParticles*4+0]));  
    }
}


//______________________________________________________________________________
// Elements building


// Builds an object to visualize the samples (positions and colors)
// The object will then be accessed and updated each frame (advectColors(...))
// 3D working
void Simulation::buildSamples()
{
    std::cout<<"Building simulation : samples positions."<<std::endl;
    
    objectSamples->setNbVertices(nbSamples);
    
    // No indices are necessary since we draw GL_POINTS
	unsigned int * indices=NULL;
	
    // Sends the data into buffers on the GPU
    objectSamples->sendPrimitives(samples, indices);
    objectSamples->sendColors(colors, true);
}


// Builds an object to visualize the pressures
// The object will then be accessed and updated each frame (project(...))
// 3D working
void Simulation::buildPressures()
{
    std::cout<<"Building simulation : pressures."<<std::endl;
    
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
    std::cout<<"Building simulation : cell forces."<<std::endl;
    
    objectForces->setNbVertices(nbSamples*2);
    objectForces->setNbIndices(objectForces->getNbVertices());
    
    float vertices[objectForces->getNbVertices()*4];
    float colors[objectForces->getNbVertices()*4];
    unsigned int indices[objectForces->getNbIndices()];
    
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
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
    // Sends the data into buffers on the GPU
    objectForces->sendPrimitives(vertices, indices);
    objectForces->sendColors(colors);
}


// Builds an object to visualize the types
// 3D working
void Simulation::buildTypes()
{
    std::cout<<"Building simulation : samples types."<<std::endl;
    
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
// 3D bugged
void Simulation::buildGrid()
{
    std::cout<<"Building simulation : grid borders."<<std::endl;
    unsigned int tmpNbSamplesZ=nbSamplesZ;
    if (surface) tmpNbSamplesZ=0;

    objectGrid->setNbVertices((nbSamplesX+1)*(nbSamplesY+1)*(tmpNbSamplesZ+1));
    objectGrid->setNbIndices((nbBordersX+nbBordersY+nbBordersZ)*2);
    float vertices[objectGrid->getNbVertices()*4];
    unsigned int indices[objectGrid->getNbIndices()];
    
    unsigned int iIndices=0;
    for (unsigned int iZ=0 ; iZ<(tmpNbSamplesZ+1) ; ++iZ)
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
                if ((!surface) && (iZ<nbSamplesZ)) // z oriented segments
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
// 3D not implemented
void Simulation::buildVelocitiesBorders()
{
    std::cout<<"Building simulation : cell borders velocities."<<std::endl;
    
    objectVelocitiesBorders->setNbVertices(nbBordersX+nbBordersY+nbBordersZ);
	unsigned int * indices=NULL;
	
    float vertices[objectVelocitiesBorders->getNbVertices()*4];
    float colors[objectVelocitiesBorders->getNbVertices()*4];

	unsigned int iBorders=0;
	unsigned int iBordersX, iBordersY;
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
	{
	    for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
	    {   
	        iBordersX=iY*(nbSamplesX+1)+iX;
            // Velocity x components are sampled on at the middle 
            // of the yz oriented cell borders (y in 2D)
            vertices[iBorders*4+0]=offset[0]+iX*h;
            vertices[iBorders*4+1]=offset[1]+offsetInCell+iY*h;
            vertices[iBorders*4+2]=offset[2];
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
    
    for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
	{
	    for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
	    {   
            iBordersY=iY*nbSamplesX+iX;
            // Velocity y components are sampled on at the middle 
            // of the xz oriented cell borders (x in 2D)
            vertices[iBorders*4+0]=offset[0]+offsetInCell+iX*h;
            vertices[iBorders*4+1]=offset[1]+iY*h;
            vertices[iBorders*4+2]=offset[2];
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
    // Sends the data into buffers on the GPU
    objectVelocitiesBorders->sendPrimitives(vertices, indices);
    objectVelocitiesBorders->sendColors(colors, true);
}


// Cancels normal velocities on solid boundaries
// 3D not implemented
void Simulation::enforceVelocitiesBorders()
{
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;

            // If solid cell
            if (types[iSamples]==2)   
            { 
                //  normal velocity component must be null
                velocitiesX[indexVelocitiesXLeft]  =0.0;
                velocitiesX[indexVelocitiesXRight] =0.0;
                velocitiesY[indexVelocitiesYBottom]=0.0;
                velocitiesY[indexVelocitiesYTop]   =0.0;
            }
            
            if (solidWalls)
            {
                if (iX==0)              velocitiesX[indexVelocitiesXLeft]  =0.0;
                if (iX==(nbSamplesX-1)) velocitiesX[indexVelocitiesXRight] =0.0;
                if (iY==0)              velocitiesY[indexVelocitiesYBottom]=0.0;
                if (iY==(nbSamplesY-1)) velocitiesY[indexVelocitiesYTop]   =0.0;
            }
        }
    }
}


// Builds an object to visualize the velocity vectors interpolated at cell centers
// The object will then be accessed and updated each frame (advectVelocities(...))
// 3D not implemented
void Simulation::buildVelocitiesCenters()
{
    std::cout<<"Building simulation : cell centers velocities."<<std::endl;
    
    objectVelocitiesCenters->setNbVertices(nbSamples*2);
    objectVelocitiesCenters->setNbIndices(objectVelocitiesCenters->getNbVertices());
    
    float vertices[objectVelocitiesCenters->getNbVertices()*4];
    float colors[objectVelocitiesCenters->getNbVertices()*4];
    unsigned int indices[objectVelocitiesCenters->getNbIndices()];
    
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
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
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;

            float velocity[]={0.0, 0.0, 0.0, 0.0};
            // X component is interpolated from left and right
            velocity[0]=(velocitiesX[indexVelocitiesXLeft]
                        +velocitiesX[indexVelocitiesXRight])/2.0;
            // Y component is interpolated from bottom and top           
            velocity[1]=(velocitiesY[indexVelocitiesYBottom]
                        +velocitiesY[indexVelocitiesYTop])/2.0;
                        
            // The vector is printed inversed and transparent at the end (for effect)
            vertices[index*4+0]=samples[iSamples*4+0]-velocity[0]*vectorScale;
            vertices[index*4+1]=samples[iSamples*4+1]-velocity[1]*vectorScale;
            vertices[index*4+2]=samples[iSamples*4+2]-velocity[2]*vectorScale;
            vertices[index*4+3]=1.0;
            colors[index*4+0]=0.0;
            colors[index*4+1]=0.0;
            colors[index*4+2]=0.0;
            colors[index*4+3]=1.0;
            indices[index]=index;
        }  
    }
    // Sends the data into buffers on the GPU
    objectVelocitiesCenters->sendPrimitives(vertices, indices, true);
    objectVelocitiesCenters->sendColors(colors);
}


// Builds an object to visualize the particles 
// The object will then be accessed and updated each frame (advectParticles(...))
// 3D working theorically
void Simulation::buildParticles()
{
    std::cout<<"Building simulation : particles."<<std::endl;
    
    objectParticles->setNbVertices(nbParticles);
	unsigned int * indices=NULL;
	
    // Sends the data into buffers on the GPU
    objectParticles->sendPrimitives(particles, indices, true);
    objectParticles->sendColors(particleColors, true);
}


// Builds an object to visualize the particles interpolated velocities
// The object will then be accessed and updated each frame (advectParticles(...))
// 3D not implemented
void Simulation::buildParticleVelocities()
{
    std::cout<<"Building simulation : particles velocities."<<std::endl;
    
    objectParticleVelocities->setNbVertices(nbParticles*2);
    objectParticleVelocities->setNbIndices(objectParticleVelocities->getNbVertices());
    
	float colors[objectParticleVelocities->getNbVertices()*2*4];
	unsigned int indices[objectParticleVelocities->getNbIndices()];
	
	for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
	{
	    indices[iParticles*2+0]=iParticles*2+0;
	    indices[iParticles*2+1]=iParticles*2+1;
	    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
	    {
	        // Vector start on particle location
	        particleVelocities[(iParticles*2+0)*4+iCoord]=particles[iParticles*4+iCoord];
	        // Vector starts with the particle interpolated color
	        colors[(iParticles*2+0)*4+iCoord]=particleColors[iParticles*4+iCoord];
	        
	        // Bilinear inteprolation of velocity
	        float particleVelocity[]={0.0, 0.0, 0.0, 0.0};
	        if (surface)
            {
                interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), particleVelocity);
            }else{}
            
            // The vector is printed inversed and transparent at the end (for effect)
	        particleVelocities[(iParticles*2+1)*4+iCoord]=particles[iParticles*4+iCoord]-(particleVelocity[iCoord]*vectorScale);
	        colors[(iParticles*2+1)*4+iCoord]=0.0;
	    }
	}
    // Sends the data into buffers on the GPU
    objectParticleVelocities->sendPrimitives(particleVelocities, indices, true);
    objectParticleVelocities->sendColors(colors, true);
}


//______________________________________________________________________________
// Interpolation


// Interpolates (bi/tri)linearly a data field sampled on cell centers
// 3D not implemented
void Simulation::interpolateFromCenters(const float * const data, const float * const position, float * const result)
{
    float xCorner=(position[0]-offset[0])/h;
    float yCorner=(position[1]-offset[1])/h;
    if ((xCorner<=0.0) || (xCorner>=nbSamplesX) || (yCorner<=0.0) || (yCorner>=nbSamplesY))
        return;
        
    // Left-bottom closest sample
    float x=xCorner-0.5;
    float y=yCorner-0.5;
    
    // Distances to left-bottom closest sample
    float toLeft  =x-floor(x);
    float toBottom=y-floor(y);    
    
    // Corresponding indices on x/y axis
    int iX=int(floor(x));
    int iY=int(floor(y));
    
    // Four neighbours indices
    int x0y0= iY   *nbSamplesX+iX  ;
    int x1y0= iY   *nbSamplesX+iX+1;
    int x0y1=(iY+1)*nbSamplesX+iX  ;
    int x1y1=(iY+1)*nbSamplesX+iX+1;
    if (x<0.0)              { x0y0=x1y0; x0y1=x1y1;}
    if (x>(nbSamplesX-1.0)) { x1y0=x0y0; x1y1=x0y1;}
    if (y<0.0)              { x0y0=x0y1; x1y0=x1y1;}
    if (y>(nbSamplesY-1.0)) { x0y1=x0y0; x1y1=x1y0;}
    // Bilinear interpolation (4 components) from the four neighbours
    biLinearInterpolation(4, &(data[x0y0*4]), &(data[x1y0*4]), &(data[x0y1*4]), &(data[x1y1*4]), toLeft, toBottom, result);
}


// Interpolates (bi/tri)linearly a data field sampled on cell borders per components
// 3D not implemented
void Simulation::interpolateFromBorders(const float * const dataX, const float * const dataY, const float * const dataZ, const float * const position, float * const result)
{
    float xCorner=(position[0]-offset[0])/h;
    float yCorner=(position[1]-offset[1])/h;
    if ((xCorner<=0.0) || (xCorner>=nbSamplesX) || (yCorner<=0.0) || (yCorner>=nbSamplesY))
        return;

    float x=xCorner;
    float y=yCorner-0.5;
    // Distances to left closest sample and corresponding y  
    float toLeft  =x-floor(x);
    float toBottom=y-floor(y);
    
    // Corresponding indices on x/y axis
    int indexLeft  =int(floor(x));
    int indexBottom=int(floor(y)); 

    // Corresponding index in dataX
    int indexDataXLeft=indexBottom*(nbSamplesX+1)+indexLeft;
    
    // Four dataX neighbours indices
    int x0y0=indexDataXLeft;
    int x1y0=indexDataXLeft+1;
    int x0y1=indexDataXLeft +(nbSamplesX+1);
    int x1y1=indexDataXLeft+(nbSamplesX+1)+1;
    if (y<0.0)              { x0y0=x0y1; x1y0=x1y1;}
    if (y>(nbSamplesY-1.0)) { x0y1=x0y0; x1y1=x1y0;}
    // Bilinear interpolation (1 x components) from the four neighbours
    biLinearInterpolation(1, &(dataX[x0y0]), &(dataX[x1y0]), &(dataX[x0y1]), &(dataX[x1y1]), toLeft, toBottom, &(result[0]));
   
    x=xCorner-0.5;
    y=yCorner;
    // Distances to bottom closest sample and corresponding x   
    toLeft  =x-floor(x);
    toBottom=y-floor(y);
    
    // Corresponding indices on x/y axis     
    indexBottom=int(floor(y));
    indexLeft  =int(floor(x));
    
    // Corresponding index in dataY
    int indexDataYBottom=indexBottom*nbSamplesX+indexLeft;
    // Four dataY neighbours indices   
    x0y0=indexDataYBottom;
    x1y0=indexDataYBottom+1;
    x0y1=indexDataYBottom+nbSamplesX;
    x1y1=indexDataYBottom+nbSamplesX+1;
    if (x<0.0)              { x0y0=x1y0; x0y1=x1y1;}
    if (x>(nbSamplesX-1.0)) { x1y0=x0y0; x1y1=x0y1;}
    // Bilinear interpolation (1 y components) from the four neighbours
    biLinearInterpolation(1, &(dataY[x0y0]), &(dataY[x1y0]), &(dataY[x0y1]), &(dataY[x1y1]), toLeft, toBottom, &(result[1]));
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
// 3D not implemented
unsigned int Simulation::type(int iX, int iY, int iZ)
{
    if ( (iX<0) || (iX>=(int)nbSamplesX) || (iY<0) || (iY>=(int)nbSamplesY) )
    {
        if (solidWalls) return 2;
        else return 1;
    }
    return types[iY*nbSamplesX+iX];
}


// Sets integer indices of the cell position pos overlaps
// 3D not implemented
void Simulation::getCell(const float * const pos, int * const iX, int * const iY, int * const iZ)
{
    (*iX)=int(floor((pos[0]-offset[0])/h));
    (*iY)=int(floor((pos[1]-offset[1])/h));
}


//______________________________________________________________________________
// System resolution


// Evaluates Modified Incomplete Choleski (level 0) preconditioner
void Simulation::MICPreconditioner(const double * const Adiag, const double * const Aright, const double * const Atop, Eigen::VectorXd  & precon)
{
    double tau=0.97;
    double safetyConstant=0.25;
    for (int iY=0 ; iY<(int)nbSamplesY ; ++iY)
    {
        for (int iX=0 ; iX<(int)nbSamplesX ; ++iX)
        {
            int iS=iY*nbSamplesX+iX;

            int iSLeft=iY*nbSamplesX+(iX-1);
            int iSBottom=(iY-1)*nbSamplesX+iX;
            double e=Adiag[iS];
            if (iX>0)
            {
                e-=pow(Aright[iSLeft]*(precon)[iSLeft], 2)
                  +tau*(Aright[iSLeft]  *Atop[iSLeft]  *pow((precon)[iSLeft]  , 2));
            }
            if (iY>0)
            {                  
                e-=pow(Atop[iSBottom]*(precon)[iSBottom], 2)
                  +tau*(Atop  [iSBottom]*Aright[iSBottom]*pow((precon)[iSBottom], 2));
            }  
            if (e<(safetyConstant*Adiag[iS])) e=Adiag[iS];
            if (e!=0.0) precon[iS]=1.0/sqrt(e);
        }
    }
}
    

// Multiplies preconditionneur to r and stores it in z
void Simulation::applyPreconditioner(const double * const Aright, const double * const Atop, const Eigen::VectorXd & precon, const Eigen::VectorXd & r, Eigen::VectorXd & z)
{
    double t=0.0;
    Eigen::VectorXd q=Eigen::VectorXd::Zero(nbSamples);
    for (int iY=0 ; iY<(int)nbSamplesY ; ++iY)
    {
        for (int iX=0 ; iX<(int)nbSamplesX ; ++iX)
        {
            int iS=iY*nbSamplesX+iX;
            int iSLeft=iY*nbSamplesX+(iX-1);
            int iSBottom=(iY-1)*nbSamplesX+iX;
            t=r[iS];
            if (iX>0) t-=Aright[iSLeft]*precon[iSLeft]  *q[iSLeft];
            if (iY>0) t-=Atop[iSBottom]*precon[iSBottom]*q[iSBottom];
            q[iS]=t*precon[iS];
        }
    }
    for (int iY=((int)nbSamplesY-1) ; iY>=0; --iY)
    {
        for (int iX=((int)nbSamplesX-1) ; iX>=0 ; --iX)
        {
            int iS=iY*nbSamplesX+iX;
            int iSRight=iY*nbSamplesX+(iX+1);
            int iSTop=(iY+1)*nbSamplesX+iX;
            t=q[iS];
            if (iX<((int)nbSamplesX-1)) t-=Aright[iS]*precon[iS]*z[iSRight];
            if (iY<((int)nbSamplesY-1)) t-=Atop[iS]  *precon[iS]*z[iSTop];
            z[iS]=t*precon[iS];
        }
    }
}


// Multiplies sparse matrix to v and stores to result
void Simulation::multiplySparseMatrix(const double * const Adiag, const double * const Aright, const double * const Atop, const Eigen::VectorXd & v, Eigen::VectorXd & result)
{
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
            result[iSamples]=Adiag[iSamples]*v[iSamples];
            
            if (iX<(nbSamplesX-1)) // Right
                result[iSamples]+=Aright[iSamples]           *v[iSamples+1         ];
            if (iY<(nbSamplesY-1)) // Top
                result[iSamples]+=Atop  [iSamples]           *v[iSamples+nbSamplesX];
            if (iX>0)                 // Left
                result[iSamples]+=Aright[iSamples-1]         *v[iSamples-1         ];
            if (iY>0)                 // Bottom
                result[iSamples]+=Atop  [iSamples-nbSamplesX]*v[iSamples-nbSamplesX];            
        }
    }
}


// Fills pressures from sparse matrix and rhs=b
void Simulation::conjugateGradient(const double * const Adiag, const double * const Aright, const double * const Atop, const double * const b)
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
    MICPreconditioner(Adiag, Aright, Atop, precon);

    // Apply preconditionner from r to z
    Eigen::VectorXd z=Eigen::VectorXd::Zero(nbSamples);
    applyPreconditioner(Aright, Atop, precon, r, z);
    // Copy z to s
    Eigen::VectorXd s(z);
    double sigma=z.dot(r);
    unsigned int iLast=maxIterations;
    for (unsigned int i=0 ; i<maxIterations ; ++i)
    {
        // Multiply matrix A to s to get z;
        multiplySparseMatrix(Adiag, Aright, Atop, s, z);
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
            applyPreconditioner(Aright, Atop, precon, r, z);
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
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;

            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;

            // Sets Right Hand Side        
            divergences[iSamples]=scale*( velocitiesX[indexVelocitiesXRight]
                                         -velocitiesX[indexVelocitiesXLeft]
                                         +velocitiesY[indexVelocitiesYTop]
                                         -velocitiesY[indexVelocitiesYBottom] );
        }
    }
}


//______________________________________________________________________________
// Main steps


// Modifies the color samples so that color derivative is null (particles moving through the grid keep their color)
// 3D not implemented
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectColors(float dt)
{
    float newColors[nbSamples*4];
    
    float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
    float virtualParticlePos[]={0.0, 0.0, 0.0, 1.0};
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
            
            // Interpolates the velocity at cell center (sampleVelocity)
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                              +velocitiesX[indexVelocitiesXRight])/2.0;
            sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                              +velocitiesY[indexVelocitiesYTop])/2.0;
            
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
            if (type(iXvirtual, iYvirtual, 0)==0)             
                interpolateFromCenters(colors, virtualParticlePos, &(newColors[iSamples*4+0]));
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
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
        {
            for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
            {
                unsigned int iSamples=iY*nbSamplesX+iX;
                
                float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
                
                // Interpolates the velocity at cell center (sampleVelocity)
                unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
                unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
                unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
                unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
                sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                                  +velocitiesX[indexVelocitiesXRight])/2.0;
                sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                                  +velocitiesY[indexVelocitiesYTop])/2.0;
                                  
                for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                {
                    vertices[(iSamples*2+0)*4+iCoord]=samples[iSamples*4+iCoord];
                    vertices[(iSamples*2+1)*4+iCoord]=samples[iSamples*4+iCoord]-sampleVelocity[iCoord]*vectorScale;
                }
            }
        }
        objectVelocitiesCenters->updateVertices(vertices, true);
    }
    
    if (objectVelocitiesBorders!=NULL)
    {
        float colors[objectVelocitiesBorders->getNbVertices()*4];
        unsigned int iBorders=0;
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
	    {
	        for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
	        {   
                float normVelocityX=velocitiesX[iY*(nbSamplesX+1)+iX];
                if (normVelocityX<0.0) normVelocityX=-normVelocityX;
                colors[iBorders*4+0]=normVelocityX*10.0;//new
                colors[iBorders*4+1]=0.0; 
                colors[iBorders*4+2]=0.0; 
                colors[iBorders*4+3]=1.0;
                iBorders++;
            }
        }
        
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
	    {
	        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
	        {   
                float normVelocityY=velocitiesY[iY*nbSamplesX+iX];
                if (normVelocityY<0.0) normVelocityY=-normVelocityY;
                colors[iBorders*4+0]=0.0;
                colors[iBorders*4+1]=normVelocityY*10.0;//new
                colors[iBorders*4+2]=0.0; 
                colors[iBorders*4+3]=1.0;
                iBorders++;
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
    
    // Uses the new centers velocity to update the separated velocity components
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            
            float coefLeft=0.5;
            float coefRight=0.5;
            float coefBottom=0.5;
            float coefTop=0.5;
            if (iX==0)              coefLeft  =1.0;
            if (iX==(nbSamplesX-1)) coefRight =1.0;
            if (iY==0)              coefBottom=1.0;
            if (iY==(nbSamplesY-1)) coefTop   =1.0;            
            
            velocitiesX[indexVelocitiesXLeft]  +=coefLeft  *centeredVel[iSamples*4+0];
            velocitiesX[indexVelocitiesXRight] +=coefRight *centeredVel[iSamples*4+0];
            velocitiesY[indexVelocitiesYBottom]+=coefBottom*centeredVel[iSamples*4+1];
            velocitiesY[indexVelocitiesYTop]   +=coefTop   *centeredVel[iSamples*4+1];
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
        for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
	    {
	        for (unsigned int iX=0 ; iX<(nbSamplesX+1) ; ++iX)
	        {   
                float normVelocityX=velocitiesX[iY*(nbSamplesX+1)+iX];
               // if (normVelocityX<0.0) normVelocityX=-normVelocityX;
                colors[iBorders*4+0]=normVelocityX;
                colors[iBorders*4+1]=0.0; 
                colors[iBorders*4+2]=-normVelocityX; 
                colors[iBorders*4+3]=1.0;
                iBorders++;
            }
        }
        
        for (unsigned int iY=0 ; iY<(nbSamplesY+1) ; ++iY)
	    {
	        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
	        {   
                float normVelocityY=velocitiesY[iY*nbSamplesX+iX];
                //if (normVelocityY<0.0) normVelocityY=-normVelocityY;
                colors[iBorders*4+0]=0.0;
                colors[iBorders*4+1]=normVelocityY;
                colors[iBorders*4+2]=-normVelocityY; 
                colors[iBorders*4+3]=1.0;
                iBorders++;
            }
        }
        objectVelocitiesBorders->updateColors(colors, true);
    }
}


// Modifies the velocity samples so that velocity derivative is null (particles moving through the grid keep their velocity)
// 3D not implemented
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectVelocities(float dt)
{
    float newVelocities[nbSamples*4];
    
    float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
    float virtualParticlePos[]={0.0, 0.0, 0.0, 1.0};
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;
            
            // Interpolates the velocity at cell center (sampleVelocity)
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                              +velocitiesX[indexVelocitiesXRight])/2.0;
            sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                              +velocitiesY[indexVelocitiesYTop])/2.0;
            
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
            newVelocities[iSamples*4+3]=sampleVelocity[3];
            
            int iXvirtual, iYvirtual, iZvirtual;
            getCell(virtualParticlePos, &iXvirtual, &iYvirtual, &iZvirtual);
            if (type(iXvirtual, iYvirtual, 0)==0)
                interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, virtualParticlePos, &(newVelocities[iSamples*4+0]));
        }  
    }
    updateVelocitiesFromCenters(newVelocities);
}


// Change velocities according to forces
void Simulation::applyForces(float dt)
{
    float newVelocities[nbSamples*4];
    
    // Uses the new centers velocity to update the separated velocity components
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;   
                        
            float sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
            
            // Interpolates the velocity at cell center (sampleVelocity)
            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]
                              +velocitiesX[indexVelocitiesXRight])/2.0;
            sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]
                              +velocitiesY[indexVelocitiesYTop])/2.0;
                                  
            for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                newVelocities[iSamples*4+iCoord]=sampleVelocity[iCoord];
                
            if (type(iX, iY, 0)==0)
            {        
                for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
                    newVelocities[iSamples*4+iCoord]+=dt*forces[iSamples*4+iCoord];
            }
        }
    }        
    updateVelocitiesFromCenters(newVelocities);
}


// Projection : substracting the pressure gradient
// Pressures evaluated so that incompressibility and boundary conditions are enforced 
void Simulation::project(float dt)
{
    // To evaluate Right Hand Side of equation
    double scaleRhs=1.0/h;
    double rhs[nbSamples];
    
    // To evaluate A matrix
    double scaleA=dt/(density*h*h);
    double Adiag[nbSamples];
    double Aright[nbSamples];
    double Atop[nbSamples];
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
    {
        rhs[iSamples]=0.0;
        Adiag[iSamples]=0.0;
        Aright[iSamples]=0.0;
        Atop[iSamples]=0.0;
    }
    
    double velocitySolids=0.0;
    setDivergences(rhs);

    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;

            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            
            unsigned int typeSample=type(int(iX), int(iY), 0);
            unsigned int typeNeighbourLeft  =type(int(iX)-1, int(iY),   0);
            unsigned int typeNeighbourRight =type(int(iX)+1, int(iY),   0);
            unsigned int typeNeighbourBottom=type(int(iX),   int(iY)-1, 0);
            unsigned int typeNeighbourTop   =type(int(iX),   int(iY)+1, 0);   
            
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

                // Sets matrix A compact form
                if (typeNeighbourLeft!=2)   Adiag[iSamples]+=scaleA;
                if (typeNeighbourRight!=2)  Adiag[iSamples]+=scaleA;
                if (typeNeighbourBottom!=2) Adiag[iSamples]+=scaleA;
                if (typeNeighbourTop!=2)    Adiag[iSamples]+=scaleA;
                    
                if (typeNeighbourRight==0) Aright[iSamples]=-scaleA;
                if (typeNeighbourTop==0)   Atop  [iSamples]=-scaleA;
            }
            else
            {
                rhs[iSamples]=0.0;
            }
        }
    }
    // Resolves system with preconditionned conjugate gradient method
    conjugateGradient(Adiag, Aright, Atop, rhs);
    
    // Substracts pressure gradient
    float scale=dt/(density*h);
    for (unsigned int iY=0 ; iY<nbSamplesY ; ++iY)
    {
        for (unsigned int iX=0 ; iX<nbSamplesX ; ++iX)
        {
            unsigned int iSamples=iY*nbSamplesX+iX;

            unsigned int indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            unsigned int indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            unsigned int indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            unsigned int indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;
            
            if (types[iSamples]==0)
            {
                velocitiesX[indexVelocitiesXLeft]  -=scale*pressures[iSamples];
                velocitiesX[indexVelocitiesXRight] +=scale*pressures[iSamples];
                velocitiesY[indexVelocitiesYBottom]-=scale*pressures[iSamples];
                velocitiesY[indexVelocitiesYTop]   +=scale*pressures[iSamples];
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
// 3D working theorically
void Simulation::advectParticles(float dt)
{
    for (unsigned int iSamples=0 ; iSamples<nbSamples ; ++iSamples)
        if (types[iSamples]!=2) types[iSamples]=1;
        
    //float velocityColors[nbParticles*4*2];
    for (unsigned int iParticles=0 ; iParticles<nbParticles ; ++iParticles)
	{
	    int iX, iY, iZ;
        getCell(&(particles[iParticles*4+0]), &iX, &iY, &iZ);
            if (type(iX, iY, 0)!=2)
                types[iY*nbSamplesX+iX]=0;
        
        float knownVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), knownVelocity);

        // What is the new position after dt ?
        integrate(1.0,                            // direction
                &(particles[iParticles*4+0]),     // startValue
                  knownVelocity,                  // startVariation 
                  dt,                             // dt
                &(particles[iParticles*4+0])); // resultValue*/ 

        
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
    }

    // Updates the corresponding data if stored on GPU
    if (objectParticles!=NULL)
    {
        objectParticles->updateVertices(particles, true);
        //objectParticles->updateColors(particleColors, true);
    }
    if (objectParticleVelocities!=NULL)
    {
         objectParticleVelocities->updateVertices(particleVelocities, true);
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
Simulation::Simulation(Scene & scene,
                       const Shaders & shaders, 
                       const unsigned int defaultShader,
                       const unsigned int spriteShader,
                       const bool surface,
                       const float size,
                       const bool solidWalls,
                       const float density,
                       const float viscosity,
                       const unsigned int nbSamplesX,
                       const unsigned int nbSamplesY,
                       unsigned int nbSamplesZ,
                       const unsigned int nbParticlesCoef)
                       : scene(scene), 
                         shaders(shaders),
                         defaultShader(defaultShader),
                         spriteShader(spriteShader),
                         surface(surface),
                         size(size),
                         solidWalls(solidWalls),
                         density(density),
                         viscosity(viscosity),
                         nbSamplesX(nbSamplesX),
                         nbSamplesY(nbSamplesY),
                         nbSamplesZ(nbSamplesZ),
                         nbParticlesCoef(nbParticlesCoef)
{
        // About the data stored on centers of the MAC grid
    // Samples number on Z axis
    if (this->surface) this->nbSamplesZ=1;
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
           
           
        // About the data stored on borders of the MAC grid
    // Number of left/right cell borders
    this->nbBordersX=(this->nbSamplesX+1)*this->nbSamplesY*this->nbSamplesZ;
    // Number of bottom/top cell borders
    this->nbBordersY=this->nbSamplesX*(this->nbSamplesY+1)*this->nbSamplesZ;
    // Number of back/front cell borders
    this->nbBordersZ=this->nbSamplesX*this->nbSamplesY*(this->nbSamplesZ+1);
    if (this->surface) this->nbBordersZ=0;
    // X velocity components, sampled on left/right borders
    this->velocitiesX=new float[this->nbBordersX];
    // Y velocity components, sampled on bottom/top borders
    this->velocitiesY=new float[this->nbBordersY];
    // Z velocity components, sampled on back/front borders
    this->velocitiesZ=new float[this->nbBordersZ];
    
    
        // About the particles used for visualization (advected on the whole grid)
    // Number of rendered particles (=nbSamples*[4,8]^nbParticlesCoef)
    this->nbParticles=this->nbSamples*(unsigned int)pow(8.0, this->nbParticlesCoef);
    if (this->surface) this->nbParticles=this->nbSamples*(unsigned int)pow(4.0, this->nbParticlesCoef);
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
    if (this->surface) this->offset[2]=this->size/2.0;
    // Offset from left/bottom/back of grid cell
    this->offsetInCell=this->h/2.0;
    // Scale coefficient to display vector lengths
    this->vectorScale=1.0/2.0;
        
        // Simulation parameters
    // Vector of acceleration due to gravity
    this->g[0]= 0.0; this->g[1]=-9.81; this->g[2]=0.0;


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
    if (this->surface) std::cout<<"Simulation : surface of ["<<this->nbSamplesX<<"x"<<this->nbSamplesY<<"]"<<std::endl;
    else         std::cout<<"Simulation : volume of ["<<this->nbSamplesX<<"*"<<this->nbSamplesY<<"*"<<this->nbSamplesZ<<"]."<<std::endl;
    std::cout<<"Simulation : "<<this->nbSamples<<" samples."<<std::endl;
    std::cout<<"Simulation : "<<this->nbParticles<<" particles."<<std::endl;


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
    unsigned int starTextureID=loadTexture("../textures/star.ppm");
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
    return particleVelocitiesID;
}


//______________________________________________________________________________
// Main actions


// Simulation update
void Simulation::update()
{
    float dt=1.0/30.0;
    
    // Colors semi-lagrangian advection (and potential visualization update)
    this->advectColors(dt);
    
    // Velocities semi-lagrangian advection (and potential visualization update)
    this->advectVelocities(dt);
    
    // Extern forces are applied
    this->applyForces(dt);

    // Velocities are uptated to satisfy incompressiblity and limit conditons
    this->project(dt);
}


// Particle visualization for the simulation
void Simulation::render()
{
    float dt=1.0/30.0;
    
    // Moves particles according to current velocity field 
    // Updates potential visualization
    if ((objectParticles!=NULL) || (objectParticleVelocities!=NULL))
        this->advectParticles(dt);
}
