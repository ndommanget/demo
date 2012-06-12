// Application.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Application.hpp"
#include "Tools.hpp"
#include "Builders.hpp"

// OpenGL
#include "glew/glew.h"
#include <SFML/OpenGL.hpp>



//______________________________________________________________________________
// Initialization


// Sets the application parameters and does all the initialisation
void Application::init()
{
    this->simulation=NULL;

    // False as long as we don't want to quit the application
    this->done=false;

    // Frame Counter
    this->cntFrame=0;

    // Simulation update dt
    this->dtSimulation=1.0/100.0;

    // Spherical camera
    this->sphericalCamera=false;

    this->speed=1.0;
    for (unsigned int i=0 ; i<6 ; ++i)
    {
        this->moveTimes[i]=0.0;
        this->movePressed[i]=false;
    }
    
    // Window size initialization (for windowed mode)
    this->windowedWidth=1280;
    this->windowedHeight=800;
    this->minWindowSize=300;
    this->fullScreen=false; 
	
    // Mouse position and scroll data initilaization
    // Positions : floats, origin in center, sides at +/-1 or more
    this->xMousePosition=0.0; 
    this->yMousePosition=0.0;
    this->xMouseOffset=0.0;
    this->yMouseOffset=0.0;
    this->xMouseClick=0.0;
    this->yMouseClick=0.0;
    this->scroll=0;
	
    // Initialisation of SFML and creation of OpenGL context
    initSFMLOpenGL();
    
    // Customize a few OpenGL and SFML states (after context creation)
    customizeStates();

    // Creations and initialisations

    this->shaders.init();
    this->camera.init();
    this->resize();
    sf::Mouse::setPosition(this->screenCenter, *this->window);
    this->moveCamera();

    this->scene.init();
}


// Inits SFML and OpenGL context
void Application::initSFMLOpenGL()
{
    // Context parameters :
    // depthBits / stencilBits / antialiasingLevel / major version / minor version  
    sf::ContextSettings wantedSettings(16, 0, 4, 3, 2);

    // Window and context creation
    if (this->fullScreen)
        this->window=new sf::Window(sf::VideoMode::getDesktopMode(),
                                    this->title,
                                    sf::Style::Fullscreen,
                                    wantedSettings);
    else
        this->window=new sf::Window(sf::VideoMode(windowedWidth, windowedHeight, 32),                       this->title, 
                                    sf::Style::Resize | sf::Style::Close, 
                                    wantedSettings);


    // Actual size of the inside of the window
    this->width=this->window->getSize().x;
    this->height=this->window->getSize().y;

    // Context verification
    /*const sf::ContextSettings& obtainedSettings=this->window->getSettings();
    std::cout<<"DepthBits="<<obtainedSettings.depthBits<<std::endl;
    std::cout<<"StencilBits="<<obtainedSettings.stencilBits<<std::endl;
    std::cout<<"AntialiasingLevel="<<obtainedSettings.antialiasingLevel<<std::endl;
    std::cout<<"Version="<<obtainedSettings.majorVersion<<"."<<obtainedSettings.minorVersion<<std::endl;*/

    // Limits the framerate to 40 frames per second
    //this->window->setFramerateLimit(1);
    // Limits the framerate to the display framerate
    this->window->setVerticalSyncEnabled(true);
}


// Customizes a few OpenGL states to fit the application's needs
void Application::customizeStates()
{
    // Glew initialisation : to register all available extentions
    glewExperimental=true; // Glew 1.70 is not yet right for 3.2 core
    GLenum glewError=glewInit();
    if (glewError!=GLEW_OK)
		std::cout<<"GLEW Error : "<<glewGetErrorString(glewError)<<std::endl;
    glGetError(); // A GL_INVALID_ENUM error is generated with or without the 3.2 core fix

    // Puts the window top-left corner at the following coordinates 
    window->setPosition(sf::Vector2i(0, 0));

    // Initialization of the mouse position in the middle of the window	
    this->screenCenter.x=this->width/2;
    this->screenCenter.y=this->height/2;
    if (!sphericalCamera)
    {
        window->setMouseCursorVisible(false);
        sf::Mouse::setPosition(this->screenCenter, *this->window);
    }

    // Background color
    glClearColor(0.0, 0.0, 0.0, 1.0);

    // Sets the with of the lines
    glEnable(GL_LINE_SMOOTH);
    //glLineWidth(2.0); // Not working on Mac OS X in 3.2 core

    glPointSize(10.0);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // Blending function (used for sprites)
    //glBlendFunc(GL_ONE, GL_ONE);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Disables culling
    //glDisable(GL_CULL_FACE);
    glEnable(GL_CULL_FACE);
}


//______________________________________________________________________________
// User interface


// Listens to events during the whole time of the application
// and distributes corresponding tasks
void Application::handleEvent(const sf::Event & event)
{
    bool tooSmall=false;
    switch(event.type) 
    {             
        // Key presses
        case sf::Event::KeyPressed:
            handleKeyEvent(event, true);
            break;
                
        // Key releases    
        case sf::Event::KeyReleased:
            handleKeyEvent(event, false);
            break;                
                
        // Mouse button pressed              
        case  sf::Event::MouseButtonPressed: 
            if ((event.mouseButton.button==sf::Mouse::Left) && sphericalCamera)
            {
                this->xMouseClick=-2.0*(screenCenter.x-event.mouseButton.x)/(float)this->width;
                this->yMouseClick=2.0*(screenCenter.y-event.mouseButton.y)/(float)this->height;
            }
            break;

        // Mouse button released              
        case  sf::Event::MouseButtonReleased: 
            if ((event.mouseButton.button==sf::Mouse::Left) && (sphericalCamera))
            {
                this->moveCamera();
                this->xMouseOffset=this->xMousePosition;
                this->yMouseOffset=this->yMousePosition;
            }
            break;
        
        // Scrollwheel
        case  sf::Event::MouseWheelMoved:
                this->scroll+=event.mouseWheel.delta; // Doesn't work right on MacOS, because delta is an int and the delta of one tick is inferior to 1.
                this->moveCamera();
            break;    
                
        // Window resize                  
        case  sf::Event::Resized:
            this->windowedWidth=this->window->getSize().x;
            this->windowedHeight=this->window->getSize().y;
            if (this->windowedWidth<this->minWindowSize)
            {
                this->windowedWidth=this->minWindowSize;
                tooSmall=true;
            }
            if (this->windowedHeight<this->minWindowSize)
            {
                this->windowedHeight=this->minWindowSize;
                tooSmall=true;
            }
            if (tooSmall)
                this->window->setSize(sf::Vector2u(this->windowedWidth, this->windowedHeight));
            resize();
            break;                      
        
        // Quit event
        case sf::Event::Closed:
            //std::cout<<"Closing"<<std::endl;
            this->done=true;
            break;
                
        default:
            break;
    }
}


// Distributes task for the keyboard events 
// down is true when the key is pressed (false when released)
void Application::handleKeyEvent(const sf::Event & event, bool down)
{
    if (down)
    {
        switch(event.key.code)
        {
            case sf::Keyboard::Escape:
                this->done=true;
            break;
            
            case sf::Keyboard::C :
                this->switchCameraMode();
            break;

            case sf::Keyboard::P :
                this->camera.switchCameraProjection();
            break;            
            
            case sf::Keyboard::W :
                this->switchWireframe();
            break;
            
            case sf::Keyboard::F :
                this->printFPS();
            break;
            
            case sf::Keyboard::F5 :
                this->switchFullScreen();
            break;
            
            case sf::Keyboard::Space :
                //if (this->simulation!=NULL) this->simulation->update();
            break;              
            
            default:
            break;
        }
    }

    bool move=false;
    switch(event.key.code)
    {
        // Camera controls
        // Keys are associated with x, y, z direction for move.     
        
        // On x axis :
        case sf::Keyboard::Right :
        case sf::Keyboard::D :
            //std::cout<<"Right."<<std::endl;
            moveTimes[0]=frameClock.getElapsedTime().asSeconds();
            movePressed[0]=down;
            move=true;
        break;
        case sf::Keyboard::Left :
        case sf::Keyboard::Q :
            //std::cout<<"Left."<<std::endl;
            moveTimes[1]=frameClock.getElapsedTime().asSeconds();
            movePressed[1]=down;
            move=true;
        break;
        
        // On y axis :
        case sf::Keyboard::PageUp :
            //std::cout<<"Up."<<std::endl;
            moveTimes[2]=frameClock.getElapsedTime().asSeconds();
            movePressed[2]=down;
            move=true;
        break;
        case sf::Keyboard::PageDown :
            //std::cout<<"Down."<<std::endl;
            moveTimes[3]=frameClock.getElapsedTime().asSeconds();
            movePressed[3]=down;
            move=true;
        break; 

        // On z axis :
        case sf::Keyboard::Down :
        case sf::Keyboard::S :
            //std::cout<<"Backward."<<std::endl;
            moveTimes[4]=frameClock.getElapsedTime().asSeconds();
            movePressed[4]=down;
            move=true;
        break; 
        case sf::Keyboard::Up :
        case sf::Keyboard::Z :
            //std::cout<<"Forward."<<std::endl;
            moveTimes[5]=frameClock.getElapsedTime().asSeconds();
            movePressed[5]=down;
            move=true;
        break;

        default:
        break;
    }
    if (move) this->moveCamera();
}


//______________________________________________________________________________
// Camera animation


// Moves camera
void Application::moveCamera()
{
    const sf::Vector2i& mousePosition=sf::Mouse::getPosition(*this->window);
    float xRelMousePosition=-2.0*(screenCenter.x-mousePosition.x)/(float)this->width;
    float yRelMousePosition=2.0*(screenCenter.y-mousePosition.y)/(float)this->height; 

    this->xMousePosition=this->xMouseOffset+xRelMousePosition-this->xMouseClick;
    this->yMousePosition=this->yMouseOffset+yRelMousePosition-this->yMouseClick;

    if ((xRelMousePosition>0.9) || (yRelMousePosition>0.9)
    || (xRelMousePosition<-0.9) || (yRelMousePosition<-0.9))  
    {
        this->xMouseOffset+=xRelMousePosition;
        this->yMouseOffset+=yRelMousePosition;
        sf::Mouse::setPosition(this->screenCenter, *this->window);
    } 

    float elapsedTimes[6];
    for (unsigned int i=0 ; i<6 ; ++i)
    {
        elapsedTimes[i]=0.0;
        if (movePressed[i])
            elapsedTimes[i]=this->frameDuration.asSeconds()-this->moveTimes[i];
    }
    float moves[3];
    moves[0]=this->speed*(elapsedTimes[0]-elapsedTimes[1]);
    moves[1]=this->speed*(elapsedTimes[2]-elapsedTimes[3]);
    moves[2]=this->speed*(elapsedTimes[4]-elapsedTimes[5]); 

    float angleForWindowWidth=M_PI;
    float angleForWindowHeight=M_PI/2.0;
    float angles[2];
    angles[0]=xMousePosition*angleForWindowWidth;
    angles[1]=yMousePosition*angleForWindowHeight;

    float radius=1.0-scroll*0.05;
    this->camera.move(this->sphericalCamera, moves, angles, radius);
}


//______________________________________________________________________________
// Main steps


// Animates camera and scene each frame
void Application::animate()
{
    if ((!sphericalCamera) 
    ||  (sphericalCamera && sf::Mouse::isButtonPressed(sf::Mouse::Left)))
        this->moveCamera();

    for (unsigned int i=0 ; i<6 ; ++i) this->moveTimes[i]=0.0;

    //this->moveScene();
    this->rotateCamera();
}


// Stores an image each frame
void Application::storeFrame() const
{
    std::string fileName="../images/render";
    char iChar[]={0, 0, 0, 0, 0};

    sprintf(iChar, "%d", this->cntFrame);
    if (this->cntFrame<10)    fileName+="0";
    if (this->cntFrame<100)   fileName+="0";
	if (this->cntFrame<1000)  fileName+="0";
    if (this->cntFrame<10000) fileName+="0";
	fileName+=iChar;
	fileName+=".ppm";

	saveFrameBufferPPM(fileName, this->width, this->height);
    //std::cout<<fileName<<" stored"<<std::endl;
}


// Clears the current frame buffer (the image on screen) 
// draws the scene in the other available frame buffer (double buffering)
// and prints the new filled frame buffer on screen
void Application::renderFrame() const
{
    // Clears the window with current clearing color, clears also the depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
    // Draws this->scene
    this->scene.draw();

    // Draws this->simulation
    if (this->simulation!=NULL) this->simulation->render();

    //this->storeFrame();
    
    // Performs the buffer swap between the current shown buffer, 
    // and the one we just worked on (contains swapBuffers ?)
    this->window->display();

    // Reports any possible glError
    printGLErrors();
}


//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Default constructor
Application::Application() : shaders(), camera(this->shaders), scene(this->shaders, this->camera)
{
    this->init();
}


// Cleans before the application can be closed
Application::~Application()
{
    printGLErrors(); std::cout<<std::endl;

    delete window;

    if (simulation!=NULL)  delete simulation;    
}


//______________________________________________________________________________
// Application modes


// Turns ON/OFF wireframe mode
void Application::switchWireframe()
{
    std::cout<<"Wireframe switch."<<std::endl;
    int wireframe[2];
    glGetIntegerv(GL_POLYGON_MODE, wireframe);
    if (wireframe[0]==GL_FILL)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


// Adapts the projection to the new size of the rendering region
void Application::resize()
{
    this->width=this->window->getSize().x;
    this->height=this->window->getSize().y;
    //std::cout<<"Window resize  : ["<<this->width<<","<<this->height<<"]"<<std::endl;

    this->screenCenter.x=this->width/2;
    this->screenCenter.y=this->height/2;
    if (!sphericalCamera)
        sf::Mouse::setPosition(this->screenCenter, *this->window);
    
    // Viewport transformation update to fit initial window size
    glViewport(0, 0, this->width, this->height);

    // Projection transformation update
    // Keeping the angle constant
    float fovy=M_PI/3.0;
    float aspectRatio=this->width/float(this->height);
    this->camera.setPerspectiveFromAngle(fovy, aspectRatio);
}


// Turns ON/OFF fullScreen mode
void Application::switchFullScreen()
{
    std::cout<<"Screen mode switch."<<std::endl;
    if (this->fullScreen)
    {
        window->create(sf::VideoMode(windowedWidth, windowedHeight, 32), this->title, sf::Style::Resize | sf::Style::Close, window->getSettings());
        this->fullScreen=false;
    }
    else
    {    window->create(sf::VideoMode::getDesktopMode(), this->title, sf::Style::Fullscreen, window->getSettings());
        this->fullScreen=true;
    }
    resize();
}


// Swiches between spherical camera mode and FPS camera mode
void Application::switchCameraMode()
{
    std::cout<<"Camera mode switch."<<std::endl;
    if (this->sphericalCamera)
    {
        // Mouse won't be seeable
        window->setMouseCursorVisible(false);
        this->sphericalCamera=false;
        sf::Mouse::setPosition(this->screenCenter, *this->window);
    }
    else
    {
        // Mouse will be seeable
        window->setMouseCursorVisible(true);
        this->sphericalCamera=true;
        this->xMousePosition=0.0;
        this->yMousePosition=0.0;
    }
}


// Prints the last frame duration and corresponding framerate
void Application::printFPS() const
{
    std::cout<<"Frame duration : "<<this->frameDuration.asMilliseconds()<<" ms ("<<1.0/this->frameDuration.asSeconds()<<" FPS)."<<std::endl;
}


//______________________________________________________________________________
// Main loop


// Loops as long as the application is opened
void Application::loop()
{
    sf::Event event;
    sf::Clock globalClock;
    while (!this->done)
    {
        // Checks events
        while (this->window->pollEvent(event))
            this->handleEvent(event);

        // Animates camera and scene
        this->animate(); 
    
        // Updates the simulation
        if (this->simulation!=NULL) this->simulation->update();

        // Renders
        this->renderFrame();

        this->cntFrame++;
        this->frameDuration=this->frameClock.restart();
        this->globalDuration=globalClock.getElapsedTime();
    } 
}


//______________________________________________________________________________
// Application content


// Sets the title
void Application::setTitle(std::string title)
{
    this->title="Fluid simulation";
    std::cout<<std::endl<<std::endl<<"________Fluids Simulation________"<<std::endl;
}


// Sets a test scene (houses, sun...)
void Application::setTestScene()
{
    // Objects creation, building and storage

    // Axis
    /*Object * objectL0=new Object(GL_LINES);
    buildAxis(*objectL0);
    unsigned int storedObjectL0ID=scene.storeObject(objectL0);*/

    // An horizontal square plane
    Object * objectT0=new Object(GL_TRIANGLES);
    buildPlane(*objectT0, 200.0, 200.0);
    unsigned int storedObjectT0ID=scene.storeObject(objectT0);

    // A sphere
    Object * objectT1=new Object(GL_TRIANGLES);
    float radius=0.3; unsigned int discLong=40; unsigned int discLat=20;
    buildSphere_TrSmoothNonRed(*objectT1, radius, discLat, discLong);
    unsigned int storedObjectT1ID=scene.storeObject(objectT1);

    bool smoothObjectFlag;
    std::string fileName;
    // A monkey head (smooth edges)
    /*Object * objectT2=new Object(GL_TRIANGLES);
    smoothObjectFlag=true;
    fileName="../objs/monkeyHead.obj";
    buildObjectGeometryFromOBJ(*objectT2, fileName, smoothObjectFlag);
    unsigned int storedObjectT2ID=scene.storeObject(objectT2);*/

    // A textured house (hard edges)
    Object * objectT3=new Object(GL_TRIANGLES);
    smoothObjectFlag=false;
    fileName="../objs/house.obj";
    buildObjectGeometryFromOBJ(*objectT3, fileName, smoothObjectFlag);
    unsigned int storedObjectT3ID=scene.storeObject(objectT3);


    // Objects to draw

    // One axis
    //unsigned int axisID=scene.addObjectToDraw(storedObjectL0ID); 
       
    // The earth plane
    unsigned int planeID=scene.addObjectToDraw(storedObjectT0ID);
    //float dryEarth[]={0.85, 0.84, 0.80, 1.0};
    float pink[]={1.0, 0.5, 0.5, 1.0};
    scene.setDrawnObjectColor(planeID, pink);
    
    // The sun
    unsigned int sphereID=scene.addObjectToDraw(storedObjectT1ID);
    float T1[16]; float t1[3]={0.0, 10.0, 0.0}; setToTranslate(T1, t1);
    scene.setDrawnObjectModel(sphereID, T1);
    float white[]={1.0, 1.0, 1.0, 1.0};
    scene.setDrawnObjectColor(sphereID, white);
    
    // One Orange monkey head
    /*unsigned int monkeyHeadID=scene.addObjectToDraw(storedObjectT2ID);
    float orange[]={0.7, 0.4, 0.1, 1.0};
    float T2[16]; float t2[3]={0.0, 0.5, 0.5}; setToTranslate(T2, t2);
    float S2[16]; float s2[3]={0.5, 0.5, 0.5}; setToScale(S2, s2);
    multMatrixBtoMatrixA(T2, S2);
    scene.setDrawnObjectModel(monkeyHeadID, T2);
    scene.setDrawnObjectColor(monkeyHeadID, orange);*/

    // A field of textured houses
    unsigned int houseTextureDiffuseID=loadTexture("../textures/house_diffuse.ppm");
    unsigned int houseTextureSpecularID=loadTexture("../textures/house_spec.ppm");
    unsigned int fieldObjectID; float T3[16]; float t3[3];
    int disc=5; float step=4.0/float(disc); float x=-2.0; float z=-2.0;
    for (int iX=-disc ; iX<=disc ; ++iX)
    {
        x=iX*step;
        for (int iZ=-disc ; iZ<=disc ; ++iZ)
        {
            z=iZ*step;
            fieldObjectID=scene.addObjectToDraw(storedObjectT3ID);
            t3[0]=iX; t3[1]=0.5; t3[2]=iZ; setToTranslate(T3, t3);
            scene.setDrawnObjectModel(fieldObjectID, T3);
            scene.setDrawnObjectShader(fieldObjectID, shaders.getLightingTexturingShader());
            scene.setDrawnObjectTextureID(fieldObjectID, 0, houseTextureDiffuseID);
            scene.setDrawnObjectTextureID(fieldObjectID, 1, houseTextureSpecularID);
        }
    }
    std::cout<<std::endl;

    // Errors checking before the loop
    printGLErrors(); std::cout<<std::endl;
}


// Animates an object
void Application::moveScene()
{
    float r=1.0; float s[]={r, r, r};

    float S[16]; setToScale(S, s);

    float d=10.0; float t[]={0.0, d, 0.0};
    float T[16]; setToTranslate(T, t);
    
    float dt=this->globalDuration.asSeconds();
    float rotationalSpeed=M_PI/6.0; // radians per seconds
    float a=rotationalSpeed*dt;
    float axis[]={0.0, 0.0, 1.0};
    float R[16]; setToRotate(R, a, axis); 
    float modelSun[16]; 
    setToIdentity(modelSun);
    multMatrixBtoMatrixA(modelSun, R);
    multMatrixBtoMatrixA(modelSun, T);
    multMatrixBtoMatrixA(modelSun, S);

    unsigned int sunID=1;
    this->scene.setDrawnObjectModel(sunID, modelSun);
    float lightPosition[]={d*modelSun[4], d*modelSun[5], d*modelSun[6], 1.0};
    this->scene.setLight(lightPosition, 1.0);
}


// Sets a simple textured floor
void Application::setSimpleFloor()
{
    // An horizontal square plane
    Object * objectT0=new Object(GL_TRIANGLES);
    buildPlane(*objectT0, 200.0, 200.0);
    unsigned int storedObjectT0ID=scene.storeObject(objectT0);

    // The earth plane
    unsigned int floorID=scene.addObjectToDraw(storedObjectT0ID);

    // Dark wooden
    unsigned int floorTextureDiffuseID=loadTexture("../textures/tileable_wood_texture_01_by_goodtextures-d31qde8Diffuse.ppm", true);
    unsigned int floorTextureSpecularID=loadTexture("../textures/tileable_wood_texture_01_by_goodtextures-d31qde8Specular.ppm", true);

    scene.setDrawnObjectShader(floorID, shaders.getLightingTexturingShader());
    scene.setDrawnObjectTextureID(floorID, 0, floorTextureDiffuseID);
    scene.setDrawnObjectTextureID(floorID, 1, floorTextureSpecularID);

    float l[]={-1.0, 5.0, 0.0, 1.0};
    float lightPower=1.0;
    scene.setLight(l, lightPower);

    /*Object * objectT1=new Object(GL_TRIANGLES);
    float radius=0.1; unsigned int discLong=40; unsigned int discLat=20;
    buildSphere_TrSmoothNonRed(*objectT1, radius, discLat, discLong);
    unsigned int storedObjectT1ID=scene.storeObject(objectT1);
    unsigned int sphereID=scene.addObjectToDraw(storedObjectT1ID);
    float T1[16]; float t1[3]={l[0], l[1], l[2]}; setToTranslate(T1, t1);
    scene.setDrawnObjectModel(sphereID, T1);
    float white[]={1.0, 1.0, 1.0, 1.0};
    scene.setDrawnObjectColor(sphereID, white);*/

}


// Animates a rotational camera
void Application::rotateCamera()
{
    //float dt=this->globalDuration.asSeconds();
    float dt=this->dtSimulation*cntFrame;
    float rotationalSpeed=M_PI/3.0; // radians per seconds
    float a=rotationalSpeed*dt-M_PI/2.0;

    float dist=3.0;
    float c[]={dist*cos(a), 1.0, -dist*sin(a)}; // Camera position    
    float aim[]={0.0, 1.0, 0.0}; // Where we look
    float up[]={0.0, 1.0, 0.0}; // Vector pointing over the camera
    this->camera.lookAt(c, aim, up);
}


// Set the simulation to happen
void Application::setSimulation()
{    
    const float size=3.0;
    const float density=1000.0;
    const float viscosity=0.0;
    unsigned int nbSamplesOneDir=40; 

    const unsigned int nbSamplesX=nbSamplesOneDir;

    //const unsigned int nbSamplesY=nbSamplesOneDir*this->height/this->width;
    const unsigned int nbSamplesY=nbSamplesOneDir/2;

    const unsigned int nbSamplesZ=nbSamplesOneDir/3;

    const unsigned int nbParticlesCoef=2; // On surface, nbSamples*4^nbParticlesCoefs particles 
                              // In volume,  nbSamples*6^nbParticlesCoefs particles
    const bool solidWalls=true;

    this->simulation=new Simulation(shaders, camera, scene, shaders.getLightingShader(), shaders.getSpriteShader(), size, solidWalls, density, viscosity, this->dtSimulation, nbSamplesX, nbSamplesY, nbSamplesZ, nbParticlesCoef);

    unsigned int drawnIDs[10];
    unsigned int i=0;
    drawnIDs[i++]=this->simulation->drawSolids();
    //drawnIDs[i++]=this->simulation->drawSamples();
    //drawnIDs[i++]=this->simulation->drawPressures();
    //drawnIDs[i++]=this->simulation->drawForces();
    //drawnIDs[i++]=this->simulation->drawTypes();
    //drawnIDs[i++]=this->simulation->drawGrid();
    //drawnIDs[i++]=this->simulation->drawVelocitiesBorders();
    //drawnIDs[i++]=this->simulation->drawVelocitiesCenters();
    drawnIDs[i++]=this->simulation->drawParticles();
    //drawnIDs[i++]=this->simulation->drawParticlesVelocities();
    unsigned int n=i;

    float t[]={0.0, 1.0, 0.0};
    float T[16]; setToTranslate(T, t);

    for (unsigned int i=0 ; i<n ; ++i)
        this->scene.setDrawnObjectModel(drawnIDs[i], T);
    this->simulation->setModel(T);
}
