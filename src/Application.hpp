// Application.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __APPLICATION_HPP__
#define __APPLICATION_HPP__


#include "Scene.hpp"
#include "Simulation.hpp"

// SFML
#include <SFML/Window.hpp>



// All the initialisation of states and events for SDL and OpenGL
class Application
{
    private:

        // Window
        sf::Window * window;            // SFML Window, GL context info
        std::string title;              // Window title
        bool done;                      // True when window closed

        // Window Sizes
	unsigned int windowedWidth;     // Rendering region size on x (when not fullscreen)
	unsigned int windowedHeight;    // Rendering region size on y (when not fullscreen)
        unsigned int width;             // Actual rendering region on x
        unsigned int height;            // Actual rendering region on y
        unsigned int minWindowSize;     // Minimum rendering region size
        bool fullScreen;                // Application in full screen mode

        // Time data
        sf::Time globalDuration;        // global time duration (since loop)
        sf::Clock frameClock;           // Clock to measure frame durations
        sf::Time frameDuration;         // Frame time duration
        unsigned int cntFrame;          // Frames counter
        float dtSimulation;             // Delta of time between simulation updates

        // Mouse data
        sf::Vector2i screenCenter;      // Center position of the screen
        float xMousePosition;           // Mouse position on x
        float yMousePosition;           // Mouse position on y
        float xMouseOffset;             // Offset for mouse move on x
        float yMouseOffset;             // Offset for mouse move on y
        float xMouseClick;              // Mouse position on click on x
        float yMouseClick;              // Mouse position on click on y          
        int scroll;                     // Mouse scroll value
        bool mouse;                     // True if mouse is seeable

        // Camera data
        bool sphericalCamera;           // Spherical camera (false : FPS)
        float speed;
        float moveTimes[6];             // Elapsed time until ctl key pressed
        bool movePressed[6];            // Ctl key pressed (x,-x,y,-y,z,-z)

        // Additional data
        Shaders shaders;                // Shaders to render the scene
        Camera camera;                  // Camera to watch the scene
        Scene scene;                    // Scene to draw
        Simulation * simulation;        // Simulation to update and draw



        // Initialization
        void init();
        void initSFMLOpenGL();
        void customizeStates();

        // User interface
        void handleEvent(const sf::Event & event);
        void handleKeyEvent(const sf::Event & event, bool down);

        // Camera animation
        void moveCamera();

        // Main steps
        void animate();
        void storeFrame() const;
        void renderFrame() const;

    public:

        // Constructors / destructor
        Application();
        ~Application();

        // Application modes
        void switchWireframe();
        void resize();
        void switchFullScreen();
        void switchCameraMode();
        void printFPS() const;

        // Main loop
        void loop();

        // Application content
        void setTitle(std::string title);
        void setTestScene();
        void moveScene();
        void setSimpleFloor();
        void rotateCamera();
        void setSimulation(); 
};


#endif //__APPLICATION_HPP__ 
