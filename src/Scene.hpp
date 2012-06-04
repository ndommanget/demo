// Scene.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __SCENE_HPP__
#define __SCENE_HPP__


#include "Camera.hpp"
#include "Shaders.hpp"
#include "Object.hpp"



// A container for objects
class Scene
{
    private:

        // Scene maximum size
        unsigned int maxStoredObjects;          // Maximum number of storable objects
        unsigned int maxDrawnObjects;           // Maximum number of drawable objects
        
        // Default drawing set-up
        float defaultColor[4];                  // Default color
        float defaultModel[16];                 // Default transformation matrix
        unsigned int defaultShader;             // Default shader
        unsigned int defaultTextureID;          // Default textureID

        // Light
        float lightPosition[4];                 // Position of the light used in shader
        float lightPower;                       // Power of the light used in shader
        
        // Shaders
        Shaders & shaders;                      // Shaders used to render the scene
        // Camera
        Camera & camera;                        // Camera used to watch the scene


        // Library of potential objects
        const Object ** storedObjects;          // Library of Objects to use from GPU
        unsigned int nbStoredObjects;           // Number of stored objects
        
        // Drawn instances of objects
        unsigned int * drawnObjects;            // Drawn object library indices
        float * drawnObjectsColors;             // Drawn object colors
        float * drawnObjectsModels;             // Drawn object transformation matrices
        unsigned int * drawnObjectsShaders;     // Drawn object shader ID 
        unsigned int * drawnObjectsTexture0IDs; // Drawn object shader ID on unit 0
        unsigned int * drawnObjectsTexture1IDs; // Drawn object shader ID on unit 1
        unsigned int nbDrawnObjects;


        
        // Drawing preparation
        void setAppearance(unsigned int indexDrawnObject) const;

    public:

        // Constructors / destructor
        Scene(Shaders & shaders, Camera & camera);
        ~Scene();

        // Initialization
        void init();

        // Building
        unsigned int storeObject(const Object * const object);
        unsigned int addObjectToDraw(unsigned int indexStoredObject);

        // Default set-up
        void setDefaultColor(const float * const defaultColor);
        void setDefaultModel(const float * const defaultModel);
        void setDefaultShader(unsigned int defaultShader);
        void setDefaultTextureID(unsigned int defaultTextureID);

        // Drawn objects set-up
        void setDrawnObjectColor(unsigned int indexDrawnObject, const float * const color);
        void setDrawnObjectModel(unsigned int indexDrawnObject, const float * const model);
        void setDrawnObjectShader(unsigned int indexDrawnObject, unsigned int shader);
        void setDrawnObjectTextureID(unsigned int indexDrawnObject, unsigned int textureUnit, unsigned int textureID);  

        // Light set-up
        void setLight(const float * const position, float power);

        //Drawing
        void draw() const;
};


#endif // __SCENE_HPP__
