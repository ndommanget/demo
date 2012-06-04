// Shaders.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __SHADERS_HPP__
#define __SHADERS_HPP__


// OpenGL
#include "glew/glew.h"
#include <SFML/OpenGL.hpp>

// C++
#include <string>
#include <vector>



// Shaders to render
class Shaders
{
    private:

        // General data
        unsigned int shaderIDs[3];              // Shaders OpenGL IDs
        bool transparent[3];                    // Transparent shader ?

        // Uniforms locations in shaders
        int materialAmbientLocs[3];             // Material ambient color location
        int materialDiffuseLocs[3];             // Material diffuse color location
        int materialSpecularLocs[3];            // Material specular color location
        int materialKaLocs[3];                  // Material ambient coef location
        int materialKdLocs[3];                  // Material diffuse coef location
        int materialKsLocs[3];                  // Material specular coef location
        int materialShininessLocs[3];           // Material shininess coef location
        int lightPositionLocs[3];               // Light position location
        int lightPowerLocs[3];                  // Light power location
        int modelLocs[3];                       // Model matrix location 
        int cameraViewLocs[3];                  // Camera view matrix location
        int cameraCLocs[3];                     // Camera center location
        int cameraProjectionLocs[3];            // Camera projection matrix location         
        int filledDataLocs[3];                  // Filled data flags location
        int textureUnitDiffuseLocs[3];          // Diffuse texture unit location
        int textureUnitSpecularlocs[3];         // Specular texture unit location



        // Compilation
        std::string * loadFile(const std::string & fileName) const;
        void printShaderLog(unsigned int shaderId) const;
        unsigned int loadProgram(const std::vector<std::string> & files) const;

    public:

        // Constructors / destructor
        Shaders();
        ~Shaders();

        // Initialization
        void init();

        // Setters / getters
        bool isTransparent(unsigned int iShader) const;
        unsigned int getLightingShader() const;
        unsigned int getLightingTexturingShader() const;
        unsigned int getSpriteShader() const;

        // Data settings in all the shaders
        void setLight(const float * const position, float power) const;
        void setMaterial(const float * const ambient, const float * const diffuse, const float * const specular, float ka, float kd, float ks, float shininess) const;
        void setTextureUnits() const;
        void setModel(const float * const model) const;
        void setView(const float * const view, const float * const c) const;
        void setProjection(const float * const projection) const;
        void changeMaterialColor(const float * const color) const;
        void setFilledData(bool primitives, bool normals, bool uvs, bool colors) const;

        // Data settings in one shader
        void setLight(unsigned int iShaders, const float * const position, float power) const;
        void setMaterial(unsigned int iShaders, const float * const ambient, const float * const diffuse, const float * const specular, float ka, float kd, float ks, float shininess) const;
        void setTextureUnits(unsigned int iShaders) const;
        void setModel(unsigned int iShaders, const float * const model) const;
        void setView(unsigned int iShaders, const float * const view, const float * const c) const;
        void setProjection(unsigned int iShaders, const float * const projection) const;
        void changeMaterialColor(unsigned int iShaders, const float * const color) const;
        void setFilledData(unsigned int iShaders, bool primitives, bool normals, bool uvs, bool colors) const;
};


#endif // __SHADERS_HPP__
