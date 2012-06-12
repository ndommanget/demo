// Object.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __OBJECT_HPP__
#define __OBJECT_HPP__


// OpenGL
#include "glew/glew.h"
#include <SFML/OpenGL.hpp>



// An object made of solid triangles
class Object
{
    private:

        // OpenGL IDs
	unsigned int vaoId;          // Vertex array objet ID
	unsigned int vboId;          // Vertices buffer ID
	unsigned int normalsVboId;   // Normals buffer ID
	unsigned int uvsVboId;       // Normals buffer ID
	unsigned int colorsVboId;    // Normals buffer ID
	unsigned int indicesVboId;   // Indices buffer ID

        // Size
	unsigned int nbVertices;     // Defined vertices number
	unsigned int nbIndices;      // Indices number

        // Filled informations
        bool primitives;             // Are the primitives filled ?
        bool normals;                // Are the normals filled ?
        bool colors;                 // Are the colors filled ?
        bool uvs;                    // Are the uv coordinates filled ?

        // Primitive type
	const GLenum primitivesType; // GL_TRIANGLES or GL_LINES



        // Initialization
        void init();

    public:

        // Constructors / destructor
        Object(const GLenum primitivesType=GL_TRIANGLES);
        ~Object();

        // Setters / getters
        void setNbVertices(const unsigned int nbVertices);
        void setNbIndices(const unsigned int nbIndices);
        unsigned int getNbVertices() const;
        unsigned int getNbIndices() const;
        bool hasPrimitives() const;
        bool hasNormals() const;
        bool hasUvs() const;
        bool hasColors() const;

        // Sending data on GPU
        void sendPrimitives(const float * const vertices, const unsigned int * const indices, bool dynamicVertices=false, bool dynamicIndices=false);
        void sendNormals(const float * const  normals, bool dynamic=false);
        void sendUvs(const float * const  uvs, bool dynamic=false);
        void sendColors(const float * const  colors, bool dynamic=false);

        // Updating data on GPU
        void updateVertices(const float * const vertices, bool dynamic=false) const;
        void updateIndices(const unsigned int * const indices, bool dynamic=false) const;
        void updateNormals(const float * const normals, bool dynamic=false) const;
        void updateUvs(const float * const uvs, bool dynamic=false) const;
        void updateColors(const float * const colors, bool dynamic=false) const;       

        // Drawing
        void drawObject() const;
};


#endif //__OBJECT_HPP__
