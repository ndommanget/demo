// Object.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Object.hpp"
#include "Tools.hpp"



//______________________________________________________________________________
// Initialization


// Inits parameters for the object, 
// and values necessary for the storage of the object on the GPU
void Object::init()
{
    this->nbVertices=0;
    this->nbIndices=0;

    this->primitives=false;
    this->normals=false;
    this->uvs=false;
    this->colors=false;

	// Creation of ids for the buffers on GPU.
	// We store them in the structure for clarity
    	// Creates a VAO id to handle the vao for objectTr
    glGenVertexArrays(1, &(this->vaoId));
    	// Creates a VBO id for a VBO to store the vertices
    glGenBuffers(1, &(this->vboId));
    	// Creates a VBO id for a VBO to store the normals
    glGenBuffers(1, &(this->normalsVboId));
    	// Creates a VBO id for a VBO to store the colors
    glGenBuffers(1, &(this->colorsVboId)); 
    	// Creates a VBO id for a VBO to store the uv coordinates (for textures)
    glGenBuffers(1, &(this->uvsVboId));        
    	// Creates a VBO id for a VBO to store the indices of the vertices
    glGenBuffers(1, &(this->indicesVboId));
}


//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Default constructor
Object::Object(const GLenum primitivesType) : primitivesType(primitivesType)
{
    this->init();
}


// Cleans memory for Object
Object::~Object()
{
    // Cleans by deleting the buffers
    glDeleteBuffers(1, &(this->vboId));
    glDeleteBuffers(1, &(this->normalsVboId));
    glDeleteBuffers(1, &(this->colorsVboId));
    glDeleteBuffers(1, &(this->uvsVboId));
    glDeleteBuffers(1, &(this->indicesVboId));
    glDeleteVertexArrays(1, &(this->vaoId));
}


//______________________________________________________________________________
// Setters / getters


// Sets the number of vertices
void Object::setNbVertices(const unsigned int nbVertices)
{
    this->nbVertices=nbVertices;
}


// Sets the number of indices
void Object::setNbIndices(const unsigned int nbIndices)
{
    this->nbIndices=nbIndices;
}


// Returns the number of vertices
unsigned int Object::getNbVertices() const
{
    return nbVertices;
}


// Returns the number of indices
unsigned int Object::getNbIndices() const
{
    return nbIndices;
}


// Returns if primitives were provided
bool Object::hasPrimitives() const
{
    return this->primitives;
}


// Returns if normals were provided
bool Object::hasNormals() const
{
    return this->normals;
}


// Returns if uv coordinates were provided
bool Object::hasUvs() const
{
    return this->uvs;
}


// Returns if vertex colors were provided
bool Object::hasColors() const
{
    return this->colors;
}



//______________________________________________________________________________
// Sending data on GPU


// Sends an array of vertices and an array of indices in buffers on the GPU
void Object::sendPrimitives(const float * const vertices, const unsigned int * const indices, bool dynamicVertices, bool dynamicIndices)
{
    // Binds the vao
    glBindVertexArray(this->vaoId);
    
    // Deactivates the generic attributes in shader 
    // (specific with the use of our shader tool)
    GLenum attributePosition=0; // First attribute in the shader is position
    glEnableVertexAttribArray(attributePosition);

    // Binds the vbo "objectTr->vboId" as a buffer for vertices
    glBindBuffer(GL_ARRAY_BUFFER, this->vboId);
    // if "dynamic", fills the vbo with the colors data, in a GPU memory special for "dynamic draw" access, else "static draw" (data should be used to draw frequently, but not be changed often)
    if (dynamicVertices) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), vertices, GL_DYNAMIC_DRAW);
    else                 glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), vertices, GL_STATIC_DRAW);
    // Specifies from where and how to read the data in the bound buffer
    glVertexAttribPointer(attributePosition, 4, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);

    if (this->primitivesType!=GL_POINTS)
    {
        // Binds the other vbo
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->indicesVboId); 
        // Loads up the indices of the vertices
        if (dynamicIndices) glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nbIndices*sizeof(unsigned int), indices, GL_DYNAMIC_DRAW);
        else                glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nbIndices*sizeof(unsigned int), indices, GL_STATIC_DRAW);
    }
    
    // Unbinds the vao
    glBindVertexArray(0);
    
    this->primitives=true;
}


// Sends an array of normals in buffers on the GPU
void Object::sendNormals(const float * const normals, bool dynamic)
{
    // Binds the vao
    glBindVertexArray(this->vaoId);
    
    // Activates the generic attribute of normal in shader
    GLenum attributeNormal=1; // Second attribute in the shader is the normal
    glEnableVertexAttribArray(attributeNormal);

    // Binds the vbo "object->normalsVboId" as a buffer for normals
    glBindBuffer(GL_ARRAY_BUFFER, this->normalsVboId);
    // if "dynamic", fills the vbo with the colors data, in a GPU memory special for "dynamic draw" access, else "static draw" (data should be used to draw frequently, but not be changed often)
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*3*sizeof(float), normals, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*3*sizeof(float), normals, GL_STATIC_DRAW);
    // Specifies from where and how to read the data in the bound buffer
    glVertexAttribPointer(attributeNormal, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);

    // Unbinds the vao
    glBindVertexArray(0);
    
    this->normals=true;
}


// Sends an array of colors in buffers on the GPU
void Object::sendUvs(const float * const uvs, bool dynamic)
{
     // Binds the vao
    glBindVertexArray(this->vaoId);
    
    // Activates the generic attribute of normal in shader     
    GLenum attributeUvs=2; // Third attribute in the shader is the texture coordinates
    glEnableVertexAttribArray(attributeUvs);

    // Binds the vbo "object->colorsVboId" as a buffer for normals
    glBindBuffer(GL_ARRAY_BUFFER, this->uvsVboId);
    // if "dynamic", fills the vbo with the colors data, in a GPU memory special for "dynamic draw" access, else "static draw" (data should be used to draw frequently, but not be changed often)
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*2*sizeof(float), uvs, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*2*sizeof(float), uvs, GL_STATIC_DRAW);
    // Specifies from where and how to read the data in the bound buffer
    glVertexAttribPointer(attributeUvs, 2, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);

    // Unbinds the vao
    glBindVertexArray(0);
    
    this->uvs=true;
}


// Sends an array of colors in buffers on the GPU
void Object::sendColors(const float * const colors, bool dynamic)
{
    // Binds the vao
    glBindVertexArray(this->vaoId);
    
    // Activates the generic attribute of normal in shader
    GLenum attributeColor=3; // Fourth attribute in the shader is the color
    glEnableVertexAttribArray(attributeColor);

    // Binds the vbo "object->colorsVboId" as a buffer for normals
    glBindBuffer(GL_ARRAY_BUFFER, this->colorsVboId);
    // if "dynamic", fills the vbo with the colors data, in a GPU memory special for "dynamic draw" access, else "static draw" (data should be used to draw frequently, but not be changed often)
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), colors, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), colors, GL_STATIC_DRAW);
    // Specifies from where and how to read the data in the bound buffer
    glVertexAttribPointer(attributeColor, 4, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);

    // Unbinds the vao
    glBindVertexArray(0);
    
    this->colors=true;
}


//______________________________________________________________________________
// Updating data on GPU


// Updates vertices positions on GPU
void Object::updateVertices(const float * const vertices, bool dynamic) const
{
    glBindBuffer(GL_ARRAY_BUFFER, this->vboId);
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), vertices, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), vertices, GL_STATIC_DRAW);
}


// Updates primitives indices on GPU
void Object::updateIndices(const float * const indices, bool dynamic) const
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->indicesVboId); 
    if (dynamic) glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nbIndices*sizeof(unsigned int), indices, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nbIndices*sizeof(unsigned int), indices, GL_STATIC_DRAW);
}


// Updates normals on GPU
void Object::updateNormals(const float * const normals, bool dynamic) const
{
    glBindBuffer(GL_ARRAY_BUFFER, this->normalsVboId);
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*3*sizeof(float), normals, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*3*sizeof(float), normals, GL_STATIC_DRAW);
}


// Updates Uvs on GPU
void Object::updateUvs(const float * const uvs, bool dynamic) const
{
    glBindBuffer(GL_ARRAY_BUFFER, this->uvsVboId);
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*2*sizeof(float), uvs, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*2*sizeof(float), uvs, GL_STATIC_DRAW);
}


// Updates colors on GPU
void Object::updateColors(const float * const colors, bool dynamic) const
{
    glBindBuffer(GL_ARRAY_BUFFER, this->colorsVboId);
    if (dynamic) glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), colors, GL_DYNAMIC_DRAW);
    else         glBufferData(GL_ARRAY_BUFFER, this->nbVertices*4*sizeof(float), colors, GL_STATIC_DRAW);
} 


//______________________________________________________________________________
// Drawing


// Draw the object with the vao
void Object::drawObject() const
{
	// Binds the vao to draw (the vbos binds are stored in the vao, no need to redo them)
    glBindVertexArray(this->vaoId);
    
    // Draw the elements stored by the vao
    if (this->primitivesType==GL_POINTS)
        glDrawArrays(this->primitivesType, 0, this->nbVertices);
    else      
        glDrawElements(this->primitivesType, this->nbIndices, GL_UNSIGNED_INT, 0);
}
