// Shaders.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Shaders.hpp"
#include "Tools.hpp"

// Stream
#include <fstream>
#include <sstream>



//______________________________________________________________________________
// Compilation


// Writes content of a text file [fileName] in a returned std::string
std::string * Shaders::loadFile(const std::string &fileName) const
{
    std::string* result = new std::string();
    std::ifstream file(fileName.c_str());
    if (!file) 
    {
        std::cerr<<"Cannot open file"<<fileName<<std::endl;
        throw std::exception();
    }
    std::string line;
    while (getline(file, line)) 
    {
        *result+=line;
        *result+='\n';
    }
    file.close();
    return result;
}


// Prints info about potential problems which occured during compilation
void Shaders::printShaderLog(unsigned int shaderId) const
{
    int logLength;
    glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength>0)
    {
        char log[logLength];
        glGetShaderInfoLog(shaderId, logLength, &logLength, log);
        std::cout << std::string(log);
    }
}


// Loads files [files] in string table, builds program object, compile shaders in shaders objects, links and return the id program object (with the executable code)
unsigned int Shaders::loadProgram(const std::vector<std::string> &files) const
{
    // Creates a program object which id is returned
    unsigned int programId=glCreateProgram();

    glBindAttribLocation(programId, 0, "vertexPosition");
    glBindAttribLocation(programId, 1, "vertexNormal");
    glBindAttribLocation(programId, 2, "vertexUvs");
    glBindAttribLocation(programId, 3, "vertexColor");

    // Creates a vertex shader object which id is returned
    unsigned int vertexShaderId=glCreateShader(GL_VERTEX_SHADER);
    // Creates a fragment shader object which id is returned
    unsigned int fragmentShaderId=glCreateShader(GL_FRAGMENT_SHADER);
    glAttachShader(programId, vertexShaderId);
    glAttachShader(programId, fragmentShaderId);

    unsigned int n=files.size();
    std::string ** strs=new std::string*[n];
    const char ** lines=new const char*[n+2];
    //std::cout<<"Loading program "<<files[n - 1]<<"..."<<std::endl;
    
    bool geo=false;
    // For every file :
    for (unsigned int i=0 ; i<n ; ++i)
    {
        // Gets std::string of file
        std::string* s=loadFile(files[i]);
        // Stores it [std::string] in strs[i]
        strs[i]=s;
        // Stores it [char] in lines[i+2]
        lines[i+2]=s->c_str();
        
        // If _GEOMETRY_ is in file, geometrey shader : geo=true
        // strstr(a, b)-> finds firt occurence of b in a
        if (strstr(lines[i+2], "_GEOMETRY_")!=NULL) 
            geo=true;
    }

    lines[0]="#version 150\n";
    lines[1]="#define _VERTEX_\n";
    // Loads shader source chars in shader object
    glShaderSource(vertexShaderId, n+2, lines, NULL);
    // Compiles the loaded shader source code
    glCompileShader(vertexShaderId);
    // Prints compilation potential problems
    printShaderLog(vertexShaderId);

    if (geo) 
    {
        // Creates a geometry shader object which id is returned
        unsigned int geometryShaderId=glCreateShader(GL_GEOMETRY_SHADER_EXT);
        glAttachShader(programId, geometryShaderId);
        // Adds this text before the source text
        lines[1]="#define _GEOMETRY_\n";
        glShaderSource(geometryShaderId, n+2, lines, NULL);
        glCompileShader(geometryShaderId);
        printShaderLog(geometryShaderId);
        // Specifies type of primitives accessible in geometry shader
        glProgramParameteriEXT(programId, GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES);
    }

    lines[1]="#define _FRAGMENT_\n";
    glShaderSource(fragmentShaderId, n+2, lines, NULL);
    glCompileShader(fragmentShaderId);
    printShaderLog(fragmentShaderId);

    // Links the program object to build the executable code
    glLinkProgram(programId);

    for (unsigned int i=0 ; i<n ; ++i) 
        delete strs[i];

    delete[] strs;
    delete[] lines;

    // Returns the ID of the program object
    return programId;
}


//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Default constructor
Shaders::Shaders()
{

}


// Cleans memory for Shaders
Shaders::~Shaders()
{

}


//______________________________________________________________________________
// Initializes


// Compiles, stores and inializes the shaders (needs a GL context)
void Shaders::init()
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
    {
        this->shaderIDs[iShaders]=0;
        this->transparent[iShaders]=false;
    }

    // Uniforms common basis
    std::vector<std::string> files;
    files.push_back("../shaders/shaderTools.glsl");

    // Regular lighting shader
    files.push_back("../shaders/lightingShader.glsl");
    this->shaderIDs[0]=loadProgram(files); files.pop_back();
    
    // Texture and lighting shader
    files.push_back("../shaders/lightingTexturingShader.glsl");
    this->shaderIDs[1]=loadProgram(files); files.pop_back();
    
    // Sprites shader
    files.push_back("../shaders/spriteShader.glsl");
    this->shaderIDs[2]=loadProgram(files); files.pop_back();
    this->transparent[2]=true;


    // Initializes the uniforms locations in shaders
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
    {
        this->materialAmbientLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.ambient");
        this->materialDiffuseLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.diffuse");
        this->materialSpecularLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.specular");
        this->materialKaLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.ka");
        this->materialKdLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.kd");
        this->materialKsLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.ks");
        this->materialShininessLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "material.shininess");
        this->lightPositionLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "light.position");
        this->lightPowerLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "light.power");
        this->modelLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "model");
        this->cameraViewLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "view");
        this->cameraCLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "center");
        this->cameraProjectionLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "projection");
        this->filledDataLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "filledData");
        this->textureUnitDiffuseLocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "textureUnitDiffuse");
        this->textureUnitSpecularlocs[iShaders]=glGetUniformLocation(shaderIDs[iShaders], "textureUnitSpecular");
    }

    // Main material settings
    float ambient[]={1.0, 1.0, 1.0, 1.0}; float ka=0.01;
    float diffuse[]={1.0, 1.0, 1.0, 1.0}; float kd=1.0;
    float specular[]={1.0, 1.0, 1.0, 1.0}; float ks=2.0; float shininess=5.0;
    this->setMaterial(ambient, diffuse, specular, ka, kd, ks, shininess);
     
    // Texture settings
    this->setTextureUnits(1);
    this->setTextureUnits(2);
}


//______________________________________________________________________________
// Setters / getters


// Is the shader transprent ?
bool Shaders::isTransparent(unsigned int iShader) const
{
    return transparent[iShader];
}


// Provides the lighting shader
unsigned int Shaders::getLightingShader() const
{
    return 0;
}


// Provides the lighting texturing shader
unsigned int Shaders::getLightingTexturingShader() const
{
    return 1;
}


// Provides the sprite shader
unsigned int Shaders::getSpriteShader() const
{
    return 2;
}


//______________________________________________________________________________
// Data settings in all the shaders


// Passes the light to all the shaders
void Shaders::setLight(const float * const position, float power) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setLight(iShaders, position, power);
}

// Passes the material to all the shaders
void Shaders::setMaterial(const float * const ambient, const float * const diffuse, const float * const specular, float ka, float kd, float ks, float shininess) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setMaterial(iShaders, ambient, diffuse, specular, ka, kd, ks, shininess);
}


// Passes the units of the textures to use for diffuse and specular as a uniform samplers to all the shaders
void Shaders::setTextureUnits() const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setTextureUnits(iShaders);
}


// Passes the model matrix to all the shaders
void Shaders::setModel(const float * const model) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setModel(iShaders, model);
}


// Passes the view matrix to all the shaders
void Shaders::setView(const float * const view, const float * const c) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setView(iShaders, view, c);
}


// Passes the projection matrix to all the shaders
void Shaders::setProjection(const float * const projection) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->setProjection(iShaders, projection);
}


// Changes color in the material in all the shaders (diffuse component)
void Shaders::changeMaterialColor(const float * const color) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
        this->changeMaterialColor(iShaders, color);
}


// Passes filledData to all the shader (flag saying which data is provided as attribute)
void Shaders::setFilledData(bool primitives, bool normals, bool uvs, bool colors) const
{
    for (unsigned int iShaders=0 ; iShaders<3 ; ++iShaders)
       this->setFilledData(iShaders, primitives, normals, uvs, colors);
}


//______________________________________________________________________________
// Data settings in one shader


// Passes the light to one shader
void Shaders::setLight(unsigned int shader, const float * const position, float power) const
{
    glUseProgram(shaderIDs[shader]);
    glUniform4fv(this->lightPositionLocs[shader], 1, position);
    glUniform1f(this->lightPowerLocs[shader], power);
}


// Passes the material to one shader
void Shaders::setMaterial(unsigned int shader, const float * const ambient, const float * const diffuse, const float * const specular, float ka, float kd, float ks, float shininess) const
{
    glUseProgram(shaderIDs[shader]);
    glUniform4fv(this->materialAmbientLocs[shader], 1, ambient);
    glUniform4fv(this->materialDiffuseLocs[shader], 1, diffuse);
    glUniform4fv(this->materialSpecularLocs[shader], 1, specular);
    glUniform1f(this->materialKaLocs[shader], ka);
    glUniform1f(this->materialKdLocs[shader], kd);
    glUniform1f(this->materialKsLocs[shader], ks);
    glUniform1f(this->materialShininessLocs[shader], shininess);
}


// Passes the units of the textures to use for diffuse and specular as a uniform samplers to one shader
void Shaders::setTextureUnits(unsigned int shader) const
{
    glUseProgram(shaderIDs[shader]);
    glUniform1i(this->textureUnitDiffuseLocs[shader], 0);
    glUniform1i(this->textureUnitSpecularlocs[shader], 1);
}


// Passes the model matrix to one shader
void Shaders::setModel(unsigned int shader, const float * const model) const
{
    bool toTranspose=false;
    glUseProgram(shaderIDs[shader]);
    glUniformMatrix4fv(this->modelLocs[shader], 1, toTranspose, model);
}


// Passes the view matrix to one shader
void Shaders::setView(unsigned int shader, const float * const view, const float * const c) const
{
    bool toTranspose=false;
    glUseProgram(shaderIDs[shader]);
    glUniformMatrix4fv(this->cameraViewLocs[shader], 1, toTranspose, view);
    glUniform4fv(this->cameraCLocs[shader], 1, c);
}


// Passes the projection matrix to one shader
void Shaders::setProjection(unsigned int shader, const float * const projection) const
{
    bool toTranspose=false;
    glUseProgram(shaderIDs[shader]);
    glUniformMatrix4fv(this->cameraProjectionLocs[shader], 1, toTranspose, projection);
}


// Changes color in the material in one shader (diffuse component)
void Shaders::changeMaterialColor(unsigned int shader, const float * const color) const
{
    glUseProgram(shaderIDs[shader]);
    glUniform4fv(this->materialDiffuseLocs[shader], 1, color);
}


// Passes filledData to one shader (flag saying which data is provided as attribute)
void Shaders::setFilledData(unsigned int shader, bool primitives, bool normals, bool uvs, bool colors) const
{
    glUseProgram(shaderIDs[shader]);
    glUniform4i(this->filledDataLocs[shader], primitives, normals, uvs, colors);
}
