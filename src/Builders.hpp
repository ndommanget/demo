// Builders.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __BUILDERS_HPP__
#define __BUILDERS_HPP__


// C++
#include <vector>
#include <string>



//Forward declaration
class Object;

// Building functions
void buildAxis(Object & object);
void buildPlane(Object & object, float width, float length);
void buildSphere_TrSmoothNonRed(Object & object, float radius, unsigned int discLat, unsigned int discLong);
void buildWindRose(Object & object);
void buildCube(Object & object);

// Mesh conformation functions
void setNormalsFlatTr(const Object & object, const float * const vertices, const unsigned int * const indices, float * const normals);
void centerAndNormalizeMesh(const Object & object, float * const vertices);


// Obj file reading fuctions
void split(const std::string * string, std::vector<std::string> * tokens, const std::string & delim);
void readVec3(std::istringstream & line, std::vector<float> * const vertices);
void readVec2(std::istringstream & line, std::vector<float> * const vertices);
void readFace(std::istringstream & line, std::vector<unsigned int> * const indices,  std::vector<unsigned int> * const uvIndices, std::vector<unsigned int> * const normalsIndices);
void reorderUvsAndNormalsIfSmooth(const std::vector<float> * const vertices, std::vector<float> * const uvs, std::vector<float> * const normals, const std::vector<unsigned int> * const indices, const std::vector<unsigned int> * const uvIndices, const std::vector<unsigned int> * const normalIndices);
void reorderUvsAndNormalsIfNonSmooth(std::vector<float> * const vertices, std::vector<float> * const uvs, std::vector<float> * const normals, const std::vector<unsigned int> * const indices, const std::vector<unsigned int> * const uvIndices, const std::vector<unsigned int> * const normalsIndices);
void addHomogeneousToVertices(std::vector<float> * const vertices);
void conformToObject(std::vector<float> * const vertices, const std::vector<float> * const uvs, std::vector<float> * const normals);
bool buildObjectGeometryFromOBJ(Object & object, const std::string& fileName, bool smoothObject);


#endif
