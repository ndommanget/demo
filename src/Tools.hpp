// Tools.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__


// Mersenne Twister random generation
#include <../api/random/MersenneTwister.h>



// Mathematics : vectors
float distance(const float * const A, const float * const B, bool square);
float getNorm(const float * const a);
void normalize(float * const a);
float dotProduct(const float * const a, const float * const b);
void vectorProduct (const float * const a, const float * const b, float * const result);
float scalarTriple(const float * const a, const float * const b, const float * const c);
void normalFace(const float * const O, const float * const A, const float * const B, float * const normal, bool toNormalize);

// Mathematics : matrices
void setToIdentity(float * const matrix);
void setToTranslate(float * const matrix, const float * const t);
void setToScale(float * const matrix, const float * const s);
void setToRotate(float * const matrix, float angle, const float * const axis);
void setPerspective(float * const mat, float l, float r, float b, float t, float n, float f);
void multMatrixBtoMatrixA(float * const A, const float * const B);
void multMatrixToPoint(float * const pos, const float * const mat);
void getInverseMatrix(const float * const M, float * const Mm1);
void getTransposedMatrix(const float * const M, float * const Mt);


// Mathematics : interpolation
void linearInterpolation(unsigned int nbCoords, const float * const x0, const float * const x1, float toX0, float * const result);
void biLinearInterpolation(unsigned int nbCoords, const float * const x0y0, const float * const x1y0, const float * const x0y1, const float * const x1y1, float toX0, float toY0, float * const result);
void triLinearInterpolation(unsigned int nbCoords, const float * const x0y0z0, const float * const x1y0z0, const float * const x0y1z0, const float * const x1y1z0, const float * const x0y0z1, const float * const x1y0z1, const float * const x0y1z1, const float * const x1y1z1, float toX0, float toY0, float toZ0, float * const result);

// Random generation
double getRand(MTRand & mtrand, bool uniform, double min, double max);

// Printing
void printVec2(const float * const vect);
void printVec3(const float * const vect);
void printVec4(const float * const vect);
void printMat16(const float * const mat);
void printGLErrors();

// Texture loading
unsigned int loadTexture(const std::string &fileName, bool mipmapping=true);

// Images reading and writing
unsigned char * loadPPM(const std::string & filename, unsigned int * const width, unsigned int * const height);
void saveFrameBufferPPM(const std::string & fileName, unsigned int width, unsigned int height);


#endif //  __TOOLS_HPP__
