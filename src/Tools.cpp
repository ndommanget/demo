// Tools.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Tools.hpp"


// OpenGL
#include "glew/glew.h"
#include <SFML/OpenGL.hpp>

// Stream
#include <fstream>
#include <sstream>



//______________________________________________________________________________
// Mathematics : vectors


// To get the norm of a vector
float getNorm (const float * const a)
{
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}


// To normalize a vector
void normalize (float * const a)
{
	float norm=getNorm(a);
	if (norm!=0.0)
	{
	    a[0]/=norm;
	    a[1]/=norm;
	    a[2]/=norm;
	}
}


// To get the vector product
void vectorProduct (const float * const a, const float * const b, float * const result)
{
	result[0]=a[1]*b[2]-a[2]*b[1];
	result[1]=a[2]*b[0]-a[0]*b[2];
	result[2]=a[0]*b[1]-a[1]*b[0];
}


// To fill the provided normal
void normalFace(const float * const O, const float * const A, const float * const B, float * const normal, bool toNormalize)
{	
	float a[]={A[0]-O[0], A[1]-O[1], A[2]-O[2]};
	float b[]={B[0]-O[0], B[1]-O[1], B[2]-O[2]};		

	vectorProduct(&(a[0]), &(b[0]), normal);
	
    if (toNormalize)
        normalize(normal);
}


//______________________________________________________________________________
// Mathematics : matrices


// Sets the provided matrix to identity
void setToIdentity(float * const matrix)
{
    float I[]={1.0, 0.0, 0.0, 0.0, 
               0.0, 1.0, 0.0, 0.0, 
               0.0, 0.0, 1.0, 0.0, 
               0.0, 0.0, 0.0, 1.0};
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        matrix[iMatrixCoord]=I[iMatrixCoord];
}


// Sets the provided matrix to a translate matrix on vector t
void setToTranslate(float * const matrix, const float * const t)
{
    float T[]={1.0,   0.0,   0.0,   0.0,
               0.0,   1.0,   0.0,   0.0,
               0.0,   0.0,   1.0,   0.0,
               t[0],  t[1],  t[2],  1.0}; 
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        matrix[iMatrixCoord]=T[iMatrixCoord];
}


// Sets the provided matrix to a scale matrix by coeficients in s
void setToScale(float * const matrix, const float * const s)
{
    float S[]={s[0], 0.0,  0.0,  0.0,
               0.0,  s[1], 0.0,  0.0,
               0.0,  0.0,  s[2], 0.0,
               0.0,  0.0,  0.0,  1.0};  
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        matrix[iMatrixCoord]=S[iMatrixCoord];
}


// Sets the provided matrix to a rotate matrix of angle "angle", around axis "axis"
void setToRotate(float * const matrix, float angle, const float * const axis)
{
    float c=cos(angle);
    float s=sin(angle);
    float x=axis[0]; 
    float y=axis[1]; 
    float z=axis[2];

    if ((x==1.0) && (y==0.0) && (z==0.0))
    {
        float R[]={1.0, 0.0, 0.0, 0.0, 
                   0.0, c,   s,   0.0, 
                   0.0, -s,  c,   0.0, 
                   0.0, 0.0, 0.0, 1.0};
        for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
            matrix[iMatrixCoord]=R[iMatrixCoord];
    }
    else
    {
        if ((x==0.0) && (y==1.0) && (z==0.0))
        {                    
            float R[]={c,   0.0, -s,  0.0, 
                       0.0, 1.0, 0.0, 0.0, 
                       s,   0.0, c,   0.0, 
                       0.0, 0.0, 0.0, 1.0};
            for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
                matrix[iMatrixCoord]=R[iMatrixCoord];
        }
        else
        {

            if ((x==0.0) && (y==0.0) && (z==1.0))
            {                                          
                float R[]={c,   s,   0.0, 0.0, 
                           -s,  c,   0.0, 0.0, 
                           0.0, 0.0, 1.0, 0.0, 
                           0.0, 0.0, 0.0, 1.0};
                for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
                    matrix[iMatrixCoord]=R[iMatrixCoord];
            }
            else
            {
                float R[]={ (1.0-c)*(x*x-1.0) + 1.0, (1.0-c)*x*y + (z*s),     (1.0-c)*x*z - (y*s),      0.0, 
                            (1.0-c)*x*y - (z*s),     (1.0-c)*(y*y-1.0) + 1.0, (1.0-c)*y*z + (x*s),      0.0, 
                            (1.0-c)*x*z + (y*s),     (1.0-c)*y*z - (x*s),     (1.0-c)*(z*z-1.0) + 1.0,  0.0, 
                            0.0,                     0.0,                     0.0,                      1.0};
                for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
                    matrix[iMatrixCoord]=R[iMatrixCoord];
                //std::cout<<"Rotation on non standard axis."<<std::endl;
            }
        }
    }
}


// Builds a perspective projection matrix and stores it in mat
// l=left, r=right, b=bottom, t=top, n=near, f=far in the frustum
void setPerspective(float * const mat, float l, float r, float b, float t, float n, float f)
{
    float P[]={(2*n)/(r-l), 0.0, 0.0, 0.0,
               0.0, (2*n)/(t-b), 0.0, 0.0,
               (r+l)/(r-l), (t+b)/(t-b), -(f+n)/(f-n), -1.0,
               0.0, 0.0, 0.0, 0.0};
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        mat[iMatrixCoord]=P[iMatrixCoord];
}


// Does the multiplication A=A*B : all the matrices are described column-major
void multMatrixBtoMatrixA(float * const A, const float * const B)
{
    unsigned int i=0; // row index
    unsigned int j=0; // column index
    float temp[16];
    
    for (unsigned int iValue=0 ; iValue<16 ; ++iValue)
    {
        temp[iValue]=0;
        //j=iValue%4; // if raw-major
        //i=iValue/4; // if raw-major
        i=iValue%4; // if column-major
        j=iValue/4; // if column-major
        for (unsigned int k=0 ; k<4 ; ++k)
        {
            unsigned int indexik=k*4+i;
            unsigned int indexkj=j*4+k;
            temp[iValue]+=A[indexik]*B[indexkj];
        }
    }
    
    for (unsigned int iValue=0 ; iValue<16 ; ++iValue)
        A[iValue]=temp[iValue];
}


//______________________________________________________________________________
// Mathematics : interpolation


// Fills the linear interpolation of a nbCoord vector
void linearInterpolation(unsigned int nbCoords, const float * const x0, const float * const x1, float toX0, float * const result)
{
    for (unsigned int iCoord=0 ; iCoord<nbCoords ; ++iCoord)
        result[iCoord]=(1.0-toX0)*x0[iCoord]+toX0*x1[iCoord];
}


// Fills the bilinear interpolation of a nbCoord vector
void biLinearInterpolation(unsigned int nbCoords, const float * const x0y0, const float * const x1y0, const float * const x0y1, const float * const x1y1, float toX0, float toY0, float * const result)
{
    float x0[nbCoords]; linearInterpolation(nbCoords, x0y0, x0y1, toY0, x0);
    float x1[nbCoords]; linearInterpolation(nbCoords, x1y0, x1y1, toY0, x1);
    linearInterpolation(nbCoords, x0, x1, toX0, result); 
}


// Fills the trilinear interpolation of a nbCoord vector
void triLinearInterpolation(unsigned int nbCoords, const float * const x0y0z0, const float * const x1y0z0, const float * const x0y1z0, const float * const x1y1z0, const float * const x0y0z1, const float * const x1y0z1, const float * const x0y1z1, const float * const x1y1z1, float toX0, float toY0, float toZ0, float * const result)
{
    float x0[nbCoords]; biLinearInterpolation(nbCoords, x0y0z0, x0y1z0, x0y0z1, x0y1z1, toY0, toZ0, x0);
    float x1[nbCoords]; biLinearInterpolation(nbCoords, x1y0z0, x1y1z0, x1y0z1, x1y1z1, toY0, toZ0, x1);
    linearInterpolation(nbCoords, x0, x1, toX0, result);
}


//______________________________________________________________________________
// Random generation


// Generate a random double between min and max 
double getRand(MTRand & mtrand, bool uniform=true, double min=0.0, double max=1.0)
{
    double value=min;
    if (uniform) // uniform law
    {
        value=min+mtrand.rand(max-min);
    } 
    else // normal law
    {
        double variance=(max-min)/2.0;
        double standardDeviation=sqrt(variance);
        double average=min+variance;
        value=min+mtrand.randNorm(average, standardDeviation); 
    }
    return value;
}


//______________________________________________________________________________
// Printing


// Prints a vector of 2 values (textures 2D coordinates for example)
void printVec2(const float * const vect)
{
	std::cout<<"["<<vect[0]<<", "<<vect[1]<<"]"<<std::endl;
}


// Prints a vector of 3 values (a normal for example)
void printVec3(const float * const vect)
{
	std::cout<<"["<<vect[0]<<", "<<vect[1]<<", "<<vect[2]<<"]"<<std::endl;
}


// Prints a vector of 4 values (a position for example)
void printVec4(const float * const vect)
{
	std::cout<<"["<<vect[0]<<", "<<vect[1]<<", "<<vect[2]<<", "<<vect[3]<<"]"<<std::endl;
}


// Prints a matrix (written in math form, not in column-major 16 values array)
void printMat16(const float * const mat)
{
    std::cout<<std::endl;
    for (unsigned int iRow=0 ; iRow<4 ; ++iRow)
    {
        std::cout<<"| ";
        for (unsigned int iColumn=0 ; iColumn<4 ; ++iColumn)
        {
            std::cout<<mat[(iColumn*4)+iRow];
            std::cout<<" ";
        }
        std::cout<<"|"<<std::endl;
    }
}


// Prints the stack of glErrors, or nothing if no error occured
void printGLErrors()
{
    GLenum error = glGetError();
    
    // The glGetError function returns one error at a time, and then unstacks it.
    // Here we call it until all the errors are shown.
    while (error!=GL_NO_ERROR)
    {
        std::cout<<"!! GL Error : ";
        if (error==GL_INVALID_ENUM)
            std::cout<<"GL_INVALID_ENUM"<<std::endl;
        if (error==GL_INVALID_VALUE)
            std::cout<<"GL_INVALID_VALUE"<<std::endl;
        if (error==GL_INVALID_OPERATION)
            std::cout<<"GL_INVALID_OPERATION"<<std::endl;
        if (error==GL_OUT_OF_MEMORY)
            std::cout<<"GL_OUT_OF_MEMORY"<<std::endl;

        error=glGetError();

        std::cout<<std::endl;
    }
}


//______________________________________________________________________________
// Texture loading


// Loads a simple texture
unsigned int loadTexture(const std::string &fileName)
{
    unsigned int  w;
    unsigned int  h;
    // Loads the image from a ppm file to an unsigned char array
    unsigned char *data=loadPPM(fileName, &w, &h);

    // Allocates a texture id
    unsigned int textureID;
    glGenTextures(1, &textureID);
    // Selects our current texture
    glBindTexture(GL_TEXTURE_2D, textureID);

    // How to handle not normalised uvs
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // How to handle interpolation from texels to fragments
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Specifies which image will be used for this texture objet
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    
    delete [] data;

    return textureID;
}


//______________________________________________________________________________
// Images reading and writing


// Loads an image encoded in ppm format (P6)
unsigned char * loadPPM(const std::string & filename, unsigned int * const width, unsigned int * const height)
{	
    unsigned char * data;
    // Reads header
	std::ifstream file(filename.c_str());
	std::string line;
	
	// Formats
	std::getline(file, line);
	
	// Skips comments
	std::getline(file, line);
	while(line[0]=='#') std::getline(file, line);

	// Reads dimensions
	std::istringstream ist(line);
    unsigned int w, h;
	ist>>w>>h;
    *width=w;
    *height=h;

    // Reads the data
	std::getline(file, line);

    /* May be necessary on windows	   
    // 	unsigned int dataStart = file.tellg();	
    // 	file.close();
    // 	file.open(filename.c_str(),std::ios::in | std::ios::binary); 
    // 	file.seekg(dataStart);
	*/

	data=new unsigned char[w*h*3];
	file.read((char*)data, w*h*sizeof(unsigned char)*3);

	// Closes the file
	file.close();

    return data;
}


// Saves frameBuffer as a PPM
void saveFrameBufferPPM(const std::string & fileName, unsigned int width, unsigned int height)
{
	unsigned char * data=new unsigned char[width*height*4];
	unsigned char * dataRGB=new unsigned char[width*height*3];
	glReadPixels( 0,
				  0,
				  width,
				  height,
				  GL_RGBA,
				  GL_UNSIGNED_BYTE,
				  data);
	for (unsigned int i=0 ; i<height; ++i)
	{
		for (unsigned int j=0 ; j<width; ++j)
		{
			dataRGB[((height-i-1)*width+j)*3]=data[(i*width+j)*4];
			dataRGB[((height-i-1)*width+j)*3+1]=data[(i*width+j)*4+1];
			dataRGB[((height-i-1)*width+j)*3+2]=data[(i*width+j)*4+2];
		}
	}

	FILE * file=fopen(fileName.c_str(), "w");
	fprintf(file, "P6\n");
	fprintf(file, "%d %d\n", width, height);
	fprintf(file, "255\n");
	fwrite(dataRGB, width*height*3, sizeof(unsigned char), file);
	fclose(file);
	delete[] data;
}
