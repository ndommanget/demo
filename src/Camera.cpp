// Camera.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Camera.hpp"
#include "Tools.hpp"



//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Default constructor
Camera::Camera(const Shaders & shaders) : shaders(shaders)
{

}


// Cleans memory for Camera
Camera::~Camera()
{

}


//______________________________________________________________________________
// Initialization


// Inits parameters for the camera
void Camera::init()
{
    // View Data
    // Camera position and orientation  
    float c[]={0.0, 1.5, 2.2}; // Camera position    
    float aim[]={0.0, 1.5, 0.0}; // Where we look
    float up[]={0.0, 1.0, 0.0}; // Vector pointing over the camera
    lookAt(c, aim, up);

    // Projection data
    // Projection type : perspective:true / orthographic:false
    this->perspectiveProjection=true;
    // Projection frustum data
    float l=1.0;
    this->left=-l; 
    this->right=l;
    this->bottom=-l;
    this->top=l; 
    this->near=0.1;
    this->far=100.0;
    setToIdentity(this->projection);
    updateProjection();
}


//______________________________________________________________________________
// View set-up


// Builds the camera axis from position c, aim (focus) of the camera, and up vector
void Camera::lookAt(const float * const c, const float * const aim, const float * const up)
{
    for (unsigned int iCoord=0 ; iCoord<3; ++iCoord)
    {
        this->c[iCoord]=c[iCoord]; // c : camera position
        this->y[iCoord]=up[iCoord]; // y : up vector
        this->z[iCoord]=c[iCoord]-aim[iCoord]; // z : from aim to camera position
    }
    // Opposite direction of axis z
    float minusZ[]={-this->z[0], -this->z[1], -this->z[2]};
    // axis x  : from axis y and axis z
    vectorProduct(this->y, this->z, this->x);
    // axis y  : from new axis x and opposite of axis z
    vectorProduct(this->x, minusZ, this->y);
    // Normalizes all axis vectors
    normalize(this->x);
    normalize(this->y);
    normalize(this->z);

    // Builds the new view matrix
    updateView();
}


// Moves
void Camera::move(const bool spherical, const float * const moves, const float * const angles, const float radius)
{
    float cameraNewPos[3];
    for (unsigned int iCoord=0 ; iCoord<3 ; ++iCoord)
    {
        cameraNewPos[iCoord]=this->c[iCoord]
                            +this->x[iCoord]*moves[0]
                            +this->y[iCoord]*moves[1]
                            +this->z[iCoord]*moves[2];
    }

    float xAxis[]={1.0, 0.0, 0.0}; float yAxis[]={0.0, 1.0, 0.0};
    float rotateAroundX[16]; setToRotate(rotateAroundX, -angles[1], xAxis);
    float rotateAroundY[16]; setToRotate(rotateAroundY, angles[0], yAxis);
    float t[]={-cameraNewPos[0], -cameraNewPos[1], -cameraNewPos[2]};
    float translate[16]; setToTranslate(translate, t);

    setToIdentity(this->view);

    if (spherical)
    {
        float tRadius[]={0.0, 0.0, -radius};
        float translateRadius[16]; setToTranslate(translateRadius, tRadius);
        multMatrixBtoMatrixA(this->view, translateRadius);
    }

    multMatrixBtoMatrixA(this->view, rotateAroundX);
    multMatrixBtoMatrixA(this->view, rotateAroundY);
    multMatrixBtoMatrixA(this->view, translate);

    for (unsigned int iCoord=0 ; iCoord<3 ; ++iCoord)
    {
            // Updates the axis with values in view
            this->x[iCoord]=this->view[iCoord*4+0];
            this->y[iCoord]=this->view[iCoord*4+1];
            this->z[iCoord]=this->view[iCoord*4+2];
            // Updates the position of the camera c
            this->c[iCoord]=cameraNewPos[iCoord];
    }

    // Updates in shaders
    this->shaders.setView(this->view, this->c);
}


//______________________________________________________________________________
// Projection set-up


// Sores the data necassary to evaluate the porjection matrix
void Camera::setProjectionData(float left, float right, float bottom, float top, float near, float far)
{
    this->left=left;
    this->right=right;
    this->bottom=bottom;
    this->top=top;
    this->near=near;
    this->far=far;
}


// Turns camera projection from perspective to ortho or inverse
void Camera::switchCameraProjection()
{
    std::cout<<"Camera projection : ";
    if (this->perspectiveProjection)
    {
        std::cout<<"orthographic"<<std::endl;
        this->perspectiveProjection=false;
    }
    else
    {
        std::cout<<"perspective"<<std::endl;
        this->perspectiveProjection=true;
    }

    // Changes the matrix accordingly
    updateProjection();
}


// Sets the perspective similarly to gluPerspective
void Camera::setPerspectiveFromAngle(float fovy, float aspectRatio)
{
    this->top=this->near*tan(fovy/2.0);
    this->bottom=-this->top;
    this->left=this->bottom*aspectRatio;
    this->right=this->top*aspectRatio;
    
    updateProjection();
}


//______________________________________________________________________________
// Internal updates


// Updates view
void Camera::updateView()
{
    // Rotation to be aligned with correct camera axis
    float RcInv[]={this->x[0], this->y[0], this->z[0], 0.0,
                   this->x[1], this->y[1], this->z[1], 0.0,
                   this->x[2], this->y[2], this->z[2], 0.0,
                   0.0,        0.0,        0.0,        1.0};

    // Translation to be at the right distance from the scene
    float TcInv[]={1.0,         0.0,         0.0,         0.0,
                   0.0,         1.0,         0.0,         0.0,
                   0.0,         0.0,         1.0,         0.0,
                   -this->c[0], -this->c[1], -this->c[2], 1.0};

    // Initializes
    setToIdentity(this->view);
    // Rotates
    multMatrixBtoMatrixA(this->view, RcInv);  
    // Translates
    multMatrixBtoMatrixA(this->view, TcInv);

    // Updates in shaders
    this->shaders.setView(this->view, this->c);
}


// Updates the projection matrix from the data
void Camera::updateProjection()
{  
    float l=this->left;
    float r=this->right;
    float b=this->bottom;
    float t=this->top;
    float n=this->near;
    float f=this->far;

    if (this->perspectiveProjection) // Perspective projection
    {
        float P[]={(2.0*n)/(r-l), 0.0,           0.0,              0.0,
                   0.0,           (2.0*n)/(t-b), 0.0,              0.0,
                   (r+l)/(r-l),   (t+b)/(t-b),   -(f+n)/(f-n),    -1.0,
                   0.0,           0.0,           -(2.0*f*n)/(f-n), 0.0};

        for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
            this->projection[iMatrixCoord]=P[iMatrixCoord];
    }
    else // Orthographic projection
    { 
        float P[]={ 2.0/(r-l),   0.0,          0.0,          0.0,
                    0.0,         2.0/(t-b),    0.0,          0.0,
                    0.0,         0.0,          -2.0/(f-n),   0.0,
                   -(r+l)/(r-l), -(t+b)/(t-b), -(f+n)/(f-n), 1.0};

        for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
            this->projection[iMatrixCoord]=P[iMatrixCoord];
    }

    // Updates in shaders
    this->shaders.setProjection(this->projection);
}
