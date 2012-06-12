// Camera.hpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#ifndef __CAMERA_HPP__
#define __CAMERA_HPP__


#include "Shaders.hpp"



// A camera to watch the scene
class Camera
{
    private:

        // View Data
        float c[3];                 // Camera position 
        float x[3];                 // Camera axis x : right side
        float y[3];                 // Camera axis y : up
        float z[3];                 // Camera axis z : backward
        float view[16];             // View matrix

        // Projection data
        bool perspectiveProjection; // Persepective projection:true 
                                    // Orthographic projection:false
        
        float left;                 // x coord from center to left plane of frustum
        float right;                // x coord from center to right plane of frustum
        float bottom;               // y coord from center to bottom plane of frustum
        float top;                  // y coord from center to top plane of frustum
        float near;                 // z coord from c to near plane of frustum
        float far;                  // z coord from c to far plane of frustum
        float projection[16];       // Projection matrix

        // Users
        const Shaders & shaders;    // Shaders using the camera



    public:

        // Constructors / destructor
        Camera(const Shaders & shaders);
        ~Camera();

        // Initialization
        void init();

        // View set-up
        void lookAt(const float * const c, const float * const aim, const float * const up);
        void move(const bool spherical, const float * const moves, const float * const angles, float radius);

        // Data using
        void inView(float * const position) const;
        const float * const getCenter() const;

        // Projection set-up
        void setProjectionData(float left, float right, float bottom, float top, float near, float far);
        void switchCameraProjection();
        void setPerspectiveFromAngle(float fovy, float aspectRatio);

        // Internal updates
        void updateView();
        void updateProjection();
};


#endif // __CAMERA_HPP__
