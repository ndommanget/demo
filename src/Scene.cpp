// Scene.cpp
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com


#include "Scene.hpp"
#include "Tools.hpp"
#include "Object.hpp"



//______________________________________________________________________________
// Drawing preparation


// Decides what will elements drawn after this call will look like
void Scene::setAppearance(unsigned int indexDrawnObject) const
{
    // If sprite shader, transparency is expected
    if (this->shaders.isTransparent(drawnObjectsShaders[indexDrawnObject]))
    {
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
    }
    else
    {
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_BLEND);
    }

    // We use the specific values of model per object
    this->shaders.setModel(this->drawnObjectsShaders[indexDrawnObject], &(this->drawnObjectsModels[indexDrawnObject*16]));

    // We use the specific values of color per object
    this->shaders.changeMaterialColor(this->drawnObjectsShaders[indexDrawnObject], &(this->drawnObjectsColors[indexDrawnObject*4]));

    // Specifies which VBO were filled
    const Object * object=this->storedObjects[this->drawnObjects[indexDrawnObject]];
    this->shaders.setFilledData(this->drawnObjectsShaders[indexDrawnObject], object->hasPrimitives(), object->hasNormals(), object->hasUvs(), object->hasColors());

    // Selects our current texture for unit 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, this->drawnObjectsTexture0IDs[indexDrawnObject]);

    // Selects our current texture for unit 1
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, this->drawnObjectsTexture1IDs[indexDrawnObject]);
}


//______________________________________________________________________________
//______________________________________________________________________________
// Constructors / destructor


// Default constructor
Scene::Scene(Shaders & shaders, Camera & camera) : shaders(shaders), camera(camera)
{

}


// Cleans memory for Scene
Scene::~Scene()
{
    for (GLuint iStoredObjects=0 ; iStoredObjects<this->nbStoredObjects ; iStoredObjects++)
        delete this->storedObjects[iStoredObjects];

    delete [] this->storedObjects;
    delete [] this->drawnObjects;
    delete [] this->drawnObjectsColors;
    delete [] this->drawnObjectsModels;
    delete [] this->drawnObjectsShaders;
    delete [] this->drawnObjectsTexture0IDs;
    delete [] this->drawnObjectsTexture1IDs;
}


//______________________________________________________________________________
// Initialization


// Initializes all the pointers to NULL
// and sets the numbers of objects to 0
void Scene::init()
{
    // Fixed max sizes for arrays
    this->maxStoredObjects=50;
    this->maxDrawnObjects=200;
    
    // Initialisation of the numbers of objets filled
    this->nbStoredObjects=0;
    this->nbDrawnObjects=0;
    
    // Array creation (with fixed size)
    this->storedObjects=new const Object * [this->maxStoredObjects];
    this->drawnObjects=new unsigned int[this->maxDrawnObjects];
    this->drawnObjectsColors=new float[this->maxDrawnObjects*4];
    this->drawnObjectsModels=new float[this->maxDrawnObjects*16];
    this->drawnObjectsShaders=new unsigned int[this->maxDrawnObjects];
    this->drawnObjectsTexture0IDs=new unsigned int[this->maxDrawnObjects];
    this->drawnObjectsTexture1IDs=new unsigned int[this->maxDrawnObjects];
           
    // Default color
    float grey[16]={0.2, 0.2, 0.2, 1.0}; setDefaultColor(grey);
    // Default model
    float I[16]; setToIdentity(I); setDefaultModel(I);
    // Default shader ID
    setDefaultShader(1);
    // Default texture ID   
    setDefaultTextureID(1);
    
     // Initalisation of storedObjects
    for (unsigned int iStoredObjects=0 ; iStoredObjects<this->maxStoredObjects ; ++iStoredObjects)
        this->storedObjects[iStoredObjects]=NULL;
        
    // Light creation
    float lightPosition[]={0.0, 5.0, 0.0, 1.0}; float lightPower=1.0;
    this->setLight(lightPosition, lightPower);

    // Default shader
    this->setDefaultShader(0);
}


//______________________________________________________________________________
// Building


// Adds the object in the Objects library array, 
// after the last added object, and only if the array is not full, returns the index
unsigned int Scene::storeObject(const Object * const object)
{
	if (this->nbStoredObjects<this->maxStoredObjects)
	{
		this->storedObjects[this->nbStoredObjects]=object;
		//std::cout<<"Object stored (index: "<<this->nbStoredObjects<<")."<<std::endl;
		
		this->nbStoredObjects++;
	}
	else
	{
		std::cout<<"There is already "<<this->maxStoredObjects<<" objects stored."<<std::endl;
	}
	return (this->nbStoredObjects-1);
}


// Adds the index of stored object to draw in the drawnObjects array, 
// after the last added, and only if the array is not full
unsigned int  Scene::addObjectToDraw(unsigned int indexStoredObject)
{
    if (this->nbDrawnObjects<this->maxDrawnObjects)
	{
		this->drawnObjects[this->nbDrawnObjects]=indexStoredObject;
        setDrawnObjectColor(this->nbDrawnObjects, this->defaultColor);
        setDrawnObjectModel(this->nbDrawnObjects, this->defaultModel);
        setDrawnObjectShader(this->nbDrawnObjects, this->defaultShader);
        setDrawnObjectTextureID(this->nbDrawnObjects, 0, this->defaultTextureID);
        setDrawnObjectTextureID(this->nbDrawnObjects, 1, this->defaultTextureID);
        
		//std::cout<<"Object "<<this->nbDrawnObjects<<" to draw."<<std::endl;
		this->nbDrawnObjects++;
	}
	else
	{
		std::cout<<"There is already "<<this->maxDrawnObjects<<" objects drawn."<<std::endl;
	}
    return (this->nbDrawnObjects-1);
}


//______________________________________________________________________________
// Default set-up


// Sets default color
void Scene::setDefaultColor(const float * const defaultColor)
{
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        this->defaultColor[iCoord]=defaultColor[iCoord];

    for (unsigned int iDrawnObjects=0 ; iDrawnObjects<this->maxStoredObjects ; ++iDrawnObjects)
        setDrawnObjectColor(iDrawnObjects, this->defaultColor);
}


// Sets default transformation matrix
void Scene::setDefaultModel(const float * const defaultModel)
{
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        this->defaultModel[iMatrixCoord]=defaultModel[iMatrixCoord];

    for (unsigned int iDrawnObjects=0 ; iDrawnObjects<this->maxStoredObjects ; ++iDrawnObjects)
        setDrawnObjectModel(iDrawnObjects, this->defaultModel);
}


// Sets default shader
void Scene::setDefaultShader(unsigned int defaultShader)
{
    this->defaultShader=defaultShader;
}


// Sets default texture ID
void Scene::setDefaultTextureID(unsigned int defaultTextureID)
{
    this->defaultTextureID=defaultTextureID;
}


//______________________________________________________________________________
// Drawn objects set-up


// Sets the color for the drawn object of index indexDrawObject
void Scene::setDrawnObjectColor(unsigned int indexDrawnObject, const float * const color)
{
    for (unsigned int iColorComp=0 ; iColorComp<4 ; ++iColorComp)
        this->drawnObjectsColors[(indexDrawnObject*4)+iColorComp]=color[iColorComp];
}


// Sets the transformation matrix for the drawn object of index indexDrawObject
void Scene::setDrawnObjectModel(unsigned int indexDrawnObject, const float * const model)
{
    for (unsigned int iMatrixCoord=0 ; iMatrixCoord<16 ; ++iMatrixCoord)
        this->drawnObjectsModels[(indexDrawnObject*16)+iMatrixCoord]=model[iMatrixCoord];

}


// Sets the shader to use on the drawn object of index indexDrawObject
void Scene::setDrawnObjectShader(unsigned int indexDrawnObject, unsigned int shader)
{
    this->drawnObjectsShaders[indexDrawnObject]=shader;
}


// Sets the ID of the texture to use on the drawn object of index indexDrawObject
void Scene::setDrawnObjectTextureID(unsigned int indexDrawnObject, unsigned int textureUnit, unsigned int textureID) 
{
    if (textureUnit==0) this->drawnObjectsTexture0IDs[indexDrawnObject]=textureID;
    if (textureUnit==1) this->drawnObjectsTexture1IDs[indexDrawnObject]=textureID;
}
 

//______________________________________________________________________________
// Light set-up


// Sets light data to use in shader
void Scene::setLight(const float * const position, float power)
{
    for (unsigned int iCoord=0 ; iCoord<4 ; ++iCoord)
        this->lightPosition[iCoord]=position[iCoord];
    this->lightPower=power;

    this->shaders.setLight(this->lightPosition, this->lightPower);
}


//______________________________________________________________________________
//Drawing


// Draw all Objects
void Scene::draw() const
{
    
    unsigned int indexStoredObjects=0;
    for (unsigned int iDrawnObjects=0 ; iDrawnObjects<this->nbDrawnObjects; ++iDrawnObjects)
    {
        indexStoredObjects=this->drawnObjects[iDrawnObjects];
	    if ((indexStoredObjects<this->maxStoredObjects) && (this->storedObjects[indexStoredObjects]!=NULL))
	    {
	        setAppearance(iDrawnObjects);    
    	    this->storedObjects[indexStoredObjects]->drawObject();
        }
    }
}
