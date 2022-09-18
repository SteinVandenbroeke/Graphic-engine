//
// Created by stein on 7/05/2021.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "easy_image.h"
#include "ini_configuration.h"

class Light {
public:
    Light(img::Color ambientLight, img::Color diffuseLight, img::Color specularLight){
        this->ambientLight = ambientLight;
        this->diffuseLight = diffuseLight;
        this->specularLight = specularLight;
    }
    //de ambiente licht component
    img::Color ambientLight;
    //de diffuse licht component
    img::Color diffuseLight;
    //de diffuse licht component
    img::Color specularLight;
    virtual Vector3D getlVector(){return Vector3D();}
    virtual Vector3D getLocation(){cout<< "fout" << endl; return Vector3D();}
    virtual bool isInfLight() {return false;}
    virtual bool isPointLight() {return false;}
    virtual double getspotAngle() {return 360;}
};

class InfLight : public Light {
public:
    InfLight(img::Color ambientLight, img::Color diffuseLight, img::Color specularLight, Vector3D direction): Light(ambientLight, diffuseLight, specularLight){
        this->lVector = -Vector3D::normalise(direction);
    }
    //de richting waarin het
    //licht schijnt
    Vector3D getlVector() override { return lVector;}
    bool isInfLight() override {return true;}
private:
    Vector3D lVector;
};

class PointLight : public Light {
private:
    //de locatie van de puntbron
    Vector3D location1;
    double spotAngle;
public:
    PointLight(img::Color ambientLight, img::Color diffuseLight, img::Color specularLight, Vector3D location, double spotAngle): Light(ambientLight, diffuseLight, specularLight){
        this->location1 = location;
        this->spotAngle = spotAngle;
    }
    //de hoek van een spotlicht
    Vector3D getLocation() override{
        return this->location1;
    }

    double getspotAngle() override {return this->spotAngle;}
    bool isPointLight() override {return true;}
};

typedef std::list<Light*> Lights3D;


#endif //ENGINE_LIGHT_H
