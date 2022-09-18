//
// Created by stein on 18/03/2021.
//

#ifndef ENGINE_L3_H
#define ENGINE_L3_H
#include "easy_image.h"
#include "ini_configuration.h"
#include <string>
#include <fstream>
#include <iostream>
#include "l_parser.h"
#include <math.h>
#include <stack>

using namespace img;

class L3 {
protected:
    LParser::LSystem3D lsystem;
    Vector3D currentH;
    Vector3D currentL;
    Vector3D currentU;
    Vector3D currentPoint;

    Figure* F;

    stack<Vector3D>Hstack;
    stack<Vector3D>Lstack;
    stack<Vector3D>Ustack;
    stack<Vector3D>pointStack;
    double degree_to_rad(double degree);

    double angle;

    bool update_angle(char v);
    bool update_currentPoint(char v);

    void proccessChar(char v, Figure* f);
    void generateList_pv(string &Initiator, int itteration, Figure* f);
public:
    L3(string &l2d, Figure* f);
    double x_y_verhouding();
};
#endif //ENGINE_L3_H
