//
// Created by stein on 22/02/2021.
//

#ifndef ENGINE_L2_H
#define ENGINE_L2_H
#include "easy_image.h"
#include "ini_configuration.h"
#include <string>
#include <fstream>
#include <iostream>
#include "l_parser.h"
#include <math.h>
#include <stack>


using namespace img;

class L2D {
protected:
    LParser::LSystem2D lsystem;
    double curren_angle = 0;
    Point2D currentPoint;
    Color c;
    stack<Point2D>pointStack;
    stack<double>angleStack;
    double degree_to_rad(double degree);

    bool update_angle(char v);
    bool update_currentPoint(char v);

    Line2D proccessChar(char v);
    Lines2D generateList_pv(string &Initiator, int itteration);

    double maxX = 0;
    double maxY = 0;

public:
    L2D(string &l2d, Color &c);
    Lines2D generateList();
    double x_y_verhouding();
};

#endif // EvNGINE_L2_H
