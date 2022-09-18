//
// Created by stein on 22/02/2021.
//

#include "L3.h"

L3::L3(string &L3, Figure* f) {
    std::ifstream fin(L3);
    fin >> lsystem;
    fin.close();

    F = f;

    angle = degree_to_rad(lsystem.get_angle());

    currentH = Vector3D::point(1,0,0);
    currentL = Vector3D::point(0,1,0);
    currentU = Vector3D::point(0,0,1);
    string generate = lsystem.get_initiator();
    f->points.push_back(currentPoint);
    generateList_pv(generate, lsystem.get_nr_iterations(), f);
}

void L3::proccessChar(char v, Figure* f) {
    double currentXPossition = currentPoint.x + currentH.x;
    double currentYPossition = currentPoint.y + currentH.y;
    double currentZPossition = currentPoint.z + currentH.z;

    Vector3D point1 = Vector3D::point(currentXPossition, currentYPossition, currentZPossition);

    f->points.push_back(point1);
    int point1Location = f->points.size() - 1;
    int point2Location = f->points.size() - 2;
    f->faces.push_back(vector<int>({point1Location, point2Location}));
    currentPoint = point1;
}

bool L3::update_currentPoint(char v) {
    switch (v) {
        case '(':
            Hstack.push(currentH);
            Lstack.push(currentL);
            Ustack.push(currentU);
            pointStack.push(currentPoint);
            return true;
        case ')':
            currentH = Hstack.top();
            Hstack.pop();
            currentL = Lstack.top();
            Lstack.pop();
            currentU = Ustack.top();
            Ustack.pop();

            currentPoint = pointStack.top();
            pointStack.pop();
            F->points.push_back(currentPoint);
            return true;
        default:
            return false;
    }
}

bool L3::update_angle(char v) {
    Vector3D currentLTemp = currentL;
    Vector3D currentHTemp = currentH;
    Vector3D currentUTemp = currentU;
    switch (v) {
        case '+':
            currentH = currentH * cos(angle) + currentL * sin(angle);
            currentL = -currentHTemp * sin(angle) + currentL * cos(angle);
            return true;
        case '-':
            currentH = currentH * cos(-angle) + currentL * sin(-angle);
            currentL = -currentHTemp * sin(-angle) + currentL * cos(-angle);
            return true;
        case '^':
            currentH = currentH * cos(angle) + currentU * sin(angle);
            currentU = -currentHTemp * sin(angle) + currentU * cos(angle);
            return true;
        case '&':
            currentH = currentH * cos(-angle) + currentU * sin(-angle);
            currentU = -currentHTemp * sin(-angle) + currentU * cos(-angle);
            return true;
        case '\\':
            currentL = currentL * cos(angle) - currentU * sin(angle);
            currentU = currentLTemp * sin(angle) + currentU * cos(angle);
            return true;
        case '/':
            currentL = currentL * cos(-angle) - currentU * sin(-angle);
            currentU = currentLTemp * sin(-angle) + currentU * cos(-angle);
            return true;
        case '|':
            currentH = -currentH;
            currentL = -currentHTemp;
            return true;
        default:
            return false;
    }
}
void L3::generateList_pv(string &Initiator, int itteration, Figure* f) {
    Lines2D lines;
    set<char> alfabet = lsystem.get_alphabet();
    string itString = lsystem.get_initiator();
    for(int i = 0; i < Initiator.size(); i++){
        if(itteration > 0){
            if(!update_angle(Initiator[i]) && !update_currentPoint(Initiator[i])) {
                string replacement = lsystem.get_replacement(Initiator[i]);
                generateList_pv(replacement,itteration - 1, f);
            }
        }
        else if(!update_angle(Initiator[i]) && !update_currentPoint(Initiator[i])) {
            proccessChar(Initiator[i],f);
        }
    }
}

double L3::degree_to_rad(double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}