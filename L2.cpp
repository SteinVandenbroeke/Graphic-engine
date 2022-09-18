//
// Created by stein on 22/02/2021.
//

#include "L2.h"

L2D::L2D(string &l2d, Color &c) {
    std::ifstream fin(l2d);
    fin >> lsystem;
    fin.close();

    curren_angle = lsystem.get_starting_angle();
    currentPoint.set(0,0);
    this->c = c;
}

Line2D L2D::proccessChar(char v) {
    double currentXPossition = currentPoint.getX();
    double currentYPossition = currentPoint.getY();
    double draw = (double)lsystem.draw(v);
    double cosV = (double)cos(degree_to_rad(curren_angle));
    double sinV = (double)sin(degree_to_rad(curren_angle));

    double x = currentXPossition + draw * cosV;
    double y = currentYPossition + draw * sinV;
    Point2D point1(x,y);
    Line2D line(currentPoint, point1, c);
    currentPoint = point1;
    if(x > maxX) maxX = x;
    if(y > maxY) maxY = y;

    return line;
}

bool L2D::update_currentPoint(char v) {
    switch (v) {
        case '(':
            pointStack.push(currentPoint);
            angleStack.push(curren_angle);
            return true;
        case ')':
            currentPoint = pointStack.top();
            pointStack.pop();
            curren_angle = angleStack.top();
            angleStack.pop();
            return true;
        default:
            return false;
    }
}

bool L2D::update_angle(char v) {
    switch (v) {
        case '-':
            curren_angle = curren_angle - lsystem.get_angle();
            return true;
        case '+':
            curren_angle = curren_angle + lsystem.get_angle();
            return true;
        default:
            return false;
    }
}

Lines2D L2D::generateList() {
    string generate = lsystem.get_initiator();
    return generateList_pv(generate, lsystem.get_nr_iterations());
}

Lines2D L2D::generateList_pv(string &Initiator, int itteration) {
    Lines2D lines;
    set<char> alfabet = lsystem.get_alphabet();
    string itString = lsystem.get_initiator();
    for(int i = 0; i < Initiator.size(); i++){
        if(itteration > 0){
            if(!update_angle(Initiator[i]) && !update_currentPoint(Initiator[i])) {
                string replacement = lsystem.get_replacement(Initiator[i]);
                Lines2D lines1 = generateList_pv(replacement,itteration - 1);
                copy(lines1.rbegin(), lines1.rend(), front_inserter(lines));
            }
        }
        else if(!update_angle(Initiator[i]) && !update_currentPoint(Initiator[i])) {
            lines.push_back(proccessChar(Initiator[i]));
        }
    }
    return lines;
}

double L2D::x_y_verhouding() {
    cout << "x: " << maxX << endl;
    cout << "y: " << maxY << endl;
    double deling = maxX/maxY;
    cout << deling <<endl;
    return deling / 1000;
}
double L2D::degree_to_rad(double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}