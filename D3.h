//
// Created by stein on 5/03/2021.
//

#ifndef ENGINE_D3_H
#define ENGINE_D3_H
#include "vector3d.h"
#include "easy_image.h"
#include <math.h>
#include "ini_configuration.h"
#include "L3.h"

using namespace img;
class D3 {
private:
    double pi = 3.141592653;
    Figures3D * figures;
    Matrix* all;
    Matrix makeSchaalMatrix(const double scaleFactor);
    Matrix makeX_asMatrix(const double angel);
    Matrix makeY_asMatrix(const double angel);
    Matrix makeZ_asMatrix(const double angel);
    Matrix makeTranslatieMatrix(double x, double y, double z);
    double degree_to_rad(double degree);

    Matrix complete_matrix_1(const double scaleFactor, const double angel_x, const double angel_y, const double angel_z,
                      const double transmatie_x, const double transmatie_y, const double transmatie_z);
    Matrix* complete_matrix(const double scaleFactor, const double angel_x, const double angel_y, const double angel_z, const double transmatie_x, const double transmatie_y,const double transmatie_z);
    Vector3D* Transformaties();
    void applyTransformation(Figure*, const Matrix&);
    Matrix eyePointTrans(const Vector3D &eyepoint);
    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    Point2D doProjection(const Vector3D &point, const double d);
public:
    Matrix eyePointTransMatrix;
    Vector3D eyeCordinaat;
    D3();
    ~D3();
    Lines2D doProjection(const Figures3D* figuren3D);
    Figure* createFigure(const ini::Configuration &configuration, string& name);
    Figures3D* getAllFigures(const ini::Configuration &configuration);
    Figure* createKubus(Figure* fg);
    Figure* createTetrahedron(Figure* fg);
    Figure* createOctahedron(Figure* fg);
    Figure* createIcosahedron(Figure* fg);
    Figure* createDodecahedron(Figure* fg);
    Figure* createbuckyball(Figure* fg);
    Figure* createSphere(Figure *fg, const double radius, const int n);
    Figure* createCone(Figure *fg, const int n, const double h);
    Figure * createCylinder(Figure *fg, const int n, const double h);
    Figure * createTorus(Figure *fg, const double r, const double R, const int n, const int m);
    Vector3D rescalePointSphere(Vector3D &point);
    void TriangulateAll(const Figures3D* figuren3D);
    void generateFractal(Figure* fig, const int nr_iterations, const double scale, const double, const double, const double, const double, const double,const double, Matrix& eyePoint);
    void generateMergeSpons(Figure* fig,const int nr_iterations,const double rotateX,const double rotateY, const double rotateZ, const double translateX,const double translateY,const double translateZ, Matrix& eyePoint);
    list<Face> triangulate(const Face& face);
};


#endif //ENGINE_L3H
