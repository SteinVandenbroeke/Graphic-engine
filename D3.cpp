//
// Created by stein on 5/03/2021.
//

#include "D3.h"
#include <math.h>

Vector3D *D3::Transformaties() {
    pi = 2*acos(0.0);
    Matrix schalen;
    Matrix x_as;
    Matrix y_as;
    Matrix z_as;
    return NULL;
}

Matrix D3::makeSchaalMatrix(const double scaleFactor) {
    Matrix schalen;
    schalen(1,1) = scaleFactor;
    schalen(2,2) = scaleFactor;
    schalen(3,3) = scaleFactor;
    schalen(4,4) = 1;
    return schalen;
}

Matrix D3::makeX_asMatrix(const double angel) {
    Matrix x_as;
    x_as(2,2) = (double)cos(degree_to_rad(angel));
    x_as(2,3) = (double)sin(degree_to_rad(angel));
    x_as(3,2) = -(double)sin(degree_to_rad(angel));
    x_as(3,3) = (double)cos(degree_to_rad(angel));
    return x_as;
}

Matrix D3::makeY_asMatrix(const double angel) {
    Matrix y_as;
    y_as(1,1) = (double)cos(degree_to_rad(angel));
    y_as(3,1) = (double)sin(degree_to_rad(angel));
    y_as(1,3) = -(double)sin(degree_to_rad(angel));
    y_as(3,3) = (double)cos(degree_to_rad(angel));
    return y_as;
}

Matrix D3::makeZ_asMatrix(const double angel) {
    Matrix z_as;
    z_as(1,1) = (double)cos(degree_to_rad(angel));
    z_as(2,1) = -(double)sin(degree_to_rad(angel));
    z_as(1,2) = (double)sin(degree_to_rad(angel));
    z_as(2,2) = (double)cos(degree_to_rad(angel));
    return z_as;
}

Matrix D3::makeTranslatieMatrix(double x, double y, double z) {
    Matrix translatie;
    translatie(4,1) = x;
    translatie(4,2) = y;
    translatie(4,3) = z;
    return translatie;
}

Matrix* D3::complete_matrix(const double scaleFactor, const double angel_x, const double angel_y, const double angel_z,
                            const double transmatie_x, const double transmatie_y, const double transmatie_z) {
    all = new Matrix();
    *all = (makeSchaalMatrix(scaleFactor)) * (makeX_asMatrix(angel_x) * makeY_asMatrix(angel_y) * makeZ_asMatrix(angel_z)) * (makeTranslatieMatrix(transmatie_x,transmatie_y,transmatie_z));
    return all;
}

void D3::applyTransformation(Figure* f , const Matrix& m) {
    for(int i = 0; i < f->points.size(); i++){
        f->points[i] *= m;
    }
}

Figure *D3::createFigure(const ini::Configuration &configuration, string& name) {
    int pointCount = configuration[name]["nrPoints"].as_int_or_default(0);
    int lineCount = configuration[name]["nrLines"].as_int_or_default(0);
    string type = configuration[name]["type"].as_string_or_default("");
    double rotateX = configuration[name]["rotateX"].as_double_or_die();
    double rotateY = configuration[name]["rotateY"].as_double_or_die();
    double rotateZ = configuration[name]["rotateZ"].as_double_or_die();
    double scale = configuration[name]["scale"].as_double_or_die();
    vector<double> translate = configuration[name]["center"].as_double_tuple_or_die();
    vector<double> color = configuration[name]["color"].as_double_tuple_or_default(vector<double>({1,1,1}));
    vector<double> ambientReflection = configuration[name]["ambientReflection"].as_double_tuple_or_default(vector<double>({color[0],color[1],color[2]}));
    Color ambR(ambientReflection[0], ambientReflection[1], ambientReflection[2]);

    vector<double> diffuseReflection = configuration[name]["diffuseReflection"].as_double_tuple_or_default(vector<double>({1,1,1}));
    Color diffR(ambientReflection[0], ambientReflection[1], ambientReflection[2]);

    vector<double> specularReflection = configuration[name]["specularReflection"].as_double_tuple_or_default(vector<double>({0,0,0}));
    Color specR(specularReflection[0], specularReflection[1], specularReflection[2]);

    double reflectionCoefficient = configuration[name]["reflectionCoefficient"].as_double_or_default(0);



    double translateX = translate[0];
    double translateY = translate[1];
    double translateZ = translate[2];

    Figure* fg = new Figure();

    if(type.find("Cube") != std::string::npos || type == "MengerSponge"){
        createKubus(fg);
    }
    else if(type.find("Tetrahedron") != std::string::npos){
        createTetrahedron(fg);
    }
    else if(type.find("Octahedron") != std::string::npos){
        createOctahedron(fg);
    }
    else if(type.find("Icosahedron") != std::string::npos){
        createIcosahedron(fg);
    }
    else if(type.find("Dodecahedron") != std::string::npos){
        createDodecahedron(fg);
    }
    else if(type.find("BuckyBall") != std::string::npos){
        createbuckyball(fg);
    }
    else if(type == "Sphere"){
        int n = configuration[name]["n"].as_int_or_die();
        createSphere(fg,0,n);
    }
    else if(type == "Cone"){
        int n = configuration[name]["n"].as_int_or_die();
        double h = configuration[name]["height"].as_double_or_die();
        createCone(fg,n,h);
    }
    else if(type == "Cylinder"){
        int n = configuration[name]["n"].as_int_or_die();
        double h = configuration[name]["height"].as_double_or_die();
        createCylinder(fg,n,h);
    }
    else if(type == "Torus"){
        int n = configuration[name]["n"].as_int_or_die();
        int m = configuration[name]["m"].as_int_or_die();
        double r = configuration[name]["r"].as_double_or_die();
        double R = configuration[name]["R"].as_double_or_die();
        createTorus(fg,r,R,n,m);
    }
    else if(type == "3DLSystem"){
        string fileName = configuration[name]["inputfile"].as_string_or_die();
        double scale = configuration[name]["scale"].as_double_or_die();
        L3 lijnStysteem(fileName,fg);
    }

    for(int i = 0; i < pointCount; i++){
        vector<double> pointValue = configuration[name]["point" + to_string(i)].as_double_tuple_or_die();
        double x =pointValue[0];
        double y =pointValue[1];
        double z =pointValue[2];
        Vector3D* pointV = new Vector3D();
        fg->points.push_back(pointV->point(x,y,z));
    }

    for(int i = 0; i < lineCount; i++){
        vector<int> lineValue = configuration[name]["line" + to_string(i)].as_int_tuple_or_die();
        Face face;
        for(int a = 0; a < lineValue.size(); a++){
            face.point_indexes.push_back(lineValue[a]);
        }
        fg->faces.push_back(face);
    }

    Matrix* x = complete_matrix(scale,rotateX,rotateY, rotateZ,translateX,translateY,translateZ);

    fg->ambientReflection = ambR;
    fg->diffuseReflection = diffR;
    fg->specularReflection = specR;
    fg->reflectionCoefficient = reflectionCoefficient;

    if(type == "MengerSponge"){
        int iterations = configuration[name]["nrIterations"].as_int_or_die();
        generateMergeSpons(fg, iterations,rotateX,rotateY, rotateZ,translateX,translateY,translateZ, eyePointTransMatrix);
        return nullptr;
    }

    if(type.find("Fractal") != std::string::npos){
        double scale1 = configuration[name]["fractalScale"].as_double_or_die();
        int iterations = configuration[name]["nrIterations"].as_int_or_die();
        generateFractal(fg, iterations, scale1,rotateX,rotateY, rotateZ,translateX,translateY,translateZ, eyePointTransMatrix);
        return nullptr;
    }
    else{
        applyTransformation(fg, *x);
        applyTransformation(fg, eyePointTransMatrix);
        return fg;
    }
}

Figures3D* D3::getAllFigures(const ini::Configuration &configuration) {
    figures = new Figures3D();
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    for(int i = 0; i < nrFigures; i++){
        string name = "Figure" + to_string(i);
        //createFigure(configuration, name);
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D v3d;
        v3d.x=eye[0];
        v3d.y=eye[1];
        v3d.z=eye[2];
        eyeCordinaat = v3d;
        eyePointTransMatrix = eyePointTrans(v3d);

        Figure* fig = createFigure(configuration, name);
        if(fig != nullptr)
            figures->push_back(fig);
    }
    return figures;
}

Matrix D3::eyePointTrans(const Vector3D &eyepoint) {
    Matrix eyePoint;
    double theta;
    double phi;
    double r;
    toPolar(eyepoint, theta, phi, r);

    double cosTheta = (double)cos(theta);
    double sinTheta = (double)sin(theta);
    double cosPhi = (double)cos(phi);
    double sinPhi = (double)sin(phi);

    eyePoint(1,1) = -sinTheta;
    eyePoint(2,1) = cosTheta;

    eyePoint(1,2) = (-cosTheta) * cosPhi;
    eyePoint(2,2) = (-sinTheta) * cosPhi;
    eyePoint(3,2) = sinPhi;

    eyePoint(1,3) = cosTheta * sinPhi;
    eyePoint(2,3) = sinTheta * sinPhi;
    eyePoint(3,3) = cosPhi;
    eyePoint(4,3) = -r;
    return eyePoint;
}

void D3::toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt((pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2)));
    theta = atan2(point.y, point.x);
    phi = acos((point.z/r));
}

Lines2D D3::doProjection(const Figures3D* figuren3D) {
    Lines2D lines;

    for(auto it = figuren3D->begin(); it != figuren3D->end(); it++){
        for(auto it1 = (*it)->faces.begin(); it1 != (*it)->faces.end(); it1++){
            for(int i = 1; i < (it1)->point_indexes.size(); i++){
                Line2D lijn(doProjection((*it)->points[(*it1).point_indexes[i - 1]], 1),doProjection((*it)->points[(*it1).point_indexes[i]], 1), (*it)->ambientReflection);
                lines.push_back(lijn);
            }
            //einde terug naar begin!
            if((it1)->point_indexes.size() > 2){
                Line2D lijn(doProjection((*it)->points[(*it1).point_indexes[(*it1).point_indexes.size() - 1]], 1),doProjection((*it)->points[(*it1).point_indexes[0]], 1), (*it)->ambientReflection);
                lines.push_back(lijn);
            }
        }
    }
    return lines;
}

Point2D D3::doProjection(const Vector3D &point, const double d) {
    double x = (d*point.x)/(-point.z);
    double y = (d*point.y)/(-point.z);
    return Point2D(x,y, point.z);
}

void D3::TriangulateAll(const Figures3D* figuren3D) {
    for(auto it = figuren3D->begin(); it != figuren3D->end(); it++){
        list<Face> newFaces;
        for(auto it1 = (*it)->faces.begin(); it1 != (*it)->faces.end(); it1++){
            list<Face> newFace = triangulate((*it1));
            //TODO weg
            if(newFaces.size() == 0){
                newFaces = newFace;
            }
            else {
                if(newFace.size() != 0){
                    newFaces.insert(newFaces.begin(),newFace.begin(), newFace.end());
                }
                else{
                    newFaces.push_back(*it1);
                }
            }
        }
        (*it)->faces = newFaces;
    }
}

double D3::degree_to_rad(double degree)
{
    return (degree * (pi / 180));
}

D3::~D3() {
    for(auto it = figures->begin(); it != figures->end(); it++){
        delete *it;
    }
    delete figures;
    delete all;
}

Figure* D3::createKubus(Figure* fg) {
    fg->points = {Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),Vector3D::point(-1,-1,1),Vector3D::point(1,1,-1),Vector3D::point(-1,-1,-1),Vector3D::point(1,-1,1),Vector3D::point(-1,1,1)};
    fg->faces = {vector<int>({0,4,2,6}),vector<int>({4,1,7,2}),vector<int>({1,5,3,7}),vector<int>({5,0,6,3}),vector<int>({6,2,7,3}),vector<int>({0,5,1,4})};
    return fg;
}

Figure* D3::createTetrahedron(Figure* fg) {
    fg->points = {Vector3D::point(1,-1,-1),
                  Vector3D::point(-1,1,-1),
                  Vector3D::point(1,1,1),
                  Vector3D::point(-1,-1,1)};
    fg->faces = {vector<int>({0,1,2}),vector<int>({1,3,2}),vector<int>({0,3,1}),vector<int>({0,2,3})};
    return fg;
}

Figure* D3::createOctahedron(Figure* fg) {
    fg->points = {Vector3D::point(1,0,0),Vector3D::point(0,1,0),Vector3D::point(-1,0,0),Vector3D::point(0,-1,0),Vector3D::point(0,0,-1),Vector3D::point(0,0,1)};
    fg->faces = {vector<int>({0,1,5}),vector<int>({1,2,5}),vector<int>({2,3,5}),vector<int>({3,0,5}),vector<int>({1,0,4}),vector<int>({2,1,4}),vector<int>({3,2,4}),vector<int>({0,3,4})};
    return fg;
}

Figure* D3::createIcosahedron(Figure* fg) {
    fg->points.emplace_back(Vector3D::point(0,0,(sqrt(5)/2)));
    for(int i = 2; i <=6; i++){
        fg->points.emplace_back(Vector3D::point(cos(((i-2) * 2 * pi/5)),sin((i - 2) * 2 * pi / 5),0.5));
    }
    for(int i = 7; i <=11; i++){
        fg->points.emplace_back(Vector3D::point(cos(pi/5 + (i - 7) * 2 * pi/5),sin(pi/5 + (i - 7) * 2 * pi / 5), -0.5));
    }
    fg->points.emplace_back(Vector3D::point(0,0,(-sqrt(5)/2)));
    fg->faces = {vector<int>({0,1,2}),vector<int>({0,2,3}),vector<int>({0,3,4}),vector<int>({0,4,5}),vector<int>({0,5,1}),vector<int>({1,6,2}),vector<int>({2,6,7}),vector<int>({2,7,3}),vector<int>({3,7,8}),vector<int>({3,8,4}),vector<int>({4,8,9}),vector<int>({4,9,5}),vector<int>({5,9,10}),vector<int>({5,10,1}),vector<int>({1,10,6}),vector<int>({11,7,6}),vector<int>({11,8,7}),vector<int>({11,9,8}),vector<int>({11,10,9}),vector<int>({11,6,10})};
    return fg;
}

Figure *D3::createDodecahedron(Figure *fg) {
    Figure* icosahedron = new Figure();
    createIcosahedron(icosahedron);
    int teller = 0;
    for(auto it1 = icosahedron->faces.begin(); it1 != icosahedron->faces.end(); it1++){
        double newXPoint = (icosahedron->points[(*it1).point_indexes[0]].x + icosahedron->points[(*it1).point_indexes[1]].x + icosahedron->points[(*it1).point_indexes[2]].x)/3;
        double newYPoint = (icosahedron->points[(*it1).point_indexes[0]].y + icosahedron->points[(*it1).point_indexes[1]].y + icosahedron->points[(*it1).point_indexes[2]].y)/3;
        double newZPoint = (icosahedron->points[(*it1).point_indexes[0]].z + icosahedron->points[(*it1).point_indexes[1]].z + icosahedron->points[(*it1).point_indexes[2]].z)/3;
        fg->points.push_back(Vector3D::point(newXPoint, newYPoint, newZPoint));
        teller++;
    }
    //
    fg->faces = {vector<int>({0,1,2,3,4}), vector<int>({0,5,6,7,1}),vector<int>({1,7,8,9,2}),vector<int>({2,9,10,11,3}),vector<int>({3,11,12,13,4}),vector<int>({4,13,14,5,0}),vector<int>({19,18,17,16,15}),vector<int>({19,14,13,12,18}),vector<int>({18,12,11,10,17}),vector<int>({17,10,9,8,16}),vector<int>({16,8,7,6,15}),vector<int>({15,6,5,14,19})};
    delete icosahedron;
    return fg;
}

Figure *D3::createSphere(Figure *fg, const double radius, const int n) {
    Figure* icosahedron = new Figure();
    createIcosahedron(icosahedron);
    if(n == 0){
        *fg = *icosahedron;
    }
    for(int i = 0; i < n; i++) {
        fg->faces.clear();
        fg->points.clear();
        for (auto it1 = icosahedron->faces.begin(); it1 != icosahedron->faces.end(); it1++) {

            fg->points.push_back(icosahedron->points[(*it1).point_indexes[0]]);
            fg->points.push_back(icosahedron->points[(*it1).point_indexes[1]]);
            fg->points.push_back(icosahedron->points[(*it1).point_indexes[2]]);

            fg->points.push_back(Vector3D::point(((icosahedron->points[(*it1).point_indexes[0]].x + icosahedron->points[(*it1).point_indexes[1]].x)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[0]].y + icosahedron->points[(*it1).point_indexes[1]].y)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[0]].z + icosahedron->points[(*it1).point_indexes[1]].z)/2)));
            fg->points.push_back(Vector3D::point(((icosahedron->points[(*it1).point_indexes[0]].x + icosahedron->points[(*it1).point_indexes[2]].x)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[0]].y + icosahedron->points[(*it1).point_indexes[2]].y)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[0]].z + icosahedron->points[(*it1).point_indexes[2]].z)/2)));
            fg->points.push_back(Vector3D::point(((icosahedron->points[(*it1).point_indexes[1]].x + icosahedron->points[(*it1).point_indexes[2]].x)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[1]].y + icosahedron->points[(*it1).point_indexes[2]].y)/2),
                                                 ((icosahedron->points[(*it1).point_indexes[1]].z + icosahedron->points[(*it1).point_indexes[2]].z)/2)));
            int size = fg->points.size();
            int A = size - 6;
            int B = size - 5;
            int C = size - 4;
            int D = size - 3;
            int E = size - 2;
            int F = size - 1;

            fg->faces.push_back(vector<int>({A,D,E}));
            fg->faces.push_back(vector<int>({B,F,D}));
            fg->faces.push_back(vector<int>({C,E,F}));
            fg->faces.push_back(vector<int>({D,F,E}));
        }
        *icosahedron = *fg;
    }

    for(auto it = fg->points.begin(); it < fg->points.end(); it++){
        rescalePointSphere(*it);
    }

    delete icosahedron;
    return fg;
}

Vector3D D3::rescalePointSphere(Vector3D &point){
    double r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    point = Vector3D::point(point.x/r,point.y/r,point.z/r);
    return point;
}

Figure *D3::createCone(Figure *fg, const int n, const double h) {
    fg->points.push_back(Vector3D::point(0,0,h));//pn = top punt
    fg->faces.push_back(vector<int>());  //circel onderaan

    for(int i = 0; i <= n; i++){
        fg->points.push_back(Vector3D::point(cos((2*pi*i)/n),sin((2*pi*i)/n),0));
        fg->faces.push_back(vector<int>({(i),(i+1)%n,0}));//driehoek toevoegen aan faces
        (*fg->faces.begin()).point_indexes.push_back((n + 1) - i);//punt om circel te maken toevoegen aan faces[0]
    }
    fg->faces.push_back(vector<int>({((int)fg->points.size() - 3),((int)fg->points.size() - 2),0}));
    return fg;
}

Figure *D3::createCylinder(Figure *fg, const int n, const double h) {
    fg->faces.push_back(vector<int>());  //circel onderaan
    fg->faces.push_back(vector<int>());  //circel bovenaan
    for(int i = 2; i <= n + 1; i++){
        fg->points.push_back(Vector3D::point(cos((2*pi*i)/n),sin((2*pi*i)/n),0));//punt toevoegen voor onder aan
        int pos1 = fg->points.size() - 1;
        fg->points.push_back(Vector3D::point(cos((2*pi*i)/n),sin((2*pi*i)/n),h));   //punt toevoegen voor boven aan
        int pos2 = fg->points.size() - 1;
       // ,(pos2), (pos2 + 1)%n})
        if(!(fg->points.size() == 2 * n)){
            fg->faces.push_back(vector<int>({pos1 + 2, pos2 + 2,pos2, pos1}));//rechthoek toevoegen aan faces
        }
        else{
            fg->faces.push_back(vector<int>({0, 1,pos2, pos1}));//rechthoek toevoegen aan faces
        }

        (*fg->faces.begin()).point_indexes.push_back(pos1 );//punt om circel te maken toevoegen aan faces[0]
        (*std::next(fg->faces.begin())).point_indexes.push_back(pos2);//punt om circel te maken toevoegen aan faces[1]
    }
    return fg;
}

Figure *D3::createTorus(Figure *fg, const double r, const double R, const int n, const int m) {
    for(int i = 0; i < n; i++){
        double u = (2*i*pi)/n;
        for(int j = 0; j < m; j++){
            double v = (2*j*pi)/m;
            double Xuv = (R + r * cos(v))*cos(u);
            double Yuv = (R + r * cos(v))*sin(u);
            double Zuv = r * sin(v);
            fg->points.push_back(Vector3D::point(Xuv,Yuv,Zuv));

            int Pij = (i * m) + j;
            int Pij1 = (((i + 1)%n )* m) + j;
            int Pij2 = (((i + 1)%n )* m) + (j + 1)%m;
            int Pij3 = (i * m)  + (j + 1)%m;
            fg->faces.push_back(vector<int>({(Pij), (Pij1), (Pij2), (Pij3)}));
        }
    }
    return fg;
}

Figure *D3::createbuckyball(Figure *fg) {
    Figure* icosahedron = new Figure();
    createIcosahedron(icosahedron);
    int teller = 0;
    for(auto it1 = icosahedron->faces.begin(); it1 != icosahedron->faces.end(); it1++){
        /*
        fg->points.push_back(Vector3D::point(icosahedron->points[(*it1).point_indexes[0]]));
        fg->points.push_back(Vector3D::point(icosahedron->points[(*it1).point_indexes[1]]));
        fg->points.push_back(Vector3D::point(icosahedron->points[(*it1).point_indexes[2]]));
        */
        double u1 = (1.0/3.0);
        double u2 = (1.0/3.0)*2;
        fg->points.push_back(Vector3D::point((1-u1)*icosahedron->points[(*it1).point_indexes[0]].x + u1 * icosahedron->points[(*it1).point_indexes[1]].x,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[0]].y + u1 * icosahedron->points[(*it1).point_indexes[1]].y,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[0]].z + u1 * icosahedron->points[(*it1).point_indexes[1]].z
                                             ));
        fg->points.push_back(Vector3D::point((1-u2)*icosahedron->points[(*it1).point_indexes[0]].x + u2 * icosahedron->points[(*it1).point_indexes[1]].x,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[0]].y + u2 * icosahedron->points[(*it1).point_indexes[1]].y,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[0]].z + u2 * icosahedron->points[(*it1).point_indexes[1]].z
        ));


        fg->points.push_back(Vector3D::point((1-u1)*icosahedron->points[(*it1).point_indexes[1]].x + u1 * icosahedron->points[(*it1).point_indexes[2]].x,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[1]].y + u1 * icosahedron->points[(*it1).point_indexes[2]].y,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[1]].z + u1 * icosahedron->points[(*it1).point_indexes[2]].z
        ));
        fg->points.push_back(Vector3D::point((1-u2)*icosahedron->points[(*it1).point_indexes[1]].x + u2 * icosahedron->points[(*it1).point_indexes[2]].x,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[1]].y + u2 * icosahedron->points[(*it1).point_indexes[2]].y,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[1]].z + u2 * icosahedron->points[(*it1).point_indexes[2]].z
        ));


        fg->points.push_back(Vector3D::point((1-u1)*icosahedron->points[(*it1).point_indexes[2]].x + u1 * icosahedron->points[(*it1).point_indexes[0]].x,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[2]].y + u1 * icosahedron->points[(*it1).point_indexes[0]].y,
                                             (1-u1)*icosahedron->points[(*it1).point_indexes[2]].z + u1 * icosahedron->points[(*it1).point_indexes[0]].z
        ));
        fg->points.push_back(Vector3D::point((1-u2)*icosahedron->points[(*it1).point_indexes[2]].x + u2 * icosahedron->points[(*it1).point_indexes[0]].x,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[2]].y + u2 * icosahedron->points[(*it1).point_indexes[0]].y,
                                             (1-u2)*icosahedron->points[(*it1).point_indexes[2]].z + u2 * icosahedron->points[(*it1).point_indexes[0]].z
        ));
        fg->faces.push_back(vector<int>({teller + 7- 3,teller + 8- 3,teller + 3- 3,teller + 4- 3,teller + 5- 3,teller + 6 - 3}));

        teller += 6;
    }
    teller = 0;

    //Vijfhoeken
    fg->faces.push_back(vector<int>({12,18,24,0,6}));
    fg->faces.push_back(vector<int>({81,30,2,1,27}));
    fg->faces.push_back(vector<int>({4,3,33,41,8}));
    fg->faces.push_back(vector<int>({33,41,8,4,3}));
    fg->faces.push_back(vector<int>({48,54,14,10,9}));
    fg->faces.push_back(vector<int>({57,66,20,16,15}));
    fg->faces.push_back(vector<int>({69,77,26,22,21}));
    fg->faces.push_back(vector<int>({94,38,32,31,87}));
    fg->faces.push_back(vector<int>({91,50,44,40,39}));
    fg->faces.push_back(vector<int>({97,62,56,52,51}));
    fg->faces.push_back(vector<int>({103,74,68,64,63}));
    fg->faces.push_back(vector<int>({109,86,80,76,75}));
    fg->faces.push_back(vector<int>({102,108,114,90,96}));


    delete icosahedron;
    return fg;
}


D3::D3() {
    pi = 2*acos(0.0);
    Matrix schalen;
    Matrix x_as;
    Matrix y_as;
    Matrix z_as;
}

list<Face> D3::triangulate(const Face &face) {
    list<Face> faceTriangulate;
    int P0 = face.point_indexes[0];
    for(int i = 1; i < face.point_indexes.size() - 1; i++){
        faceTriangulate.push_back(Face(vector<int>({P0, face.point_indexes[i], face.point_indexes[i + 1]})));
    }
    return faceTriangulate;
}

void D3::generateFractal(Figure* fig,const int nr_iterations, const double scale,const double rotateX,const double rotateY, const double rotateZ, const double translateX,const double translateY,const double translateZ, Matrix& eyePoint) {
    if(nr_iterations > 0){
        for(int point = 0; point < fig->points.size(); point++){
            Figure* newFig = fig->copy();

            applyTransformation(newFig, complete_matrix_1(1/scale,rotateX,rotateY, rotateZ,translateX,translateY,translateZ));

            Vector3D newLocation = fig->points[point]-newFig->points[point];

            applyTransformation(newFig, complete_matrix_1(1,rotateX,rotateY, rotateZ,newLocation.x,newLocation.y,newLocation.z));

            generateFractal(newFig, nr_iterations - 1, scale,rotateX,rotateY, rotateZ,translateX,translateY,translateZ, eyePoint);
        }
        delete fig;
    }
    else{
        applyTransformation(fig, eyePoint);
        figures->push_back(fig);
    }
}

Matrix D3::complete_matrix_1(const double scaleFactor, const double angel_x, const double angel_y, const double angel_z,
                      const double transmatie_x, const double transmatie_y, const double transmatie_z) {
    Matrix all;
    all = (makeSchaalMatrix(scaleFactor)) * (makeX_asMatrix(angel_x) * makeY_asMatrix(angel_y) * makeZ_asMatrix(angel_z)) * (makeTranslatieMatrix(transmatie_x,transmatie_y,transmatie_z));
    return all;
}

void D3::generateMergeSpons(Figure* fig,const int nr_iterations,const double rotateX,const double rotateY, const double rotateZ, const double translateX,const double translateY,const double translateZ, Matrix& eyePoint) {
    if(nr_iterations > 0){
        Figure* temp = fig->copy();
        delete fig;

        Vector3D midden = (temp->points[0] + temp->points[7])/2;
        applyTransformation(temp, makeTranslatieMatrix(midden.x, midden.y, midden.z));
        applyTransformation(temp, makeSchaalMatrix(1.0/3.0));
        Vector3D midden1 = (temp->points[0] + temp->points[7])/2;
        Vector3D move1 = midden - midden1;
        applyTransformation(temp, makeTranslatieMatrix(move1.x, move1.y, move1.z));

        vector<Vector3D> alreadyTransform;
        for(int i = 0; i < temp->points.size(); i++){
            for(int i1 = 0; i1 < temp->points.size(); i1++){
                Figure* fig1 = temp->copy();
                Vector3D newLocation = temp->points[i] - temp->points[i1];

                if(!(newLocation.z == 0 && newLocation.x == 0 ||
                   newLocation.z == 0 && newLocation.y == 0 ||
                   newLocation.x == 0 && newLocation.y == 0 ) &&
                   find(alreadyTransform.begin(), alreadyTransform.end(),newLocation) == alreadyTransform.end()){
                    applyTransformation(fig1, complete_matrix_1(1,0,0, 0,newLocation.x,newLocation.y,newLocation.z));
                    generateMergeSpons(fig1, nr_iterations - 1,0,0, 0,0,0,0, eyePoint);
                    alreadyTransform.push_back(newLocation);
                }

            }
        }
        delete temp;
    }
    else{
        applyTransformation(fig, eyePoint);
        figures->push_back(fig);
    }
}
