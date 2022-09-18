#include "easy_image.h"
#include "ini_configuration.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "L2.h"
#include "D3.h"
#include <chrono>
#include <ctime>
#include <thread>
#include "Light.h"
using namespace img;



Lights3D getAllLights(const ini::Configuration &configuration, Matrix &eyePointMatrix){
    Lights3D lights;
    int nrFigures = configuration["General"]["nrLights"].as_int_or_die();
    for(int i = 0; i < nrFigures; i++){
        string name = "Light" + to_string(i);
        vector<double> ambientLight = configuration[name]["ambientLight"].as_double_tuple_or_default(vector<double>({0,0,0}));
        vector<double> diffuseLight = configuration[name]["diffuseLight"].as_double_tuple_or_default(vector<double>({0,0,0}));
        vector<double> specularLight = configuration[name]["specularLight"].as_double_tuple_or_default(vector<double>({0,0,0}));

        bool infinity = configuration[name]["infinity"].as_bool_or_default(false);

        if(infinity){
            vector<double> direction = configuration[name]["direction"].as_double_tuple_or_die();
            Vector3D directionVector;
            InfLight* infLight = new InfLight(img::Color(ambientLight[0], ambientLight[1], ambientLight[2]),
                                              img::Color(diffuseLight[0], diffuseLight[1], diffuseLight[2]),
                                              img::Color(specularLight[0], specularLight[1], specularLight[2]),
                                              (directionVector.vector(direction[0], direction[1], direction[2]) * eyePointMatrix));
            lights.push_back(infLight);
        }
        else if(configuration[name]["location"].exists()){
            vector<double> location = configuration[name]["location"].as_double_tuple_or_die();
            double spotAngle = configuration[name]["spotAngle"].as_double_or_default(360);
            PointLight* pointLight = new PointLight(img::Color(ambientLight[0], ambientLight[1], ambientLight[2]),
                                                     img::Color(diffuseLight[0], diffuseLight[1], diffuseLight[2]),
                                                    img::Color(specularLight[0], specularLight[1], specularLight[2]),
                                                        (Vector3D::point(location[0], location[1], location[2]) * eyePointMatrix), spotAngle);
            lights.push_back(pointLight);
        }
        else{
            Light* light = new Light(img::Color(ambientLight[0], ambientLight[1], ambientLight[2]),
                                     img::Color(diffuseLight[0], diffuseLight[1], diffuseLight[2]),
                                     img::Color(specularLight[0], specularLight[1], specularLight[2]));
            lights.push_back(light);
        }
    }
    return lights;
}

Lights3D getDefautlLight(){
    Lights3D lights;
    Light* ambient = new Light(img::Color(1.0, 1.0, 1.0), img::Color(0.0, 0.0, 0.0),img::Color(0.0, 0.0, 0.0));
    lights.push_back(ambient);
    return lights;
}

EasyImage draw2DLines(img::Lines2D &lines, const int size, Color &backgroundColor, bool zBufferLines = false) {
    double xmin = numeric_limits<double>::infinity(), xmax = 0, ymin = +numeric_limits<double>::infinity(), ymax = 0;
    for(Lines2D::iterator it = lines.begin(); it != lines.end(); it++) {
        double x1 = it->getPoint1().getX();
        double x2 = it->getPoint2().getX();

        double y1 = it->getPoint1().getY();
        double y2 = it->getPoint2().getY();

        if(xmin > x1) xmin = x1;
        if(xmax < x1) xmax = x1;

        if(xmin > x2) xmin = x2;
        if(xmax < x2) xmax = x2;

        if(ymin > y1) ymin = y1;
        if(ymax < y1) ymax = y1;

        if(ymin > y2) ymin = y2;
        if(ymax < y2) ymax = y2;
    }

    double xrange = xmax - xmin ;
    double yrange = ymax - ymin ;
    double imageX = size * (xrange / max(xrange, yrange));
    double imageY = size * (yrange / max(xrange, yrange));

    double d = 0.95 * (imageX/xrange);

    EasyImage image(imageX, imageY, backgroundColor);

    ZBuffer z = ZBuffer(imageX,imageY);
    double DCx = d * ((xmin+xmax)/2);
    double DCy = d * ((ymin+ymax)/2);
    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;
    for(Lines2D::iterator it = lines.begin(); it != lines.end(); it++){
        if(zBufferLines){
            image.draw_zbuf_line(z,
                                 round((it->getPoint1().getX() * d) + dx),
                                 round((it->getPoint1().getY() * d) + dy),
                                 it->getPoint1().getZ(),
                                 round((it->getPoint2().getX() * d) + dx),
                                 round((it->getPoint2().getY() * d) + dy),
                                 it->getPoint2().getZ(),
                                 it->getColor());
        }
        else{
            image.draw_line(round((it->getPoint1().getX() * d) + dx), round((it->getPoint1().getY() * d) + dy), round((it->getPoint2().getX() * d) + dx), round((it->getPoint2().getY() * d) + dy) ,it->getColor());
        }
    }
    return image;
}

EasyImage drawTriangulateFaces(D3& d3Image,Figures3D* figures3d, Lights3D& lights, const int size, Color &backgroundColor) {
    Lines2D lines = d3Image.doProjection(figures3d);
    double xmin = numeric_limits<double>::infinity(), xmax = 0, ymin = +numeric_limits<double>::infinity(), ymax = 0;
    for(Lines2D::iterator it = lines.begin(); it != lines.end(); it++) {
        double x1 = it->getPoint1().getX();
        double x2 = it->getPoint2().getX();

        double y1 = it->getPoint1().getY();
        double y2 = it->getPoint2().getY();

        if(xmin > x1) xmin = x1;
        if(xmax < x1) xmax = x1;

        if(xmin > x2) xmin = x2;
        if(xmax < x2) xmax = x2;

        if(ymin > y1) ymin = y1;
        if(ymax < y1) ymax = y1;

        if(ymin > y2) ymin = y2;
        if(ymax < y2) ymax = y2;
    }

    double xrange = xmax - xmin;
    double yrange = ymax - ymin;
    double imageX = size * (xrange / max(xrange, yrange));
    double imageY = size * (yrange / max(xrange, yrange));

    double d = 0.95 * (imageX/xrange);

    d3Image.TriangulateAll(figures3d);

    EasyImage image(imageX, imageY, backgroundColor);

    ZBuffer z = ZBuffer(imageX,imageY);
    double DCx = d * ((xmin+xmax)/2);
    double DCy = d * ((ymin+ymax)/2);
    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;
    Figures3D* figuren = figures3d;
    double red = 0.3;
    double green = 0.0;
    double blue = 0.0;

    vector<vector<double>> colors = {{1,1,1}
                                    , {0,1,0}
                                    , {0,0,1}};
    int i = 0;
    for(Figures3D::iterator it = figuren->begin(); it != figuren->end(); it++) {

        Color ambientColor(0.0,0.0,0.0);

        for(auto it1 = lights.begin(); it1 != lights.end(); it1++) {
            double redAmbient = (*it)->ambientReflection.getRedDouble() * (*it1)->ambientLight.getRedDouble();
            double greenAmbient = (*it)->ambientReflection.getGreenDouble() * (*it1)->ambientLight.getGreenDouble();
            double blueAmbient = (*it)->ambientReflection.getBlueDouble() * (*it1)->ambientLight.getBlueDouble();
            ambientColor.setRed(ambientColor.getRedDouble() + redAmbient);
            ambientColor.setGreen(ambientColor.getGreenDouble() + greenAmbient);
            ambientColor.setBlue(ambientColor.getBlueDouble() + blueAmbient);
        }

        for(std::list<Face>:: iterator it1 = (*it)->faces.begin(); it1 != (*it)->faces.end(); it1++){
            image.draw_zbuf_triag(z, (*it)->points[(*it1).point_indexes[0]],
                                  (*it)->points[(*it1).point_indexes[1]],
                                  (*it)->points[(*it1).point_indexes[2]],
                                  d,
                                  dx,
                                  dy,
                                  ambientColor,
                                  (*it)->diffuseReflection,
                                  (*it)->specularReflection,
                                  lights,
                                  (*it)->reflectionCoefficient);
        }
    }

    for(auto it = lights.begin(); it != lights.end(); it++) {
        delete (*it);
    }

    return image;
}

img::EasyImage DLSystem2(const ini::Configuration &configuration, int size){
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    vector<double> color = configuration["2DLSystem"]["color"].as_double_tuple_or_die();

    string d2LFilePad = configuration["2DLSystem"]["inputfile"].as_string_or_die();
    Color lineColor(color[0], color[1], color[2]);
    L2D l2d(d2LFilePad, lineColor);
    Lines2D lijnen = l2d.generateList();

    Color backroundColor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return draw2DLines(lijnen,size, backroundColor);
}

img::EasyImage Wireframe(const ini::Configuration &configuration, int size){
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    D3 d3Image;
    Figures3D* figures3d = d3Image.getAllFigures(configuration);
    Lines2D lijnen = d3Image.doProjection(figures3d);
    Color backroundColor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return draw2DLines(lijnen,size, backroundColor);
}

img::EasyImage ZBufferedWireframe(const ini::Configuration &configuration, int size){
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    D3 d3Image;
    Figures3D* figures3d = d3Image.getAllFigures(configuration);
    //d3Image.TriangulateAll(figures3d);
    Lines2D lijnen = d3Image.doProjection(figures3d);
    Color backroundColor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return draw2DLines(lijnen,size, backroundColor, true);
}

img::EasyImage ZBuffering(const ini::Configuration &configuration, int size){
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    D3 d3Image;

    Lights3D lights = getDefautlLight();
    Figures3D* figures3d = d3Image.getAllFigures(configuration);
    Color backroundColor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return drawTriangulateFaces(d3Image,figures3d,lights, size, backroundColor);
}

img::EasyImage LightedZBuffering(const ini::Configuration &configuration, int size){
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    D3 d3Image;
    Figures3D* figures3d = d3Image.getAllFigures(configuration);
    Lights3D lights = getAllLights(configuration, d3Image.eyePointTransMatrix);
    Color backroundColor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return drawTriangulateFaces(d3Image,figures3d,lights, size, backroundColor);
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    int size = configuration["General"]["size"].as_int_or_die();
    string type = configuration["General"]["type"].as_string_or_die();

    if(type == "2DLSystem"){
        return DLSystem2(configuration, size);
    }
    else if(type == "Wireframe"){
        return Wireframe(configuration, size);
    }
    else if(type == "ZBufferedWireframe"){
        return ZBufferedWireframe(configuration, size);
    }
    else if(type == "ZBuffering"){
        return ZBuffering(configuration, size);
    }
    else if(type == "LightedZBuffering"){
        return LightedZBuffering(configuration, size);
    }
    return img::EasyImage();
}

int main(int argc, char const* argv[])
{
    std::clock_t start;
    double duration;

    start = std::clock();

    int retVal = 0;
    try
    {
        std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string fileName;
            while (std::getline(fileIn, fileName)) {
                args.push_back(fileName);
            }
        }
        vector<thread*> activeThreads;
        for(std::string fileName : args) {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName.c_str());
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
            }
            img::EasyImage image = generate_image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    string format = ".bmp";
                    fileName.append(format);
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;
                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"printf: "<< duration <<'\n';

    return retVal;
}

