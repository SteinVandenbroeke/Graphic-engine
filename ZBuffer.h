//
// Created by stein on 22/03/2021.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H


#include <vector>
#include <limits>

class ZBuffer: public std::vector<std::vector<double>>
{
public:
    ~ZBuffer();
    ZBuffer(const int width, const int height);
    vector<std::vector<double>> v;
};


#endif //ENGINE_ZBUFFER_H
