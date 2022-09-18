//
// Created by stein on 22/03/2021.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {
    double inf = std::numeric_limits<double>::infinity();
    v= vector<std::vector<double>>(width, vector<double>(height, inf));
}

ZBuffer::~ZBuffer() {
}
