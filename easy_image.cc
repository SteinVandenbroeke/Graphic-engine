/*
 * easy_image.cc
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "easy_image.h"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "Light.h"
#include <chrono>
#include <thread>

#ifndef le32toh
#define le32toh(x) (x)
#endif

namespace
{
	//structs borrowed from wikipedia's article on the BMP file format
	struct bmpfile_magic
	{
			uint8_t magic[2];
	};

	struct bmpfile_header
	{
			uint32_t file_size;
			uint16_t reserved_1;
			uint16_t reserved_2;
			uint32_t bmp_offset;
	};
	struct bmp_header
	{
			uint32_t header_size;
			int32_t width;
			int32_t height;
			uint16_t nplanes;
			uint16_t bits_per_pixel;
			uint32_t compress_type;
			uint32_t pixel_size;
			int32_t hres;
			int32_t vres;
			uint32_t ncolors;
			uint32_t nimpcolors;
	};
	//copy-pasted from lparser.cc to allow these classes to be used independently from each other
	class enable_exceptions
	{
		private:
			std::ios& ios;
			std::ios::iostate state;
		public:
			enable_exceptions(std::ios& an_ios, std::ios::iostate exceptions) :
				ios(an_ios)
			{
				state = ios.exceptions();
				ios.exceptions(exceptions);
			}
			~enable_exceptions()
			{
				ios.exceptions(state);
			}
	};
	//helper function to convert a number (char, int, ...) to little endian
	//regardless of the endiannes of the system
	//more efficient machine-dependent functions exist, but this one is more portable
	template<typename T> T to_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			temp_storage.bytes[i] = value & 0xFF;
			value >>= 8;
		}
		return temp_storage.t;
	}

	template<typename T> T from_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;
		temp_storage.t = value;
		T retVal = 0;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			retVal = (retVal << 8) | temp_storage.bytes[sizeof(T) - i - 1];
		}
		return retVal;
	}

}
img::Color::Color() :
	blue(0), green(0), red(0)
{
}
img::Color::Color(uint8_t r, uint8_t g, uint8_t b) :
	blue(b), green(g), red(r)
{
}
img::Color::~Color()
{
}

img::Color::Color(double r, double g, double b) {
    setRed(r);
    setGreen(g);
    setBlue(b);
}

void img::Color::setRed(double code) {
    if(code > 1){
        code = 1;
    }
    red = round(code * 255);
}

void img::Color::setGreen(double code) {
    if(code > 1){
        code = 1;
    }
    green = round(code * 255);
}

void img::Color::setBlue(double code) {
    if(code > 1){
        code = 1;
    }
    blue = round(code * 255);
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(std::string const& msg) :
	message(msg)
{
}
img::UnsupportedFileTypeException::UnsupportedFileTypeException(const UnsupportedFileTypeException &original)
: std::exception(original)
, message(original.message)
{
}
img::UnsupportedFileTypeException::~UnsupportedFileTypeException() throw ()
{
}
img::UnsupportedFileTypeException& img::UnsupportedFileTypeException::operator=(UnsupportedFileTypeException const& original)
{
	this->message = original.message;
	return *this;
}
const char* img::UnsupportedFileTypeException::what() const throw ()
{
	return message.c_str();
}

img::EasyImage::EasyImage() :
	width(0), height(0), bitmap()
{
}

img::EasyImage::EasyImage(unsigned int _width, unsigned int _height, Color color) :
	width(_width), height(_height), bitmap(width * height, color)
{
}

img::EasyImage::EasyImage(EasyImage const& img) :
        width(img.width), height(img.height), bitmap(img.bitmap)
{
}

img::EasyImage::~EasyImage()
{
	bitmap.clear();
}

img::EasyImage& img::EasyImage::operator=(img::EasyImage const& img)
{
	width = img.width;
	height = img.height;
	bitmap.assign(img.bitmap.begin(),img.bitmap.end());
	return (*this);
}

unsigned int img::EasyImage::get_width() const
{
	return width;
}

unsigned int img::EasyImage::get_height() const
{
	return height;
}

void img::EasyImage::clear(Color color)
{
	for (std::vector<Color>::iterator i = bitmap.begin(); i != bitmap.end(); i++)
	{
		*i = color;
	}
}

img::Color& img::EasyImage::operator()(unsigned int x, unsigned int y)
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

img::Color const& img::EasyImage::operator()(unsigned int x, unsigned int y) const
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

void img::EasyImage::draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color)
{
	assert(x0 < this->width && y0 < this->height);
	assert(x1 < this->width && y1 < this->height);
	if (x0 == x1)
	{
		//special case for x0 == x1
		for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
		{
			(*this)(x0, i) = color;
		}
	}
	else if (y0 == y1)
	{
		//special case for y0 == y1
		for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
		{
			(*this)(i, y0) = color;
		}
	}
	else
	{
		if (x0 > x1)
		{
			//flip points if x1>x0: we want x0 to have the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
		if (-1.0 <= m && m <= 1.0)
		{
			for (unsigned int i = 0; i <= (x1 - x0); i++)
			{
				(*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
			}
		}
		else if (m > 1.0)
		{
			for (unsigned int i = 0; i <= (y1 - y0); i++)
			{
				(*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
			}
		}
		else if (m < -1.0)
		{
			for (unsigned int i = 0; i <= (y0 - y1); i++)
			{
				(*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
			}
		}
	}
}

void img::EasyImage::draw_zbuf_line(ZBuffer &f, unsigned int x0, unsigned int y0, double z0,
                                    unsigned int x1, unsigned int y1, double z1, img::Color color) {
    assert(x0 < this->width && y0 < this->height);
    assert(x1 < this->width && y1 < this->height);

    if (x0 == x1)
    {
        //special case for x0 == x1
        //150 - 145 = 5
        double min = std::min(y0, y1);
        double max = std::max(y0, y1);
        double a = max - min;
        for (unsigned int i = min; i <= max; i++)
        {
            double p = (a - (i - min))/a;
            double zi = p/z0 + (1 - p)/z1;
            if((zi) < (f.v[x0][i])){
                (*this)(x0, i) = color;
                f.v[x0][i] = zi;
            }
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1c
        double min = std::min(x0, x1);
        double max = std::max(x0, x1);
        double a = max - min;
        for (unsigned int i = min; i <= max; i++)
        {
            double p = (a - (i - min))/a;
            double zi = p/z0 + (1 - p)/z1;
            if((zi) < (f.v[i][y0])){
                (*this)(i, y0) = color;
                f.v[i][y0] = zi;
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            //std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            double a = (x1 - x0);
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double p = (a - i)/a;
                double zi = p/z0 + (1 - p)/z1;
                if((zi) < (f.v[x0 + i][(unsigned int) round(y0 + m * i)])){
                    (*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                    f.v[x0 + i][(unsigned int) round(y0 + m * i)] = zi;
                }
            }
        }
        else if (m > 1.0)
        {
            double a = (y1 - y0);
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double p = (a - i)/a;
                double zi = p/z0 + (1 - p)/z1;
                if((zi) < (f.v[(unsigned int) round(x0 + (i / m))][y0 + i])){
                    (*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                    f.v[(unsigned int) round(x0 + (i / m))][y0 + i] = zi;
                }
            }
        }
        else if (m < -1.0)
        {
            double a = (y0 - y1);
            for (unsigned int i = 0; i <= a; i++)
            {
                double p = (a - i)/a;
                double zi = p/z0 + (1 - p)/z1;
                if((zi) < (f.v[(unsigned int) round(x0 - (i / m))][y0 - i])){
                    (*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                    f.v[(unsigned int) round(x0 - (i / m))][y0 - i] = zi;
                }
            }
        }
    }
}


void img::EasyImage:: draw_zbuf_triag(ZBuffer& zBuffer, Vector3D const& A, Vector3D const& B, Vector3D const& C, double d, double dx, double dy, img::Color ambientColor, img::Color diffuseReflection,img::Color &specularLight, std::list<Light*>& lights, double ms) {
    double Axp = (d*A.x)/(-A.z)+dx;
    double Ayp = (d*A.y)/(-A.z)+dy;
    Point2D Ap = Point2D(Axp, Ayp);

    double Bxp = (d*B.x)/(-B.z)+dx;
    double Byp = (d*B.y)/(-B.z)+dy;
    Point2D Bp = Point2D(Bxp, Byp);

    double Cxp = (d*C.x)/(-C.z)+dx;
    double Cyp = (d*C.y)/(-C.z)+dy;
    Point2D Cp = Point2D(Cxp, Cyp);

    Line2D AB(Ap, Bp);
    Line2D AC(Ap, Cp);
    Line2D BC(Bp, Cp);

    double yMinTemp = min(Ayp, Byp);
    int yMin = round(min(yMinTemp, Cyp)+ 0.5);

    double yMaxTemp = max(Ayp, Byp);
    int yMax = round(max(yMaxTemp, Cyp)- 0.5);


    Vector3D u = B - A;
    Vector3D v = C - A;
    double w1 = (u.y * v.z) - (u.z * v.y);
    double w2 = (u.z * v.x) - (u.x * v.z);
    double w3 = (u.x * v.y) - (u.y * v.x);
    double k = w1 * A.x + w2 * A.y + w3 * A.z;
    double dzdx = w1/(-(d*k));
    double dzdy = w2/(-(d*k));

    Color colorAfterLights = calculateDiffuseInfLightColor(ambientColor,diffuseReflection, specularLight, lights, Vector3D::cross(u, v), ms);


    for(int yi = yMin; yi <= yMax; yi++){
        double XLab, XLac, XLbc;
        double XRab, XRac, XRbc;
        XLab = XLac = XLbc = std::numeric_limits<int>::max();
        XRab = XRac = XRbc = -std::numeric_limits<int>::max();

        if((yi - AB.getPoint1().getY())*(yi - AB.getPoint2().getY()) <= 0 && AB.getPoint1().getY() != AB.getPoint2().getY()){
            Point2D P = AB.getPoint1();
            Point2D Q = AB.getPoint2();
            double xi = Q.getX() + (P.getX() - Q.getX())*((yi-Q.getY())/(P.getY()-Q.getY()));
            XLab = XRab = xi;
        }
        if((yi - AC.getPoint1().getY())*(yi - AC.getPoint2().getY()) <= 0 && AC.getPoint1().getY() != AC.getPoint2().getY()){
            Point2D P = AC.getPoint1();
            Point2D Q = AC.getPoint2();
            double xi = Q.getX() + (P.getX() - Q.getX())*((yi-Q.getY())/(P.getY()-Q.getY()));
            XLac = XRac = xi;
        }
        if((yi - BC.getPoint1().getY())*(yi - BC.getPoint2().getY()) <= 0 && BC.getPoint1().getY() != BC.getPoint2().getY()){
            Point2D P = BC.getPoint1();
            Point2D Q = BC.getPoint2();
            double xi = Q.getX() + (P.getX() - Q.getX())*((yi-Q.getY())/(P.getY()-Q.getY()));
            XLbc = XRbc = xi;
        }
        double xlTemp = (min(XLab, XLac));
        int xl = round(min(xlTemp, XLbc) + 0.5);

        double xrTemp = (max(XRab, XRac));
        int xr = round(max(xrTemp, XRbc) - 0.5);

        double xg = (Axp + Bxp + Cxp)/3;
        double yg = (Ayp + Byp + Cyp)/3;
        double oneOverzg = 1/(3 * A.z) + 1/(3 * B.z) + 1/(3 * C.z);

        for(int x = xl; x <= xr; x++){
            double oneOverZ = 1.0001 * oneOverzg + (x - xg)*dzdx + (yi - yg)*dzdy;

            if(oneOverZ < zBuffer.v[x][yi]){
                Color colorObject = calculateDiffusePointLightColor(colorAfterLights,diffuseReflection, specularLight, lights, Vector3D::cross(u, v), x -dx,yi - dy,oneOverZ,d, ms);
                (*this)(x, yi) = colorObject;
                zBuffer.v[x][yi] = oneOverZ;
            }
        }
    }
}

std::ostream& img::operator<<(std::ostream& out, EasyImage const& image)
{

	//temporaryily enable exceptions on output stream
	enable_exceptions(out, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need:
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	uint8_t padding[] =
	{ 0, 0, 0, 0 };
	//calculate the total size of the pixel data
	unsigned int line_width = image.get_width() * 3; //3 bytes per pixel
	unsigned int line_padding = 0;
	if (line_width % 4 != 0)
	{
		line_padding = 4 - (line_width % 4);
	}
	//lines must be aligned to a multiple of 4 bytes
	line_width += line_padding;
	unsigned int pixel_size = image.get_height() * line_width;

	//start filling the headers
	magic.magic[0] = 'B';
	magic.magic[1] = 'M';

	file_header.file_size = to_little_endian(pixel_size + sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.bmp_offset = to_little_endian(sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.reserved_1 = 0;
	file_header.reserved_2 = 0;
	header.header_size = to_little_endian(sizeof(header));
	header.width = to_little_endian(image.get_width());
	header.height = to_little_endian(image.get_height());
	header.nplanes = to_little_endian(1);
	header.bits_per_pixel = to_little_endian(24);//3bytes or 24 bits per pixel
	header.compress_type = 0; //no compression
	header.pixel_size = pixel_size;
	header.hres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.vres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.ncolors = 0; //no color palette
	header.nimpcolors = 0;//no important colors

	//okay that should be all the header stuff: let's write it to the stream
	out.write((char*) &magic, sizeof(magic));
	out.write((char*) &file_header, sizeof(file_header));
	out.write((char*) &header, sizeof(header));

	//okay let's write the pixels themselves:
	//they are arranged left->right, bottom->top, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red they should be written automatically
			//in the right order
			out.write((char*) &image(j, i), 3 * sizeof(uint8_t));
		}
		if (line_padding > 0)
			out.write((char*) padding, line_padding);
	}
	//okay we should be done
	return out;
}
std::istream& img::operator>>(std::istream& in, EasyImage & image)
{
	enable_exceptions(in, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	//a temp buffer for reading the padding at the end of each line
	uint8_t padding[] =
	{ 0, 0, 0, 0 };

	//read the headers && do some sanity checks
	in.read((char*) &magic, sizeof(magic));
	if (magic.magic[0] != 'B' || magic.magic[1] != 'M')
		throw UnsupportedFileTypeException("Could not parse BMP File: invalid magic header");
	in.read((char*) &file_header, sizeof(file_header));
	in.read((char*) &header, sizeof(header));
	if (le32toh(header.pixel_size) + le32toh(file_header.bmp_offset) != le32toh(file_header.file_size))
		throw UnsupportedFileTypeException("Could not parse BMP File: file size mismatch");
	if (le32toh(header.header_size) != sizeof(header))
		throw UnsupportedFileTypeException("Could not parse BMP File: Unsupported BITMAPV5HEADER size");
	if (le32toh(header.compress_type) != 0)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only uncompressed BMP files can be parsed");
	if (le32toh(header.nplanes) != 1)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only one plane should exist in the BMP file");
	if (le32toh(header.bits_per_pixel) != 24)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only 24bit/pixel BMP's are supported");
	//if height<0 -> read top to bottom instead of bottom to top
	bool invertedLines = from_little_endian(header.height) < 0;
	image.height = std::abs(from_little_endian(header.height));
	image.width = std::abs(from_little_endian(header.width));
	unsigned int line_padding = from_little_endian(header.pixel_size) / image.height - (3 * image.width);
	//re-initialize the image bitmap
	image.bitmap.clear();
	image.bitmap.assign(image.height * image.width, Color());
	//okay let's read the pixels themselves:
	//they are arranged left->right., bottom->top if height>0, top->bottom if height<0, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red, the data read should be written in the right variables
			if (invertedLines)
			{
				//store top-to-bottom
				in.read((char*) &image(j, image.height - 1 - i), 3 * sizeof(uint8_t));
			}
			else
			{
				//store bottom-to-top
				in.read((char*) &image(j, i), 3 * sizeof(uint8_t));
			}
		}
		if (line_padding > 0)
		{
			in.read((char*) padding, line_padding);
		}
	}
	//okay we're done
	return in;
}

img::Color img::EasyImage::calculateDiffuseInfLightColor(img::Color &ambientColor,img::Color &diffuseReflection ,img::Color &specularLight , list<Light*> &lights, Vector3D w, double ms) {
    Color newColor(0.0,0.0,0.0);
    Vector3D n = Vector3D::normalise(w);
    for(auto it = lights.begin(); it != lights.end(); it++){
        double redDiffuse = 0;
        double greenDiffuse = 0;
        double blueDiffuse= 0;
        double redReflextionPoint = 0;
        double greenReflextionPoint = 0;
        double blueReflextionPoint = 0;
        if((*it)->isInfLight()) {
            double cosA = (n.x * (*it)->getlVector().x) + (n.y * (*it)->getlVector().y) + (n.z * (*it)->getlVector().z);
            if(cosA > 0){
                redDiffuse = (diffuseReflection.getRedDouble() * (*it)->diffuseLight.getRedDouble()) * cosA;
                greenDiffuse = (diffuseReflection.getGreenDouble() * (*it)->diffuseLight.getGreenDouble()) * cosA;
                blueDiffuse = (diffuseReflection.getBlueDouble() * (*it)->diffuseLight.getBlueDouble()) * cosA;
            }
        }

        newColor.setRed(newColor.getRedDouble() + redDiffuse);
        newColor.setGreen(newColor.getGreenDouble() + greenDiffuse);
        newColor.setBlue(newColor.getBlueDouble() + blueDiffuse);
    }
    newColor.setRed(newColor.getRedDouble() + ambientColor.getRedDouble());
    newColor.setGreen(newColor.getGreenDouble() + ambientColor.getGreenDouble());
    newColor.setBlue(newColor.getBlueDouble() + ambientColor.getBlueDouble());


    return newColor;
}

img::Color img::EasyImage::calculateDiffusePointLightColor(img::Color &ambientColorAnddiffuseInfColor,img::Color &diffuseReflection,img::Color &specularLight ,
                                                           list<Light*> &lights, Vector3D w, double x, double y, double z, double d, double ms) {
    Color newColor(0.0,0.0,0.0);
    Vector3D n = Vector3D::normalise(w);

    double Ze = 1/(double)z;
    double Xe = (x * (-Ze))/(d);
    double Ye = (y * (-Ze))/(d);
    Vector3D PointP = Vector3D::point(Xe, Ye, Ze);

    for(auto it = lights.begin(); it != lights.end(); it++){
        double redDiffusePoint = 0;
        double greenDiffusePoint = 0;
        double blueDiffusePoint= 0;
        double redReflextionPoint = 0;
        double greenReflextionPoint = 0;
        double blueReflextionPoint = 0;

        if((*it)->isPointLight()) {
            double As = ((*it)->getspotAngle()) * (3.14159265359/180);//0.01745329251 = pi / 180
            Vector3D location = (*it)->getLocation();
            Vector3D vectorL = Vector3D::normalise(location - PointP);

            double cosA = (n.x * vectorL.x) + (n.y * vectorL.y) + (n.z * vectorL.z);
            if(cosA > 0){
                //cout << cosA << ">" << As << endl;
                double calc = 0;
                if((*it)->getspotAngle() != 360 && cosA > cos(As)){
                    calc = 1-(1-cosA)/(1-cos(As));
                }
                else if((*it)->getspotAngle() == 360){
                    calc = cosA;
                }

                redDiffusePoint = diffuseReflection.getRedDouble() * (*it)->diffuseLight.getRedDouble() * calc;
                greenDiffusePoint = diffuseReflection.getGreenDouble() * (*it)->diffuseLight.getGreenDouble() * calc;
                blueDiffusePoint = diffuseReflection.getBlueDouble() * (*it)->diffuseLight.getBlueDouble() * calc;


                Vector3D r  = Vector3D::normalise(2 * n * cosA - vectorL);
                Vector3D camera = Vector3D::normalise(-PointP);
                double cosB = (r.x * camera.x) + (r.y * camera.y) + (r.z * camera.z);

                if(cosB > 0){
                    // cout << specularLight.getRedDouble() << "*" <<(*it)->specularLight.getRedDouble() << "*"  << cosB << endl;
                    cosB = pow(cosB,ms);
                    redReflextionPoint = specularLight.getRedDouble() * (*it)->specularLight.getRedDouble() * cosB;
                    greenReflextionPoint = specularLight.getGreenDouble() * (*it)->specularLight.getGreenDouble() * cosB;
                    blueReflextionPoint = specularLight.getBlueDouble() * (*it)->specularLight.getBlueDouble() * cosB;
                }
            }
            newColor.setRed(newColor.getRedDouble() + redDiffusePoint + redReflextionPoint);
            newColor.setGreen(newColor.getGreenDouble() + greenDiffusePoint + greenReflextionPoint);
            newColor.setBlue(newColor.getBlueDouble() + blueDiffusePoint + blueReflextionPoint);
        }
        else if((*it)->isInfLight()){
            Vector3D vectorL = (*it)->getlVector();
            double cosA = (n.x * vectorL.x) + (n.y * vectorL.y) + (n.z * vectorL.z);
            if(cosA > 0) {
                //Reflextions:
                Vector3D r  = Vector3D::normalise(2 * n * cosA - vectorL);
                Vector3D camera = Vector3D::normalise(-PointP);
                double cosB = (r.x * camera.x) + (r.y * camera.y) + (r.z * camera.z);

                if(cosB > 0){
                   // cout << specularLight.getRedDouble() << "*" <<(*it)->specularLight.getRedDouble() << "*"  << cosB << endl;
                    cosB = pow(cosB,ms);
                    redReflextionPoint = specularLight.getRedDouble() * (*it)->specularLight.getRedDouble() * cosB;
                    greenReflextionPoint = specularLight.getGreenDouble() * (*it)->specularLight.getGreenDouble() * cosB;
                    blueReflextionPoint = specularLight.getBlueDouble() * (*it)->specularLight.getBlueDouble() * cosB;
                }
            }
            newColor.setRed(newColor.getRedDouble() + redDiffusePoint + redReflextionPoint);
            newColor.setGreen(newColor.getGreenDouble() + greenDiffusePoint + greenReflextionPoint);
            newColor.setBlue(newColor.getBlueDouble() + blueDiffusePoint + blueReflextionPoint);
        }
    }

    newColor.setRed(newColor.getRedDouble() + ambientColorAnddiffuseInfColor.getRedDouble());
    newColor.setGreen(newColor.getGreenDouble() + ambientColorAnddiffuseInfColor.getGreenDouble());
    newColor.setBlue(newColor.getBlueDouble() + ambientColorAnddiffuseInfColor.getBlueDouble());

    return newColor;
}

img::Figure* img::Figure::copy() {
    Figure* copyFigure = new Figure();

    copyFigure->faces = this->faces;
    /*for(auto it = this->faces.begin(); it != this->faces.end(); it++){
		copyFigure->faces.push_back(*it);
    }*/
	//copyFigure->faces = this->faces;
    for(auto it = this->points.begin(); it != this->points.end(); it++){
        copyFigure->points.push_back(Vector3D::point(*it));
    }

    //std::copy(this->faces.begin(), this->faces.end(),copyFigure->faces.begin());
    //std::copy(this->points.begin(), this->points.end(),copyFigure->points.begin());
    copyFigure->ambientReflection = this->ambientReflection;
    copyFigure->diffuseReflection = this->diffuseReflection;
    copyFigure->reflectionCoefficient = this->reflectionCoefficient;
    copyFigure->specularReflection = this->specularReflection;
    return copyFigure;
}

