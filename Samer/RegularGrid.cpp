#include "RegularGrid.h"

int cint(double x){

	double fractpart, intpart;
	fractpart = modf (x , &intpart);

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

void normalize(double* vec)
{
	double a = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	if (a == 0.0) return;
	vec[0] /= a;
	vec[1] /= a;
}

RegularGrid::RegularGrid(double* _eigvals, double* _eigevecs, int _width, int _height, double spcx, double spcy)
{
	width = _width;
	height = _height;
	spc.x = spcx;
	spc.y = spcy;
	eigvals = _eigvals;
	eigevecs = _eigevecs;
}

bool RegularGrid::IsBoundary(const int& x, const int& y)
{
	if ((x < 1) || (y < 1) || (x > width - 2) || (y > height - 2))
		return true;
	else return false;
}

bool RegularGrid::IsValid(const int& x, const int& y)
{
	if ((x < 0) || (y < 0) || (x > width - 1) || (y > height - 1))
		return false;
	else return true;
}

bool RegularGrid::InDomain(float2 spt)
{
	if ((spt.x < 0.0) || (spt.x > spc.x * (width - 1)) || (spt.y < 0.0) || (spt.y > spc.y * (height - 1)))
		return false;
	return true;
}

int RegularGrid::Coord2Addr(const int& x, const int& y)
{
	return y + height * x;
}

int RegularGrid::Coord2Addr(const int2& coord)
{
	return Coord2Addr(coord.x, coord.y);
}

int2 RegularGrid::Addr2Coord(const int& i)
{
	int y = i % height;
	int x = ((i - y) / height);

	return make_int2(x, y);
}

float2 RegularGrid::Space2Grid(float2 spt)
{
	float2 gpt;
	
	gpt.x = spt.x / spc.x;
	gpt.y = spt.y / spc.y;

	return gpt;
}

float2 RegularGrid::Grid2Space(float2 gpt)
{
	float2 spt;

	spt.x = gpt.x * spc.x;
	spt.y = gpt.y * spc.y;
	
	return spt;
}

float2 RegularGrid::getVectorAtGridPoint(int2 pt, bool ispos)
{
	int idx = Coord2Addr(pt);
	
	double l1 = eigvals[idx + 0 * width * height];
    double l2 = eigvals[idx + 1 * width * height];
    float2 xi1 = make_float2(eigevecs[idx + 0 * width * height], eigevecs[idx + 1 * width * height]);
    float2 xi2 = make_float2(eigevecs[idx + 2 * width * height], eigevecs[idx + 3 * width * height]);
	
	float2 vec;
	
	if (ispos)
	{
		vec.x = sqrt(sqrt(l2) / (sqrt(l1) + sqrt(l2))) * xi1.x + sqrt(sqrt(l1) / (sqrt(l1) + sqrt(l2))) * xi2.x;
        vec.y = sqrt(sqrt(l2) / (sqrt(l1) + sqrt(l2))) * xi1.y + sqrt(sqrt(l1) / (sqrt(l1) + sqrt(l2))) * xi2.y;
	}
	else
	{
		vec.x = sqrt(sqrt(l2) / (sqrt(l1) + sqrt(l2))) * xi1.x - sqrt(sqrt(l1) / (sqrt(l1) + sqrt(l2))) * xi2.x;
        vec.y = sqrt(sqrt(l2) / (sqrt(l1) + sqrt(l2))) * xi1.y - sqrt(sqrt(l1) / (sqrt(l1) + sqrt(l2))) * xi2.y;
	}
	
	return vec;
}

float2 RegularGrid::getVector(float2 spt, bool ispos, float2 ref, bool& failed)
{
	if (!InDomain(spt))
	{
		failed = true;
		return make_float2(1.0);
	}
	
	float2 gpt = Space2Grid(spt);
	
	// get corners
	float2 vec1 = getVectorAtGridPoint(make_int2(floor(gpt.x), floor(gpt.y)), ispos);
	float2 vec2 = getVectorAtGridPoint(make_int2(floor(gpt.x), ceil(gpt.y)), ispos);
	float2 vec3 = getVectorAtGridPoint(make_int2(ceil(gpt.x), floor(gpt.y)), ispos);
	float2 vec4 = getVectorAtGridPoint(make_int2(ceil(gpt.x), ceil(gpt.y)), ispos);
	
	// correct directions
	if (dot(vec1, vec2) < 0.0)
		vec2 = -vec2;
	if (dot(vec1, vec3) < 0.0)
		vec3 = -vec3;
	if (dot(vec1, vec4) < 0.0)
		vec4 = -vec4;
		
	// now interpolation
	float x = gpt.x - floor(gpt.x);
	float y = gpt.y - floor(gpt.y);
	
	float2 vec = vec1 * (1 - x) * (1 - y)
			   + vec2 * (1 - x) * y
			   + vec3 * x * (1 - y)
			   + vec4 * x * y;
			   
	if (dot(vec, ref) < 0)
		vec = -vec;
		
	return vec;
}

int RegularGrid::c_Coord2Addr(const int& x, const int& y)
{
	return y + c_height * x;
}

int RegularGrid::c_Coord2Addr(const int2& coord)
{
	return c_Coord2Addr(coord.x, coord.y);
}

int2 RegularGrid::c_Addr2Coord(const int& i)
{
	int y = i % c_height;
	int x = ((i - y) / c_height);

	return make_int2(x, y);
}

float2 RegularGrid::c_Space2Grid(float2 spt)
{
	float2 gpt;
	
	gpt.x = spt.x / c_spc.x;
	gpt.y = spt.y / c_spc.y;

	return gpt;
}

float2 RegularGrid::c_Grid2Space(float2 gpt)
{
	float2 spt;

	spt.x = gpt.x * c_spc.x;
	spt.y = gpt.y * c_spc.y;
	
	return spt;
}

void RegularGrid::c_SetCoarseRes(int _c_width, int _c_height)
{
	c_width = _c_width;
	c_height = _c_height;
	c_spc.x = spc.x * (width - 1) / (c_width - 1);
	c_spc.y = spc.y * (height - 1) / (c_height - 1);
}

void RegularGrid::c_SetCoarseRes(double msx, double msy)
{
	c_spc.x = spc.x * msx;
	c_spc.y = spc.y * msy;
	c_width = ceil(spc.x * (width - 1) / c_spc.x) + 1;
	c_height = ceil(spc.y * (height - 1) / c_spc.y) + 1;
}