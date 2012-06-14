#include <stdio.h>
#include <mex.h> 
#include <math.h>
#include <vector>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;

int cint(double x);
void normalize(double* vec);

class RegularGrid
{
public:
	RegularGrid(double* _eigvals, double* _eigevecs, int _width, int _height, double spcx, double spcy);
	~RegularGrid(void);


	////// Properties
	int width;
	int height;
	float2 spc;
	double* eigvals;
	double* eigevecs;

	////// Methods
	bool IsBoundary(const int& x, const int& y);
	bool IsValid(const int& x, const int& y);
	bool InDomain(float2 spt);
	int Coord2Addr(const int& x, const int& y);
	int Coord2Addr(const int2& coord);
	int2 Addr2Coord(const int& i);
	float2 Space2Grid(float2 spt);
	float2 Grid2Space(float2 gpt);
	
	float2 getVectorAtGridPoint(int2 pt, bool ispos);
	float2 getVector(float2 spt, bool ispos, float2 ref, bool& failed);
	
	
	int c_width;
	int c_height;
	float2 c_spc;
	void c_SetCoarseRes(int _c_width, int _c_height);
	void c_SetCoarseRes(double msx, double msy);
	int c_Coord2Addr(const int& x, const int& y);
	int c_Coord2Addr(const int2& coord);
	int2 c_Addr2Coord(const int& i);
	float2 c_Space2Grid(float2 spt);
	float2 c_Grid2Space(float2 gpt);
	
	vector<int> GetCellsSharedPoints(int cell1, int cell2);
	vector<int> c_GetCellsSharedPoints(int cell1, int cell2);
};

