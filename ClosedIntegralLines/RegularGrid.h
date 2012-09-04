#include <stdio.h>
#include <mex.h> 
#include <math.h>
#include <string.h>
#include <vector>
#include <set>
#include <queue>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <omp.h>

#include "vectortypes.h"
//#include <teem/nrrd.h>

#define MPI 3.141592653589793238462643383279

using namespace std;

int cint(double x);
void normalize(double* vec);
float CopySign(float x, float y);

class RegularGrid
{
public:
	RegularGrid(double* _tensors, double* _eigvals, double* _eigevecs, float4 _domain, int _width, int _height, double spcx, double spcy);
	~RegularGrid(void);


	////// Properties
	int width;
	int height;
	float2 spc;
	double* eigvals;
	double* eigevecs;
	double* tensors;
	float4 domain;
	
	////// Methods
	bool IsBoundary(const int& x, const int& y);
	bool IsValid(const int& x, const int& y);
	bool InDomain(float2 spt);
	int Coord2Addr(const int& x, const int& y);
	int Coord2Addr(const int2& coord);
	int2 Addr2Coord(const int& i);
	float2 Space2Grid(float2 spt);
	float2 Grid2Space(float2 gpt);
	float2 Cell2Space(int cell);
	
	bool getEigenSystem(float2 spt, float2& lambda, float2& xi1, float2& xi2);
	float2 getVectorXi1Ref(float2 spt, bool ispos, float2& ref, bool& failed);
	float2 getVectorVecRef(float2 spt, bool ispos, float2 ref, bool& failed);
	float2 VecRef2XiRef(float2 pt, bool ispos, float2 vecref);
	
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
	float2 c_Cell2Space(int cell);
	
	vector<int> GetCellsSharedPoints(int cell1, int cell2);
	vector<int> c_GetCellsSharedPoints(int cell1, int cell2);
	
	void WriteEtaToNrrd();
	void LoadTensorFromNrrd();
	
	
	void GetDoubleGyreCGTensor();
};

