#include "RegularGrid.h"

int cint(double x){

	double fractpart, intpart;
	fractpart = modf (x , &intpart);

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

float CopySign(float x, float y)
{
	if (y >= 0)
		return fabs(x);
	else
		return fabs(x) * -1;
}

void normalize(double* vec)
{
	double a = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	if (a == 0.0) return;
	vec[0] /= a;
	vec[1] /= a;
}

void Multiply2x2(double* a, double* b, double* c)
{
	c[0] = a[0]*b[0] + a[1]*b[2];
	c[1] = a[0]*b[1] + a[1]*b[3];
	c[2] = a[2]*b[0] + a[3]*b[2];
	c[3] = a[2]*b[1] + a[3]*b[3];
}

void Invert2x2(double* a, double* c)
{
	float det = a[0]*a[3] - a[1]*a[2];
	c[0] = a[3]/det;
	c[1] = -a[1]/det;
	c[2] = -a[2]/det;
	c[3] = a[0]/det;
}

void GetTensorFromEigenSystem(double* a, double* eval, double* evec)
{
	double Q[4];
	Q[0] = evec[0];
	Q[2] = evec[1];
	Q[1] = evec[2];
	Q[3] = evec[3];
	
	double Qi[4];
	Invert2x2(Q, Qi);
	
	double M[4];
	M[0] = eval[0];
	M[3] = eval[1];
	M[1] = M[2] = 0.0;
	
	double tmp[4];
	Multiply2x2(Q, M, tmp);
	Multiply2x2(tmp, Qi, a);
}

/*void writeNrrd(void* data, const string& filename, int data_type, const vector< size_t >& dims, const vector< double >& spacing)
{
	Nrrd *nout = nrrdNew();

	if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
		//cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		//exit(-1);
	}
	if (spacing.size() == dims.size()) {
		nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
	}
	if (nrrdSave(filename.c_str(), nout, NULL)) {
		//cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		//exit(-1);
	}
}

void writeRawFile2D(float* data,const char *filename, int width, int height)
{
	vector< double > spacing;
	spacing.push_back(1.0);
	spacing.push_back(1.0);
	spacing.push_back(1.0);

	vector< size_t > dims;
	dims.push_back(2);
	dims.push_back(width);
	dims.push_back(height);

	string file_string(filename);
	writeNrrd(data, file_string, nrrdTypeFloat, dims, spacing);

	printf("Write '%s'\n", filename);
}*/

RegularGrid::RegularGrid(double* _tensors, double* _eigvals, double* _eigevecs, float4 _domain, int _width, int _height, double spcx, double spcy)
{
	width = _width;
	height = _height;
	spc.x = spcx;
	spc.y = spcy;
	eigvals = _eigvals;
	eigevecs = _eigevecs;
	tensors = _tensors;
	domain = _domain;
	
	double* tmp = (double*) malloc(4 * width * height * sizeof(double));
	for (int i = 0; i < width * height; i++)
	{
		tmp[4 * i + 0] = tensors[i + 0 * width * height];
		tmp[4 * i + 1] = tensors[i + 1 * width * height];
		tmp[4 * i + 2] = tensors[i + 2 * width * height];
		tmp[4 * i + 3] = tensors[i + 3 * width * height];
	}
	memcpy(tensors, tmp, 4 * width * height * sizeof(double));
	free(tmp);
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
	if ((spt.x < domain.x) || (spt.x > domain.y) || (spt.y < domain.z) || (spt.y > domain.w))
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
	
	gpt.x = (spt.x - domain.x) / spc.x;
	gpt.y = (spt.y - domain.z) / spc.y;

	return gpt;
}

float2 RegularGrid::Grid2Space(float2 gpt)
{
	float2 spt;

	spt.x = domain.x + gpt.x * spc.x;
	spt.y = domain.z + gpt.y * spc.y;
	
	return spt;
}

bool RegularGrid::getEigenSystem(float2 spt, float2& lambda, float2& xi1, float2& xi2)
{
	float2 gpt = Space2Grid(spt);
	int2 l = make_int2(floor(gpt.x), floor(gpt.y));
	int2 u = make_int2(ceil(gpt.x), ceil(gpt.y));
	
	if ((!IsValid(l.x, l.y)) || (!IsValid(u.x, u.y)))// || (IsBoundary(l.x, l.y)) || (IsBoundary(u.x, u.y)))
	{
		//printf("Out of boundary.\n");
		return false;
	}

	double x = gpt.x - l.x;
	double y = gpt.y - l.y;

	// get the corners
	double* V00 = &(tensors[4 * Coord2Addr(l.x, l.y)]);
	double* V01 = &(tensors[4 * Coord2Addr(l.x, u.y)]);
	double* V10 = &(tensors[4 * Coord2Addr(u.x, l.y)]);
	double* V11 = &(tensors[4 * Coord2Addr(u.x, u.y)]);

	// do the interpolation
	double ten[4];
	ten[0] =  V00[0] * (1 - x) * (1 - y) +
			  V01[0] * (1 - x) * y +
			  V10[0] * x * (1 - y) +
			  V11[0] * x * y;
	ten[1] =  V00[1] * (1 - x) * (1 - y) +
			  V01[1] * (1 - x) * y +
			  V10[1] * x * (1 - y) +
			  V11[1] * x * y;
	ten[2] =  V00[2] * (1 - x) * (1 - y) +
			  V01[2] * (1 - x) * y +
			  V10[2] * x * (1 - y) +
			  V11[2] * x * y;
	ten[3] =  V00[3] * (1 - x) * (1 - y) +
			  V01[3] * (1 - x) * y +
			  V10[3] * x * (1 - y) +
			  V11[3] * x * y;

	// find the eigensystem
	/*gsl_matrix_view m = gsl_matrix_view_array(ten, 2, 2);
	gsl_vector* evalg = gsl_vector_alloc(2);
	gsl_matrix* evecg = gsl_matrix_alloc(2, 2);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(2);
	
	gsl_eigen_symmv(&m.matrix, evalg, evecg, w);
	gsl_eigen_symmv_sort(evalg, evecg, GSL_EIGEN_SORT_VAL_ASC);
	lambda.x = gsl_vector_get(evalg, 0);
	xi1.x = gsl_matrix_get(evecg, 0, 0);
	xi1.y = gsl_matrix_get(evecg, 1, 0);
	xi1 = normalize(xi1);
	lambda.y = gsl_vector_get(evalg, 1);
	xi2.x = gsl_matrix_get(evecg, 0, 1);
	xi2.y = gsl_matrix_get(evecg, 1, 1);
	xi2 = normalize(xi2);
	
	gsl_eigen_symmv_free(w);
	gsl_vector_free(evalg);
	gsl_matrix_free(evecg);*/
	
	// find the eigensystem
	double T = ten[0] + ten[3];
	double D = ten[0] * ten[3] - ten[1] * ten[2];
	lambda.x = T/2.0 + sqrt(T*T/4 - D);
	lambda.y = T/2.0 - sqrt(T*T/4 - D);
	if (ten[2] != 0.0)
	{
		xi1 = make_float2(lambda.x - ten[3], ten[2]);
		xi2 = make_float2(lambda.y - ten[3], ten[2]);
	}
	else if (ten[1] != 0.0)
	{
		xi1 = make_float2(ten[1], lambda.x - ten[0]);
		xi2 = make_float2(ten[1], lambda.y - ten[0]);
	}
	else
	{
		xi1 = make_float2(1,0);
		xi2 = make_float2(0,1);
	}
	if (lambda.x > lambda.y)
	{
		float2 tmp;
		tmp = xi1;
		xi1 = xi2;
		xi2 = tmp;
		
		tmp = lambda;
		lambda.x = tmp.y;
		lambda.y = tmp.x;
	}
	xi1 = normalize(xi1);
	xi2 = normalize(xi2);

	// note: as in page 5 of haller paper
	// both eigen values are positive
	// and the tensor is positive definite
	
	if ((lambda.x < 0.0) || (lambda.y < 0.0))
	{
		mexPrintf("Negative eval %lf %lf\n", lambda.x, lambda.y);
		return false;
	}
	
	return true;
}

float2 RegularGrid::getVectorXi1Ref(float2 spt, bool ispos, float2& ref, bool& failed)
{
	if (!InDomain(spt))
	{
		failed = true;
		return make_float2(1.0);
	}
	
	// get the eigen system
	float2 lambda, xi1, xi2;
	bool res = getEigenSystem(spt, lambda, xi1, xi2);
	if (res == false)
	{
		failed = true;
		return make_float2(1.0);
	}
	
	// adjust xi1 to reference and set new reference
	if (dot(xi1, ref) < 0.0)
		xi1 = -xi1;
	ref = xi1;
	
	// adjust xi2 with xi1
	double cs = cos(-MPI/2.0);
	double sn = sin(-MPI/2.0);
	float2 rxi1;
	rxi1.x = xi1.x * cs - xi1.y * sn; 
	rxi1.y = xi1.x * sn + xi1.y * cs;
	if (dot(rxi1, xi2) < 0.0)
	{
		xi2 = -xi2;
	}

	// compute the eta field
	float2 vec;
	double l1 = lambda.x;
	double l2 = lambda.y;
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

float2 RegularGrid::getVectorVecRef(float2 spt, bool ispos, float2 ref, bool& failed)
{
	if (!InDomain(spt))
	{
		failed = true;
		return make_float2(1.0);
	}
	
	// get the eigen system
	float2 lambda, xi1, xi2;
	bool res = getEigenSystem(spt, lambda, xi1, xi2);
	if (res == false)
	{
		failed = true;
		return make_float2(1.0);
	}
	
	// adjust xi2 with xi1
	double cs = cos(-MPI/2.0);
	double sn = sin(-MPI/2.0);
	float2 rxi1;
	rxi1.x = xi1.x * cs - xi1.y * sn; 
	rxi1.y = xi1.x * sn + xi1.y * cs;
	if (dot(rxi1, xi2) < 0.0)
	{
		xi2 = -xi2;
	}

	// compute the eta field
	float2 vec;
	double l1 = lambda.x;
	double l2 = lambda.y;
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
	
	// adjust vec to reference and set new reference
	if (dot(vec, ref) < 0.0)
		vec = -vec;
	
	return vec;
}

float2 RegularGrid::VecRef2XiRef(float2 pt, bool ispos, float2 vecref)
{
	// get the eigen system
	float2 lambda, xi1, xi2;
	bool res = getEigenSystem(pt, lambda, xi1, xi2);
	if (res == false)
	{
		return make_float2(1.0);
	}
	
	// adjust xi2 with xi1
	double cs = cos(-MPI/2.0);
	double sn = sin(-MPI/2.0);
	float2 rxi1;
	rxi1.x = xi1.x * cs - xi1.y * sn; 
	rxi1.y = xi1.x * sn + xi1.y * cs;
	if (dot(rxi1, xi2) < 0.0)
	{
		xi2 = -xi2;
	}

	// compute the eta field
	float2 vec;
	double l1 = lambda.x;
	double l2 = lambda.y;
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
	
	// adjust vec to reference and set new reference
	if (dot(vec, vecref) < 0.0)
		xi1 = -xi1;
	
	return xi1;

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
	
	gpt.x = (spt.x - domain.x) / c_spc.x;
	gpt.y = (spt.y - domain.z) / c_spc.y;

	return gpt;
}

float2 RegularGrid::c_Grid2Space(float2 gpt)
{
	float2 spt;

	spt.x = domain.x + gpt.x * c_spc.x;
	spt.y = domain.z + gpt.y * c_spc.y;
	
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

vector<int> RegularGrid::GetCellsSharedPoints(int cell1, int cell2)
{
	vector<int> sp;
	
	int c1[4];
	c1[0] = cell1;
	c1[1] = cell1 + 1;
	c1[2] = cell1 + height;
	c1[3] = cell1 + height + 1;
	
	int c2[4];
	c2[0] = cell2;
	c2[1] = cell2 + 1;
	c2[2] = cell2 + height;
	c2[3] = cell2 + height + 1;

	for (int r = 0; r < 4; r++)
	{
		for (int s = 0; s < 4; s++)
		{
			if (c1[r] == c2[s])
				sp.push_back(c1[r]);
		}
	}
	
	return sp;
}

vector<int> RegularGrid::c_GetCellsSharedPoints(int cell1, int cell2)
{
	vector<int> sp;
	
	int c1[4];
	c1[0] = cell1;
	c1[1] = cell1 + 1;
	c1[2] = cell1 + c_height;
	c1[3] = cell1 + c_height + 1;
	
	int c2[4];
	c2[0] = cell2;
	c2[1] = cell2 + 1;
	c2[2] = cell2 + c_height;
	c2[3] = cell2 + c_height + 1;

	for (int r = 0; r < 4; r++)
	{
		for (int s = 0; s < 4; s++)
		{
			if (c1[r] == c2[s])
				sp.push_back(c1[r]);
		}
	}
	
	return sp;
}

float2 RegularGrid::Cell2Space(int cell)
{
	int2 c = Addr2Coord(cell);
	float2 gpt = make_float2(c.x, c.y);
	return Grid2Space(gpt);
}

float2 RegularGrid::c_Cell2Space(int cell)
{
	int2 c = c_Addr2Coord(cell);
	float2 gpt = make_float2(c.x, c.y);
	return c_Grid2Space(gpt);
}

void RegularGrid::WriteEtaToNrrd()
{
	float2* data = (float2*) malloc(width * height * sizeof(float2));
	
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			bool failed = false;
			float2 onef2 = make_float2(1.0);
			float2 vec = getVectorXi1Ref(Grid2Space(make_float2(x, y)), true, onef2, failed);
			data[x + width * y] = vec;
		}
	}
	
	// write the file
	FILE * pFile;
	pFile = fopen ("E:\\Chaos2Rost\\Projects\\Tensor2DTopology\\doubleGyre.nrrd","w");
	fprintf (pFile, "NRRD0001\n");
	fprintf (pFile, "# Complete NRRD file format specification at:\n");
	fprintf (pFile, "# http://teem.sourceforge.net/nrrd/format.html\n");
	fprintf (pFile, "type: float\n");
	fprintf (pFile, "dimension: 3\n");
	fprintf (pFile, "sizes: 2 %d %d\n", width, height);
	fprintf (pFile, "spacings: nan %f %f\n", spc.x, spc.y);
	fprintf (pFile, "endian: little\n");
	fprintf (pFile, "encoding: ascii\n");
	for (int i = 0; i < width * height; i++)
	{
		if ((i % 4) == 0)
			fprintf(pFile, "\n");
		fprintf (pFile, "%f %f ", data[i].x, data[i].y);
	}
	fclose (pFile);
	
	// write the file
	pFile = fopen ("E:\\Chaos2Rost\\Projects\\Tensor2DTopology\\doubleGyreTen.nrrd","w");
	fprintf (pFile, "NRRD0001\n");
	fprintf (pFile, "# Complete NRRD file format specification at:\n");
	fprintf (pFile, "# http://teem.sourceforge.net/nrrd/format.html\n");
	fprintf (pFile, "type: float\n");
	fprintf (pFile, "dimension: 3\n");
	fprintf (pFile, "sizes: 4 %d %d\n", width, height);
	fprintf (pFile, "spacings: nan %f %f\n", spc.x, spc.y);
	fprintf (pFile, "endian: little\n");
	fprintf (pFile, "encoding: ascii\n");
	for (int i = 0; i < width * height; i++)
	{
		if ((i % 2) == 0)
			fprintf(pFile, "\n");
		fprintf (pFile, "%lf %lf %lf %lf ", tensors[4 * i + 0], tensors[4 * i + 1], tensors[4 * i + 2], tensors[4 * i + 3]);
	}
	fclose (pFile);
}

void RegularGrid::LoadTensorFromNrrd()
{
	// write the file
	FILE * pFile;
	char buf[1000];
	pFile = fopen ("E:\\Chaos2Rost\\Projects\\Tensor2DTopology\\doubleGyreTen.nrrd","w");
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	fgets(buf, 1000, pFile); 
	for (int i = 0; i < width * height; i++)
	{
		if ((i % 2) == 0)
			fgets(buf, 1000, pFile);
		fscanf (pFile, "%lf %lf %lf %lf ", &(tensors[4 * i + 0]), &(tensors[4 * i + 1]), &(tensors[4 * i + 2]), &(tensors[4 * i + 3]));
	}
	fclose (pFile);
}

















double2 GetDoubleGyreFlowV(double2 spt, double t)
{
	// The Double Gyre Dataset
	/*double epsilon = 0.1;
	double a = 0.1;
    double omega = 2.0 * MPI / 10.0;
	double esot = epsilon * sin(omega * t);
	double forcing = esot * spt.x * spt.x + (1.0 - 2.0 * esot) * spt.x;
	double2 v;
	v.x = -MPI * a * sin(MPI * forcing) * cos(MPI * spt.y);
	v.y = MPI * a * cos(MPI * forcing) * sin(MPI * spt.y) * (2.0 * esot * spt.x + 1.0 - 2.0 * esot);*/
	
	// The traveling Wave Dataset
	/*double omega = 1;
    double amplitude = 1;
    double waveNumber = 1;
    double speed = .5;
    double forcingAmplitude = .25;
	double forcing = 0.0;//sin(omega*t);
	double2 v;
	v.x = speed - amplitude*sin(waveNumber*spt.x)*cos(spt.y) - forcingAmplitude*forcing;
	v.y = amplitude*waveNumber*cos(waveNumber*spt.x)*sin(spt.y);*/

	//mexPrintf("%lf %lf\n", v.x, v.y);
	
	// calling matlab for the flow function
	mxArray *fd_prhs[3];
	fd_prhs[0] = dataset_flow;
	fd_prhs[1] = mxCreateDoubleScalar(t);
	fd_prhs[2] = mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);
	*(mxGetPr(fd_prhs[2]) + 0) = spt.x;
	*(mxGetPr(fd_prhs[2]) + 1) = spt.y;
	mxArray *fd_plhs[1];
	fd_plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
	int result = mexCallMATLAB(1, fd_plhs, 3, fd_prhs, "flow_derivative_func");
	double2 v;
	v.x = *(mxGetPr(fd_plhs[0]) + 0);
	v.y = *(mxGetPr(fd_plhs[0]) + 1);
	
	//printf("%d %lf %lf\n", result, v.x, v.y);
	
	return v;
}

int funcGetDoubleGyreFlow(double t, const double y[], double f[],  void *params)
{
	double2 yp = make_double2(y[0], y[1]);
	double2 v = GetDoubleGyreFlowV(yp, t);
	f[0] = v.x;
	f[1] = v.y;
	
	return GSL_SUCCESS;
}

void IntegrateDoubleGyreFlowGSL(double2& y, float st, float et, double mu)
{
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);
	gsl_odeiv_system sys = {funcGetDoubleGyreFlow, NULL, 2, NULL};

	double t = st, t1 = et;
	double h = mu;
	double ydp[2] = { y.x, y.y };
	
	//gsl_ieee_env_setup();
	while (t < t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		if (status != GSL_SUCCESS)
		{
			break;
		}
	}
	y = make_double2(ydp[0], ydp[1]);

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}

void RegularGrid::GetDoubleGyreCGTensor()
{
	omp_set_num_threads(16);

	mexPrintf("Start tensor computation\n");
	
	// compute the flow map
	double2* flowmap = (double2*) malloc(width * height * sizeof(double2)); 
	/*#pragma omp parallel for schedule(dynamic,25)
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			float2 spt = Grid2Space(i, j);
			flowmap[Coord2Addr(i,j)] = make_double2(spt.x, spt.y);
			IntegrateDoubleGyreFlowGSL(&(flowmap[Coord2Addr(i,j)]), 0.0, 20.0, 1e-9);
		}
	}*/
	
	// compute the derivative of the flow map
	double4* derflowmap = (double4*) malloc(width * height * sizeof(double4)); 
	double msp = 0.25 * min(spc.x, spc.y);
	for (int i = 0; i < width; i++)
	{
		if ((i % (width / 10)) == 0)
		{
			mexPrintf("%.1f%% ", 100.0 * i / float(width - 1)); 
			mexEvalString("drawnow;");
		}
		
		double lgnt = 20.0;
		//#pragma omp parallel for schedule(dynamic,25)
		for (int j = 0; j < height; j++)
		{
			float2 spt = Grid2Space(make_float2(i, j));
			double2 ptx0 = make_double2(spt.x-msp, spt.y);
			double2 ptx1 = make_double2(spt.x+msp, spt.y);
			double2 pty0 = make_double2(spt.x, spt.y-msp);
			double2 pty1 = make_double2(spt.x, spt.y+msp);
			
			IntegrateDoubleGyreFlowGSL(ptx0, 0.0, lgnt, 1e-10);
			IntegrateDoubleGyreFlowGSL(ptx1, 0.0, lgnt, 1e-10);
			IntegrateDoubleGyreFlowGSL(pty0, 0.0, lgnt, 1e-10);
			IntegrateDoubleGyreFlowGSL(pty1, 0.0, lgnt, 1e-10);
			
			int idx = Coord2Addr(i, j);
			derflowmap[idx].x = 0.5 * (ptx1.x - ptx0.x) / msp;
			derflowmap[idx].y = 0.5 * (pty1.x - pty0.x) / msp;
			derflowmap[idx].z = 0.5 * (ptx1.y - ptx0.y) / msp;
			derflowmap[idx].w = 0.5 * (pty1.y - pty0.y) / msp;
			
			//mexPrintf("%lf %lf %lf %lf\n", derflowmap[idx].x, derflowmap[idx].y, derflowmap[idx].z, derflowmap[idx].w);
		}

	}
	mexPrintf("100%%\n");
	
	// compute the tensor field
	for (int i = 0; i < width; i++)
	{
		if ((i % (width / 10)) == 0)
		{
			mexPrintf("%.1f%% ", 100.0 * i / float(width - 1)); 
			mexEvalString("drawnow;");
		}
		
		#pragma omp parallel for schedule(dynamic,25)
		for (int j = 0; j < height; j++)
		{
			int idx = Coord2Addr(i, j);
			tensors[4 * idx + 0] = derflowmap[idx].x * derflowmap[idx].x + derflowmap[idx].z * derflowmap[idx].z;
			tensors[4 * idx + 1] = derflowmap[idx].x * derflowmap[idx].y + derflowmap[idx].z * derflowmap[idx].w;
			tensors[4 * idx + 2] = derflowmap[idx].y * derflowmap[idx].x + derflowmap[idx].w * derflowmap[idx].z;
			tensors[4 * idx + 3] = derflowmap[idx].y * derflowmap[idx].y + derflowmap[idx].w * derflowmap[idx].w;
			//mexPrintf("%lf %lf %lf \n", tensors[4 * idx], tensors[4 * idx + 1], tensors[4 * idx + 2]);
		}
	}
	mexPrintf("100%%\n");
	
	free(flowmap);
	free(derflowmap);
	mexPrintf("End tensor computation\n\n\n");
	mexEvalString("drawnow;");
}