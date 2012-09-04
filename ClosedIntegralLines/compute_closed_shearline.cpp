/* $Revision: 1.8.6.4 $ */
/*=========================================================
 * convec.c
 * example for illustrating how to use pass complex data 
 * from MATLAB to C and back again
 *
 * convolves  two complex input vectors
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 *=======================================================*/
#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <set>
#include <map>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
	 
#include "RegularGrid.h"
#include <omp.h>

using namespace std;

RegularGrid* grid;

int cycleiterc = 2; // not less than 3
double intg_lng = 40; // maximum integration length for detection of cycles, WARNING: Overwritten in the code based on dimensions
double error_s = 1e-7; // error use for the integrtion
double error_b = 1e-7; // same as previous one
double tole = 0.0;//0.004; // tolerance distance between end poitns of a cycle to ignore it as proof of exit


float2 GetFlowAt(float time, float2 spt, float2& ref, bool& failed, bool ispos)
{
    return grid->getVectorXi1Ref(spt, ispos, ref, failed);
}

bool CycleFound(vector<int>& cells)
{
	// find occurances
	int c2 = cells.size() - 1;
	int c0, c1;
	int i;
	
	// find c1
	i = c2 - 1;
	while (cells[i] == cells[c2])
		i--;
	while (cells[i] != cells[c2])
		i--;
	c1 = i;
	
	// find c0
	i = c1 - 1;
	while (cells[i] == cells[c2])
		i--;
	while (cells[i] != cells[c2])
		i--;
	c0 = i;
	
	// find first sequence of cells
	vector<int> l1;
	for (i = c0; i < c1; i++)
	{
		if ((l1.empty()) || (l1.back() != cells[i]))
			l1.push_back(cells[i]);
	}
	
	// find second sequence of cells
	vector<int> l2;
	for (i = c1; i < c2; i++)
	{
		if ((l2.empty()) || (l2.back() != cells[i]))
			l2.push_back(cells[i]);
	}
	
	if (l1.size() != l2.size())
	{
		return false;
	}
	
	for (i = 0; i < l1.size(); i++)
	{
		if (l1[i] != l2[i])
			return false;
	}
	
	return true;
}

int inSegment( float2 P, float2 SP0, float2 SP1)
{
    if (SP0.x != SP1.x) {    // S is not vertical
        if (SP0.x <= P.x && P.x <= SP1.x)
            return 1;
        if (SP0.x >= P.x && P.x >= SP1.x)
            return 1;
    }
    else {    // S is vertical, so test y coordinate
        if (SP0.y <= P.y && P.y <= SP1.y)
            return 1;
        if (SP0.y >= P.y && P.y >= SP1.y)
            return 1;
    }
    return 0;
}

#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product (2D)
int intersect2D_Segments( float2 S1P0, float2 S1P1, float2 S2P0, float2 S2P1, float2& I0, float2& I1 )
{
	//mexPrintf("Check intersection\n");
	//mexPrintf("%f %f\n", S1P0.x, S1P0.y);
	//mexPrintf("%f %f\n", S1P1.x, S1P1.y);
	//mexPrintf("%f %f\n", S2P0.x, S2P0.y);
	//mexPrintf("%f %f\n", S2P1.x, S2P1.y);
	
	
    float2    u = S1P1 - S1P0;
    float2    v = S2P1 - S2P0;
    float2    w = S1P0 - S2P0;
    float     D = perp(u,v);

    // test if they are parallel (includes either being a point)
    if (fabs(D) < 0.0) {          // S1 and S2 are parallel
        if (perp(u,w) != 0 || perp(v,w) != 0) {
            return 0;                   // they are NOT collinear
        }
        // they are collinear or degenerate
        // check if they are degenerate points
        float du = dot(u,u);
        float dv = dot(v,v);
        if (du==0 && dv==0) {           // both segments are points
            if ((S1P0.x != S2P0.x) || (S1P0.y != S2P0.y))         // they are distinct points
                return 0;
            I0 = S1P0;                // they are the same point
            return 1;
        }
        if (du==0) {                    // S1 is a single point
            if (inSegment(S1P0, S2P0, S2P1) == 0)  // but is not in S2
                return 0;
            I0 = S1P0;
            return 1;
        }
        if (dv==0) {                    // S2 a single point
            if (inSegment(S2P0, S1P0, S1P1) == 0)  // but is not in S1
                return 0;
            I0 = S2P0;
            return 1;
        }
        // they are collinear segments - get overlap (or not)
        float t0, t1;                   // endpoints of S1 in eqn for S2
        float2 w2 = S1P1 - S2P0;
        if (v.x != 0) {
                t0 = w.x / v.x;
                t1 = w2.x / v.x;
        }
        else {
                t0 = w.y / v.y;
                t1 = w2.y / v.y;
        }
        if (t0 > t1) {                  // must have t0 smaller than t1
                float t=t0; t0=t1; t1=t;    // swap if not
        }
        if (t0 > 1 || t1 < 0) {
            return 0;     // NO overlap
        }
        t0 = t0<0? 0 : t0;              // clip to min 0
        t1 = t1>1? 1 : t1;              // clip to max 1
        if (t0 == t1) {                 // intersect is a point
            I0 = S2P0 + t0 * v;
            return 1;
        }

        // they overlap in a valid subsegment
        I0 = S2P0 + t0 * v;
        I1 = S2P0 + t1 * v;
        return 2;
    }

    // the segments are skew and may intersect in a point
    // get the intersect parameter for S1
    float     sI = perp(v,w) / D;
    if (sI < 0 || sI > 1)               // no intersect with S1
        return 0;

    // get the intersect parameter for S2
    float     tI = perp(u,w) / D;
    if (tI < 0 || tI > 1)               // no intersect with S2
        return 0;

    I0 = S1P0 + sI * u;               // compute S1 intersect point
    return 1;
}

struct paramsDP{
	bool ispos;
	bool isforward;
	float2 ref;
};

int funcGetFlow(double t, const double y[], double f[],  void *params)
{
	paramsDP* p = (paramsDP *)params;
	float2 yp = make_float2(y[0], y[1]);
	bool failed = false;
	float2 v = GetFlowAt(t, yp, p->ref, failed, p->ispos);
	
	// output
	if (!p->isforward) v = -v;
	f[0] = v.x;
	f[1] = v.y;
		
	if (failed) 
	{
		return !GSL_SUCCESS;
	}
	return GSL_SUCCESS;
}

void IntegrateDP(bool coarse, float2* y, float st, float et, double mu, bool ispos, bool isforward, vector<float2>& shearline, vector<int>& cells, int& checks, int& results, float2 ref)
{
	ref = grid->VecRef2XiRef(*y, ispos, ref);
	
	int cell = -1, ocell = -1;
	
	int size = coarse? grid->c_width * grid->c_height : grid->width * grid->height;
	uchar* visited = (uchar*) malloc(size * sizeof(uchar));
	memset(visited, 0, size * sizeof(uchar));
	
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);

	paramsDP params = {ispos, isforward, ref};
	gsl_odeiv_system sys = {funcGetFlow, NULL, 2, &params};

	double t = st, t1 = et;
	double h = mu;
	double ydp[2] = { (*y).x, (*y).y };
	float2 oy;
	
	//mexPrintf("\n");
	while (t < t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		oy = *y;
		*y = make_float2(ydp[0], ydp[1]);
		
		if (status != GSL_SUCCESS)
		   break;
		   
		// add the current point
		shearline.push_back(*y);
		
		// get the current point cell
		float2 gpt = coarse? grid->c_Space2Grid(*y) : grid->Space2Grid(*y);
		int2 igpt = make_int2(floor(gpt.x), floor(gpt.y));
		ocell = cell;
		cell = coarse? grid->c_Coord2Addr(igpt) : grid->Coord2Addr(igpt);
		cells.push_back(cell);
		if (cell != ocell)
		{
			visited[cell] ++;
			//mexPrintf("%d ", cell);
		}
		
		// now do the checks:
		if (checks & 1)
		{
			if (visited[cell] > cycleiterc)
			{
				if (CycleFound(cells)) // second cycle detected
				{
					results |= 1;
					break;
				}
			}
		}
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	free(visited);
}

void IntegrateDP2(bool coarse, float2* y, float st, float et, double mu, bool ispos, bool isforward, vector<float2>& shearline, vector<int>& cells, int& checks, int& results, float2 ref, float2 e1, float2 e2)
{
	ref = grid->VecRef2XiRef(*y, ispos, ref);
	
	int cell = -1, ocell = -1;
	
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);

	paramsDP params = {ispos, isforward, ref};
	gsl_odeiv_system sys = {funcGetFlow, NULL, 2, &params};

	double t = st, t1 = et;
	double h = mu;
	double ydp[2] = { (*y).x, (*y).y };
	float2 oy;
	
	int itr = 0;
	while (t < t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		oy = *y;
		*y = make_float2(ydp[0], ydp[1]);
		
		if (status != GSL_SUCCESS)
		{
			mexPrintf("Case 0\n");
			break;
		}
		   
		// add the current point
		shearline.push_back(*y);
		
		// get the current point cell
		bool ec = false;
		float2 gpt = coarse? grid->c_Space2Grid(*y) : grid->Space2Grid(*y);
		int2 igpt = make_int2(floor(gpt.x), floor(gpt.y));
		ocell = cell;
		cell = coarse? grid->c_Coord2Addr(igpt) : grid->Coord2Addr(igpt);
		cells.push_back(cell);
		
		// now do the checks:
		float2 dm1, dm2;
		if ((itr >= 20) && 
			//(visited[cell] >= 2) &&
			(ocell != cell) && 
			(intersect2D_Segments(oy, *y, e1, e2, dm1, dm2) == 1))
		{
			results |= 2;
			//mexPrintf("== edge %f %f %f %f\n", e1.x, e1.y, e2.x, e2.y);
			//mexPrintf("== cells %d %d\n", ocell, cell);
			//mexPrintf("== itr %d\n", itr);
			//mexPrintf("== time %f\n", t);
			//mexPrintf("== dm1 %f %f %d\n", dm1.x, dm1.y, ec? 1:0);
			//mexPrintf("== itr %d, %f %f %f %f\n", itr, oy.x, oy.y, (*y).x, (*y).y);
			//mexPrintf("Case 1\n");
			break;
		}
		itr++;
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}

void IntegrateDP3(bool coarse, float2* y, float st, float et, double mu, bool ispos, bool isforward, vector<float2>& shearline, vector<int>& cells, int& checks, int& results, float2 ref, set<int>& icells)
{
	ref = grid->VecRef2XiRef(*y, ispos, ref);
	
	int cell = -1, ocell = -1, firstcell = -1;
	
	int size = coarse? grid->c_width * grid->c_height : grid->width * grid->height;
	uchar* visited = (uchar*) malloc(size * sizeof(uchar));
	memset(visited, 0, size * sizeof(uchar));
	
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);

	paramsDP params = {ispos, isforward, ref};
	gsl_odeiv_system sys = {funcGetFlow, NULL, 2, &params};

	double t = st, t1 = et;
	double h = mu;
	double ydp[2] = { (*y).x, (*y).y };
	float2 oy = *y;
	
	int itr = 0;
	int outs = 0;
	double alpha = 1e-6;
	results = 1;
	double oedist, edist;
	oedist = edist = std::numeric_limits<double>::max();
	//mexPrintf("start\n");
	while (t < t1)
	{   
		// add the current point
		shearline.push_back(*y);
		
		// get the current point cell
		ocell = cell;
		bool consideredin = false;
		for (double aptx = (*y).x - alpha; aptx <= (*y).x + alpha; aptx += alpha)
		{
			for (double apty = (*y).y - alpha; apty <= (*y).y + alpha; apty += alpha)
			{
				float2 apt = make_float2(aptx, apty);
				float2 agpt = coarse? grid->c_Space2Grid(apt) : grid->Space2Grid(apt);
				int2 aigpt = make_int2(floor(agpt.x), floor(agpt.y));
				int acell = coarse? grid->c_Coord2Addr(aigpt) : grid->Coord2Addr(aigpt);
				
				if (icells.find(acell) != icells.end())
				{
					consideredin = true;
					cell = acell;
					break;
				}
			}
			if (consideredin == true)
				break;
		}
		cells.push_back(cell);
		if (ocell != cell)
			visited[cell]++;
		
		// find the first cell in the set
		if ((itr > 3) && (firstcell == -1) && (icells.find(cell) != icells.end()))
		{
			firstcell = cell;
			//mexPrintf(" -> first cell is %d\n", firstcell);
		}
		
		// check if another cell was crossed
		if (ocell != -1)
		{
			vector<int> spt = coarse? grid->c_GetCellsSharedPoints(ocell, cell) : grid->GetCellsSharedPoints(ocell, cell);
			if (spt.size() == 1)
			{
				int2 b1 = coarse? grid->c_Addr2Coord(ocell) : grid->Addr2Coord(ocell);
				int2 b2 = coarse? grid->c_Addr2Coord(cell) : grid->Addr2Coord(cell);
				int cell3 = coarse? grid->c_Coord2Addr(b1.x, b2.y) : grid->Coord2Addr(b1.x, b2.y);
				int cell4 = coarse? grid->c_Coord2Addr(b2.x, b1.y) : grid->Coord2Addr(b2.x, b1.y);				
				float2 dm1, dm2;
				
				// check cell3
				vector<int> spt3 = coarse? grid->c_GetCellsSharedPoints(cell, cell3) : grid->GetCellsSharedPoints(cell, cell3);
				float2 c30 = coarse? grid->c_Cell2Space(spt3[0]) : grid->Cell2Space(spt3[0]);
				float2 c31 = coarse? grid->c_Cell2Space(spt3[1]) : grid->Cell2Space(spt3[1]);
				if (((intersect2D_Segments(*y, oy, c30, c31, dm1, dm2) != 0) && (icells.find(cell3) == icells.end()))
					&& (length(dm1 - c30) > alpha)
					&& (length(dm1 - c31) > alpha))
				{
					results = 1;
					//mexPrintf("Case 0\n");
					break;
				}
				
				// check cell4
				vector<int> spt4 = coarse? grid->c_GetCellsSharedPoints(cell, cell4) : grid->GetCellsSharedPoints(cell, cell4);
				float2 c40 = coarse? grid->c_Cell2Space(spt4[0]) : grid->Cell2Space(spt4[0]);
				float2 c41 = coarse? grid->c_Cell2Space(spt4[1]) : grid->Cell2Space(spt4[1]);
				if (((intersect2D_Segments(*y, oy, c40, c41, dm1, dm2) != 0) && (icells.find(cell4) == icells.end()))
					&& (length(dm1 - c40) > alpha)
					&& (length(dm1 - c41) > alpha))
				{
					results = 1;
					//mexPrintf("Case 1\n");
					break;
				}	
			}
		}
		
		// if not in exit
		if (consideredin == false)
		{
			results = 1;
			//mexPrintf("Case 2\n");
			break;
		}
		
		// if the first cell is hit again exit
		if ((firstcell != -1) && (visited[firstcell] >= 2))
		{
			oedist = edist;
			edist = length(shearline.back() - shearline[0]);
			//mexPrintf("%lf\n", edist);
			if (edist >= oedist)
			{
				if (edist > tole)
					results = 0;
				else
					results = 1;
				break;
			}
		}
		
		// next point
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		oy = *y;
		*y = make_float2(ydp[0], ydp[1]);
		
		if (status != GSL_SUCCESS)
		{
			results = 1;
			//mexPrintf("Case 4\n");
			break;
		} 
		itr++;
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	free(visited);
}

void IntegrateDP4(bool coarse, float2* y, float st, float et, double mu, bool ispos, bool isforward, vector<float2>& shearline, vector<int>& cells, int& checks, int& results, float2 ref)
{
	ref = grid->VecRef2XiRef(*y, ispos, ref);
	
	int cell = -1, ocell = -1;
	
	int size = coarse? grid->c_width * grid->c_height : grid->width * grid->height;
	
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);

	paramsDP params = {ispos, isforward, ref};
	gsl_odeiv_system sys = {funcGetFlow, NULL, 2, &params};

	double t = st, t1 = et;
	double h = mu;
	double ydp[2] = { (*y).x, (*y).y };
	float2 oy;
	
	int diffcount = 0;
	while (t < t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		oy = *y;
		*y = make_float2(ydp[0], ydp[1]);
		
		if (status != GSL_SUCCESS)
		   break;
		   
		// add the current point
		shearline.push_back(*y);
		
		// get the current point cell
		float2 gpt = coarse? grid->c_Space2Grid(*y) : grid->Space2Grid(*y);
		int2 igpt = make_int2(floor(gpt.x), floor(gpt.y));
		ocell = cell;
		cell = coarse? grid->c_Coord2Addr(igpt) : grid->Coord2Addr(igpt);
		cells.push_back(cell);
		if (cell != ocell)
			diffcount++;
			
		// now do the checks:
		if ((checks & 1) && (diffcount > 3))
		{
			if (length(shearline.back() - shearline[0]) < 3 * min(grid->spc.x, grid->spc.y))
			{
				results |= 1;
				break;
			}
		}
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}


struct flt2_compare {
    bool operator() (const float2& lhs, const float2& rhs) const{
        if (lhs.x < rhs.x)
			return true;
		else if (lhs.x > rhs.x)
			return false;
		else
			return lhs.y < rhs.y;
    }
};

struct int2_compare {
    bool operator() (const int2& lhs, const int2& rhs) const{
        if (lhs.x < rhs.x)
			return true;
		else if (lhs.x > rhs.x)
			return false;
		else
			return lhs.y < rhs.y;
    }
};

float2 GetClosestDirection(vector<float2>& orgline, float2 pt)
{
	int closest;
	double mindist = std::numeric_limits<double>::max();
	for (int i = 0; i < orgline.size(); i++)
	{
		if (length(pt - orgline[i]) < mindist)
		{
			closest = i;
			mindist = length(pt - orgline[i]);
		}
	}
	
	int to = (closest + 1) % orgline.size();
	
	return orgline[to] - orgline[closest];
}

bool AllSetCovered(vector<int>& cells, set<int>& cellset)
{
	//return true;
	set<int> newset;
	for (int i = 0; i < cells.size(); i++)
		newset.insert(cells[i]);
	
	//mexPrintf("Test %d %d\n", newset.size(), cellset.size());
	
	if (newset.size() < 0.5 * cellset.size())
	{
		return false;
	}
	
	return true;
}

bool IsExitEdge(int e1, int e2, vector<float2>& orgline, bool ispos, set<int>& cellset, vector< vector<float2> >& boundaries)
{
	int2 crd1 = grid->c_Addr2Coord(e1);
	int2 crd2 = grid->c_Addr2Coord(e2);
	float2 ept1 = grid->c_Grid2Space(make_float2(crd1.x,crd1.y));
	float2 ept2 = grid->c_Grid2Space(make_float2(crd2.x,crd2.y));
	float2 normv = make_float2(ept1.x - ept2.x, ept2.y - ept1.y);
	
	vector<float2> edge;
	edge.push_back(ept1);
	edge.push_back(ept2);
	boundaries.push_back(edge);
	
	//mexPrintf("Edge %d %d\n", e1, e2);
	
	int stps = grid->c_spc.x / grid->spc.x;
	int ocell = -1, cell = -1;
	float2 ovec, vec;
	bool failed;
	int cz = 0, cnz = 0;
	float2 pt, opt;
	float2 refdir;
	for (int i = 0; i < stps; i++)
	{	
		// find the point
		opt = pt;
		pt = ept1 + i * (ept2 - ept1) / (stps - 1);
		
		// find the cell
		float2 gpt = grid->c_Space2Grid(pt);
		ocell = cell;
		cell = grid->c_Coord2Addr(make_int2(floor(gpt.x), floor(gpt.y)));
		
		// find the direction
		if (i == 0)
		{
			refdir = ovec = vec = GetClosestDirection(orgline, pt);
		}
		ovec = vec;
		vec = grid->getVectorVecRef(pt, ispos, vec, failed);
		
		/*if ((cell == ocell) &&
			(dot(vec, normv) * dot(ovec, normv) > 0.0))
			continue;*/
			
		if ((i > 0) && (dot(vec, normv) * dot(ovec, normv) < 0.0))
		{
			// there is a tangent point search for it
			float2 s1 = opt;
			float2 v1 = ovec;
			float2 s2 = pt;
			float2 v2 = vec;
			while (length(s1 - s2) > 0.001 * grid->spc.x)
			{
				float2 ns = 0.5 * (s1 + s2);
				float2 nv = grid->getVectorVecRef(ns, ispos, v1, failed);
				
				if (dot(nv, normv) * dot(v1, normv) > 0.0)
				{
					s1 = ns;
					v1 = nv;
				}
				else
				{
					s2 = ns;
					v2 = nv;
				}
			}
			
			// check the point v1
			int checks = 1;
			int results = 0;
			vector<float2> shearline;
			vector<int> cells;
			IntegrateDP3(true, &s1, 0, intg_lng, error_b, ispos, false, shearline, cells, checks, results, refdir, cellset);

			if ((results == 0) && AllSetCovered(cells, cellset))
			{
			boundaries.push_back(shearline);

			
				return true;
			}
		}
		
		// check for the point
		int checks = 1;
		int results = 0;
		vector<float2> shearline;
		vector<int> cells;
		IntegrateDP3(true, &pt, 0, intg_lng, error_b, ispos, false, shearline, cells, checks, results, refdir, cellset);

		if ((results == 0) && AllSetCovered(cells, cellset))
		{
		boundaries.push_back(shearline);
			
			return true;
		}
	}
	
	return false;
}

bool RealExits(vector<float2>& shearline, vector<int>& cells, bool ispos, bool isforward, vector< vector<float2> >& boundaries)
{	
	// find crossed edges and cell set
	set<int> cellset;
	set<int2, int2_compare> crossed;
	int cell1 = -1, cell2 = -1;
	for (int i = 0; i < shearline.size(); i++)
	{
		cell1 = cell2;
		cell2 = cells[i];
		cellset.insert(cell2);
		if ((cell1 != -1) && (cell2 != cell1))
		{
			vector<int> sps = grid->c_GetCellsSharedPoints(cell1, cell2);
			
			if (sps.size() == 1)
			{
				int2 b1 = grid->c_Addr2Coord(cell1);
				int2 b2 = grid->c_Addr2Coord(cell2);
				int cell3 = grid->c_Coord2Addr(b1.x, b2.y);
				int cell4 = grid->c_Coord2Addr(b2.x, b1.y);
			
				cellset.insert(cell3);
				cellset.insert(cell4);
			}
		}
	}
	
	// check that each cell has 2 edges not covered
	bool oht = false;
	for (set<int>::iterator it = cellset.begin(); it != cellset.end(); it++)
	{
		int cell = *it;
		int2 coord = grid->c_Addr2Coord(cell);
		int n1 = grid->c_Coord2Addr(coord.x, coord.y + 1);
		int n2 = grid->c_Coord2Addr(coord.x, coord.y - 1);
		int n3 = grid->c_Coord2Addr(coord.x + 1, coord.y);
		int n4 = grid->c_Coord2Addr(coord.x - 1, coord.y);
		
		int c1 = 0, c2 = 0;
		if (cellset.find(n1) == cellset.end())
			c1++;
		if (cellset.find(n2) == cellset.end())
			c1++;
		if (cellset.find(n3) == cellset.end())
			c2++;
		if (cellset.find(n4) == cellset.end())
			c2++;
			
		if (((c1 == 2) && (c2 == 0)) || ((c2 == 2) && (c1 == 0)))
		{
			oht = true;
			break;
		}
	}
	if (!oht)
		return true;
	
	// find the shared edges between all pairs of cells
	for (set<int>::iterator it1 = cellset.begin(); it1 != cellset.end(); it1++)
	{
		for (set<int>::iterator it2 = cellset.begin(); it2 != cellset.end(); it2++)
		{
			vector<int> spstmp = grid->c_GetCellsSharedPoints(*it1, *it2);
			if (spstmp.size() == 2)
				crossed.insert(make_int2(min(spstmp[0],spstmp[1]), max(spstmp[0],spstmp[1])));
		}
	}
	
	// print the set if you want
	//mexPrintf("cellset size %d\n", cellset.size());
	//mexPrintf("crossed size %d\n", crossed.size());
	int pre = -1;
	for(vector<int>::iterator itr = cells.begin(); itr != cells.end(); itr++)
	{
		if (*itr != pre)
		{
			pre = *itr;
			//mexPrintf("%d ", *itr);
		}
	}
	//mexPrintf("===\n");
	
	// now check the exits on each cell edge
	cell1 = cell2 = -1;
	for (set<int>::iterator it = cellset.begin(); it != cellset.end(); it ++)
	{
		cell2 = *it;
		
		// edges of cell2
		if (crossed.find(make_int2(min(cell2, cell2+grid->c_height),max(cell2, cell2+grid->c_height))) == crossed.end())
		{
			if (IsExitEdge(cell2, cell2+grid->c_height, shearline, ispos, cellset, boundaries))
				return true;
		}
		if (crossed.find(make_int2(min(cell2+grid->c_height, cell2+grid->c_height+1),max(cell2+grid->c_height, cell2+grid->c_height+1))) == crossed.end())
		{
			if (IsExitEdge(cell2+grid->c_height, cell2+grid->c_height+1, shearline, ispos, cellset, boundaries))
				return true;
		}
		if (crossed.find(make_int2(min(cell2+1, cell2+grid->c_height+1),max(cell2+1, cell2+grid->c_height+1))) == crossed.end())
		{
			if (IsExitEdge(cell2+1, cell2+grid->c_height+1, shearline, ispos, cellset, boundaries))
				return true;
		}
		if (crossed.find(make_int2(min(cell2+1, cell2),max(cell2+1, cell2))) == crossed.end())
		{
			if (IsExitEdge(cell2+1, cell2, shearline, ispos, cellset, boundaries))
				return true;
		}
	}
	
	return false;
}

bool FindMapConvergence(bool ispos, bool isforward, vector<float2>& s_shearline, vector<int>& s_cells, vector<float2>& o_shearline, vector<int>& o_cells)
{
	//mexPrintf("Check convergence\n");
	
	// find a cross section
	int r = s_cells.size() - 1;
	int s;
	vector<int> sps;
	while (1)
	{
		// get two cells in a sequence
		int cell2 = s_cells[r];
		int cell1;
		s = r - 1;
		while (s_cells[s] == cell2)
			s--;
		cell1 = s_cells[s];
		
		// find shared edge
		sps = grid->c_GetCellsSharedPoints(cell1, cell2);
		
		if (sps.size() == 2)
		{
			//mexPrintf("cross section cells %d %d\n", cell1, cell2);
			break;
		}
		r--;
		if (r == 0)
		{
			mexPrintf("ERRRROR\n");
			return false;
		}
	}	
		
	// cross section points
	int2 f1 = grid->c_Addr2Coord(sps[0]);
	int2 f2 = grid->c_Addr2Coord(sps[1]);
	//mexPrintf("%d %d\n", f1.x,f1.y);
	//mexPrintf("%d %d\n", f2.x,f2.y);
	float2 e1 = grid->c_Grid2Space(make_float2(f1.x,f1.y));
	float2 e2 = grid->c_Grid2Space(make_float2(f2.x,f2.y));
	double dist, odist;
	
	// find intersection of line with section
	float2 pt;
	float2 dm1, dm2;
	if (intersect2D_Segments(s_shearline[r], s_shearline[s], e1, e2, dm1, dm2) != 1)
	{
		mexPrintf("Convergence map start point is wrong\n");
		return false;
	}
	pt = dm1;
	
	// XY map for prof. Tricoche
	/*//e2 = e2 + 0.009945 * (e2 - e1);
	float2 bla = make_float2(0.5, 0.6);//e1+3*(e1-e2);
	vector<float2> shearline;
	vector<int> cells;
	shearline.push_back(e1);
	shearline.push_back(e2);
	for (int i = 0 ; i < 200; i++)
	{
		float2 pt = bla;
		int checks = 2;
		int results = 0;
		IntegrateDP2(true, &pt, 0, 30*intg_lng, error_b, ispos, true, shearline, cells, checks, results, -normalize(s_shearline[r] - s_shearline[s]), e1, e2);
		bool donotinter = false;
		if (intersect2D_Segments(shearline[shearline.size() - 2], shearline[shearline.size() - 1], e1 + e1 - e2, e2 + e2 - e1, dm1, dm2) != 1)
		{
			donotinter = true;
			mexPrintf("Something wrong\n");
			break;
		}
		
		//printf number
		float xz = length(bla - e2);
		if (dot(bla - e2, e1 - e2) < 0)
			xz *= -1;
		float yz = length(dm1 - e2);
		if (dot(dm1 - e2, e1 - e2) < 0)
			yz *= -1;
		mexPrintf("%f\t%f\t%f\n", xz, yz, length(e1 - pt));
		
		bla = shearline[shearline.size() - 1];
	}*/
	
	// iterate on items
	vector<float2> shearline;
	vector<int> cells;
	//shearline.push_back(e1);
	//shearline.push_back(e2);
	pt = 0.5 * (e1 + e2);
	while (1)
	{
		shearline.clear();
		cells.clear();
		int checks = 2;
		int results = 0;
		
		// integrate from the point until a cycle is found
		//mexPrintf("before %f %f\n", pt.x, pt.y);
		IntegrateDP2(true, &pt, 0, intg_lng, error_b, ispos, true, shearline, cells, checks, results, normalize(s_shearline[r] - s_shearline[s]), e1, e2);
		//mexPrintf("after %f %f %d\n", pt.x, pt.y, results);
		
		// if no cycle 
		if ((results & 2) == 0)
		{
			o_shearline = shearline;
			o_cells = cells;
			return false;
		}

		if (intersect2D_Segments(shearline[shearline.size() - 2], shearline[shearline.size() - 1], e1, e2, dm1, dm2) != 1)
		{
			o_shearline = shearline;
			o_cells = cells;
			return false;
		}
			
		//mexPrintf("intersections %f %f %f %f\n", dm1.x, dm1.y, dm2.x, dm2.y);
			
		odist = length(e1 - e2);
		if (length(dm1 - e1) < length(dm1 - e2))
			e2 = 0.5 * (e1 + e2);
		else
			e1 = 0.5 * (e1 + e2);
		dist = length(e1 - e2);
		
		pt = 0.5 * (e1 + e2);
		
		//mexPrintf("dist %lf\n", dist);
		if (dist > odist)
		{
			return false;
		}
		if (dist < 0.25 * grid->spc.x)
		{
			o_shearline = shearline;
			o_cells = cells;
			return true;
		}
	}
	
	o_shearline = shearline;
	o_cells = cells;
	return true;
	
}

void CheckShearlines(bool ispos, bool isforward, int& streamlineidx, vector< vector<float2> >& shearlines)
{
	set<int> mrk;
	
	int count = 0;
	for (int x = 0; x < grid->c_width; x++)
	{
		for (int y = 0; y < grid->c_height; y++)
		{
			float2 spt = grid->c_Grid2Space(make_float2(x, y));
			
			// coarse resolution check
			vector<float2> shearline;
			vector<int> cells;
			int checks = 1;
			int results = 0;
			IntegrateDP(true, &spt, 0, intg_lng, error_b, ispos, isforward, shearline, cells, checks, results, make_float2(0.0));
			if (results == 0)
				continue;
			count++;
				
			// get the closed part 
			int c2 = cells.size() - 1;
			int c1 = c2 - 1;
			while (cells[c1] == cells[c2])
				c1--;
			while (cells[c1] != cells[c2])
				c1--;
			shearline.erase(shearline.begin(), shearline.begin() + c1);
			cells.erase(cells.begin(), cells.begin() + c1);
				
			// check not in mrkd
			bool ae = false;
			vector<int> finecells;
			for (int i = 0; i < shearline.size(); i++)
			{
				float2 gpt = grid->Space2Grid(shearline[i]);
				finecells.push_back(grid->Coord2Addr(floor(gpt.x), floor(gpt.y)));
				if ((grid->IsBoundary(floor(gpt.x), floor(gpt.y))) || (grid->IsBoundary(ceil(gpt.x), ceil(gpt.y))))
				{
					ae = true;
					break;
				}
			}
			if (ae) continue;
			for (int i = 0; i < finecells.size(); i++)
			{
				if ((i > 0) && (finecells[i] == finecells[i-1]))
					continue;
				if (mrk.find(finecells[i]) != mrk.end())
				{
					ae = true;
					break;
				}
			}
			if (ae) continue;
				
			// set for the cells
			set<int> cellset;
			for (vector<int>::iterator it = cells.begin(); it != cells.end(); it++)
				cellset.insert(*it);
			if (cellset.size() <= 4)
				continue;
				
			// now check potential exits
			vector< vector<float2> > boundaries;
			if (RealExits(shearline, cells, ispos, isforward, boundaries))
				continue;
			
			// check convergence
			vector<float2> o_shearline;
			vector<int> o_cells;
			int conv = FindMapConvergence(ispos, isforward, shearline, cells, o_shearline, o_cells);
			//mexPrintf("\nConvergence %d\n", conv);
			if (!conv)
				continue;
			
			// add cells to marked
			ae = false;
			finecells.clear();
			for (int i = 0; i < o_shearline.size(); i++)
			{
				float2 gpt = grid->Space2Grid(o_shearline[i]);
				finecells.push_back(grid->Coord2Addr(floor(gpt.x), floor(gpt.y)));
			}
			for (int i = 0; i < finecells.size(); i++)
			{
				if ((i > 0) && (finecells[i] == finecells[i-1]))
					continue;
					
				if (mrk.find(finecells[i]) != mrk.end())
				{
					ae = true;
					break;
				}
			}
			if (ae)
				continue;
			for (int i = 0; i < finecells.size(); i++)
			{
				mrk.insert(finecells[i]);
			}
				
			// now make cycle closed
			if (!conv)
				shearlines.push_back(shearline);
			else
				shearlines.push_back(o_shearline);
			streamlineidx++;
			/*for (vector<vector<float2>>::iterator it = boundaries.begin(); it != boundaries.end(); it++)
			{
				shearlines.push_back(*it);
				streamlineidx++;
			}*/
				
			// add the shearline
			mexPrintf("At %d %d (%f, %f) => Shearline %d has size %d and %d distinct cells\n", x, y, spt.x, spt.y, count, shearline.size(), cellset.size());

		}
	}
}


/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	
	
	double* timespan = mxGetPr(prhs[0]);
	double* tmp00 = mxGetPr(prhs[1]);
	float4 domain = make_float4(tmp00[0], tmp00[2], tmp00[1], tmp00[3]);
	UINT64_T* tmp01 = (UINT64_T*) mxGetData(prhs[2]);
	int2 resolution = make_int2(tmp01[0], tmp01[1]);
	double* eigvals = mxGetPr(prhs[3]);
	double* eigvecs = mxGetPr(prhs[4]);
	double* tensors = mxGetPr(prhs[5]);
	int xDim, yDim;
	double *outArray;
	int colLen = 2, rowLen = 2;
	int i, j;

	// print some info
	mexPrintf("Domain %f %f %f %f\n", domain.x, domain.y, domain.z, domain.w);
	mexPrintf("Resolution is %dx%d\n", resolution.x, resolution.y);
	cycleiterc = 2;
	intg_lng = cycleiterc * 2.0 * ((domain.y - domain.x) + (domain.w - domain.z));
	mexPrintf("Integration length is %lf\n", intg_lng);
	
	// create the grid
	grid = new RegularGrid(tensors, eigvals, eigvecs, domain, resolution.x, resolution.y, (domain.y - domain.x) / (resolution.x - 1), (domain.w - domain.z) / (resolution.y - 1));
	//grid->c_SetCoarseRes(resolution.x / 40, resolution.x / 40);
	grid->c_SetCoarseRes(20,20);
	mexPrintf("coarse dimensions are %d %d\n", grid->c_width, grid->c_height);
	grid->GetDoubleGyreCGTensor();
	//grid->LoadTensorFromNrrd();
	//grid->WriteEtaToNrrd();
	mexPrintf("%f %f\n", grid->spc.x, grid->spc.y);
	
	// integrate
	int streamlineidx = 0;
	vector< vector<float2> > shearlines;
	CheckShearlines(true, true, streamlineidx, shearlines); mexPrintf("%.1f%% done\n", 25.0f); mexEvalString("drawnow;");
	CheckShearlines(true, false, streamlineidx, shearlines); mexPrintf("%.1f%% done\n", 50.0f); mexEvalString("drawnow;");
	
	// add the arrays	
	mwSize dims[1];
	dims[0] = shearlines.size();
	mxArray* cellarr = mxCreateCellMatrix(shearlines.size(), 1);
	for (int x = 0; x < shearlines.size(); x++)
	{
		vector<float2> shearline = shearlines[x];
		
		// add new return array
		rowLen = shearline.size();
		
		//Allocate memory and assign output pointer
		mxArray* arr = mxCreateDoubleMatrix(rowLen, colLen, mxREAL);
		
		//Get a pointer to the data space in our newly allocated memory
		outArray = mxGetPr(arr);
		
		//Copy matrix while multiplying each point by 2
		for(i=0;i<rowLen;i++)
		{
			outArray[i + 0 * rowLen] = shearline[i].x;
			outArray[i + 1 * rowLen] = shearline[i].y;
		}
		
		mxSetCell(cellarr, x, arr);
	}
	plhs[0] = cellarr;
	
	// integrate
	streamlineidx = 0;
	shearlines.clear();
	CheckShearlines(false, true, streamlineidx, shearlines); mexPrintf("%.1f%% done\n", 75.0f); mexEvalString("drawnow;");
	CheckShearlines(false, false, streamlineidx, shearlines); mexPrintf("%.1f%% done\n", 100.0f); mexEvalString("drawnow;");
	
	// add the arrays	
	dims[0] = shearlines.size();
	cellarr = mxCreateCellMatrix(shearlines.size(), 1);
	for (int x = 0; x < shearlines.size(); x++)
	{
		vector<float2> shearline = shearlines[x];
		
		// add new return array
		rowLen = shearline.size();
		
		//Allocate memory and assign output pointer
		mxArray* arr = mxCreateDoubleMatrix(rowLen, colLen, mxREAL);
		
		//Get a pointer to the data space in our newly allocated memory
		outArray = mxGetPr(arr);
		
		//Copy matrix while multiplying each point by 2
		for(i=0;i<rowLen;i++)
		{
			outArray[i + 0 * rowLen] = shearline[i].x;
			outArray[i + 1 * rowLen] = shearline[i].y;
		}
		
		mxSetCell(cellarr, x, arr);
	}
	plhs[1] = cellarr;
			
			
    return;
}
