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
#include <string>
#include <set>
#include <map>
#include <vector>
#include <vector_types.h>
#include <vector_functions.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
	 
#include "RegularGrid.h"
#include <omp.h>

using namespace std;

RegularGrid* grid;

int cycleiterc = 8; // not less than 3
double intg_lng = 100;
double error_s = 1e-10;
double error_b = 1e-8;

float2 GetFlowAt(float time, float2 spt, float2 ref, bool& failed, bool ispos)
{
    return grid->getVector(spt, ispos, ref, failed);
}

bool CycleFound(vector<int>& cells)
{
	// find occurances
	int c2 = cells.size() - 1;
	int c0, c1;
	int i;

	//mexPrintf("%d\n", c2);
	
	// find c1
	i = c2 - 1;
	while (cells[i] == cells[c2])
		i--;
	while (cells[i] != cells[c2])
		i--;
	c1 = i;
	
	//mexPrintf("%d\n", c1);
	
	// find c0
	i = c1 - 1;
	while (cells[i] == cells[c2])
		i--;
	while (cells[i] != cells[c2])
		i--;
	c0 = i;
		
	//mexPrintf("%d\n", c0);
	
	
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
	
	//mexPrintf("%d %d \n", l1.size(), l2.size());
	
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
	p->ref = v;
	
	//float2 gpt = grid->Space2Grid(yp);
	//mexPrintf("%f %f %f %f %f %f %f %f\n", yp.x, yp.y, v.x, v.y, p->ref.x, p->ref.y, gpt.x, gpt.y);
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
			visited[cell] ++;
		
		// now do the checks:
		if (checks & 1)
		{
			if (visited[cell] >= cycleiterc)
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

void IntegrateDP2(bool coarse, float2* y, float st, float et, double mu, bool ispos, bool isforward, vector<float2>& shearline, vector<int>& cells, int& checks, int& results, float2 e1, float2 e2)
{
	int cell = -1, ocell = -1;
	
	const gsl_odeiv_step_type* T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 2);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(mu, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(2);

	paramsDP params = {ispos, isforward, make_float2(0.0)};
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
		   break;
		   
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
			// mexPrintf("== edge %f %f %f %f\n", e1.x, e1.y, e2.x, e2.y);
			// mexPrintf("== cells %d %d\n", ocell, cell);
			// mexPrintf("== itr %d\n", itr);
			//mexPrintf("== time %f\n", t);
			//mexPrintf("== dm1 %f %f %d\n", dm1.x, dm1.y, ec? 1:0);
			//mexPrintf("== itr %d, %f %f %f %f\n", itr, oy.x, oy.y, (*y).x, (*y).y);
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
	while (t < t1)
	{   
		// add the current point
		shearline.push_back(*y);
		
		// get the current point cell
		float2 gpt = coarse? grid->c_Space2Grid(*y) : grid->Space2Grid(*y);
		int2 igpt = make_int2(floor(gpt.x), floor(gpt.y));
		ocell = cell;
		cell = coarse? grid->c_Coord2Addr(igpt) : grid->Coord2Addr(igpt);
		cells.push_back(cell);
		if (ocell != cell)
			visited[cell]++;
		// if (cell !=ocell)
		// 	mexPrintf("%d %d, ", cell, visited[cell]);
		
		// find the first cell in the set
		if ((itr > 3) && (firstcell == -1) && (icells.find(cell) != icells.end()))
		{
			firstcell = cell;
			// mexPrintf(" -> first cell ");
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
				if ((intersect2D_Segments(*y, oy, c30, c31, dm1, dm2) != 0) && (icells.find(cell3) == icells.end()))
				{
					results = 1;
					break;
				}
				
				// check cell4
				vector<int> spt4 = coarse? grid->c_GetCellsSharedPoints(cell, cell4) : grid->GetCellsSharedPoints(cell, cell4);
				float2 c40 = coarse? grid->c_Cell2Space(spt4[0]) : grid->Cell2Space(spt4[0]);
				float2 c41 = coarse? grid->c_Cell2Space(spt4[1]) : grid->Cell2Space(spt4[1]);
				if ((intersect2D_Segments(*y, oy, c40, c41, dm1, dm2) != 0) && (icells.find(cell4) == icells.end()))
				{
					results = 1;
					break;
				}	
			}
		}
		
		// if the cell is not in the set exit
		if ((itr > 3) && (icells.find(cell) == icells.end()))
		{
			results = 1;
			break;
		}
		
		// if the first cell is hit again exit
		if ((firstcell != -1) && (visited[firstcell] >= 2))
		{
			results = 0;
			break;
		}
		
		// next point
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, ydp);
		oy = *y;
		*y = make_float2(ydp[0], ydp[1]);
		
		if (status != GSL_SUCCESS)
		   break;
		itr++;
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	free(visited);
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

bool IsExitEdge(int e1, int e2, float2 ref, bool ispos, set<int>& cellset, vector<vector<float2> >& boundaries)
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
	
	int stps = grid->c_spc.x / grid->spc.x;
	int ocell = -1, cell = -1;
	float2 ovec = ref, vec = ref;
	bool failed;
	int cz = 0, cnz = 0;
	float2 pt, opt;
	for (int i = 0; i < stps; i++)
	{
		opt = pt;
		pt = ept1 + i * (ept2 - ept1) / (stps - 1);
		float2 gpt = grid->c_Space2Grid(pt);
		ocell = cell;
		cell = grid->c_Coord2Addr(make_int2(floor(gpt.x), floor(gpt.y)));
		ovec = vec;
		vec = grid->getVector(pt, ispos, vec, failed);
		
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
				float2 nv = grid->getVector(ns, ispos, v1, failed);
				
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
			
			float2 s1c = s1;
			
			// check the point v1
			int checks = 1;
			int results = 0;
			vector<float2> shearline;
			vector<int> cells;
			IntegrateDP3(true, &s1, 0, intg_lng, error_b, ispos, false, shearline, cells, checks, results, v1, cellset);
			// mexPrintf("\n");
			if (results == 0)
			{
			boundaries.push_back(shearline);
			/*vector<float2> bshearline;
			vector<int> bcells;
			IntegrateDP3(true, &s1c, 0, intg_lng, error_b, ispos, true, bshearline, bcells, checks, results, v1, cellset);
			boundaries.push_back(bshearline);*/
			
				return true;
			}
		}
			
		float2 ptc = pt;
		
		// check for the point
		int checks = 1;
		int results = 0;
		vector<float2> shearline;
		vector<int> cells;
		IntegrateDP3(true, &pt, 0, intg_lng, error_b, ispos, false, shearline, cells, checks, results, vec, cellset);
		// mexPrintf("\n");
		if (results == 0)
		{
		boundaries.push_back(shearline);
		/*vector<float2> bshearline;
		vector<int> bcells;
		IntegrateDP3(true, &ptc, 0, intg_lng, error_b, ispos, true, bshearline, bcells, checks, results, vec, cellset);
		boundaries.push_back(bshearline);*/
			
			return true;
		}
	}
	
	//mexPrintf("count is %d/%d\n", cz, cnz);
	
	return false;
}

bool RealExits(vector<float2>& shearline, vector<int>& cells, bool ispos, bool isforward, vector<vector<float2> >& boundaries)
{	
	// find crossed edges and cell set
	vector<int> acells;
	set<int> cellset;
	map<int,int> a2n;
	set<int2, int2_compare> crossed;
	int cell1 = -1, cell2 = -1;
	for (int i = 0; i < shearline.size(); i++)
	{
		cell1 = cell2;
		cell2 = cells[i];
		cellset.insert(cell2);
		if ((i > 0) && (cell2 != cell1))
		{
			vector<int> sps = grid->c_GetCellsSharedPoints(cell1, cell2);
			
			if (sps.size() == 1)
			{
				int2 b1 = grid->c_Addr2Coord(cell1);
				int2 b2 = grid->c_Addr2Coord(cell2);
				int cell3 = grid->c_Coord2Addr(b1.x, b2.y);
				int cell4 = grid->c_Coord2Addr(b2.x, b1.y);
				
				a2n.insert(pair<int,int>(acells.size(), i));
				acells.push_back(cell3);
				a2n.insert(pair<int,int>(acells.size(), i));
				acells.push_back(cell4);
				
				cellset.insert(cell3);
				cellset.insert(cell4);
			}
			else if (sps.size() == 2)
			{
				crossed.insert(make_int2(min(sps[0],sps[1]), max(sps[0],sps[1])));
			}
		}
		
		if (cell2 != cell1)
		{
			a2n.insert(pair<int,int>(acells.size(), i));
			acells.push_back(cell2);
		}
	}
	
	// find the shared edges between all pairs of cells
	for (set<int>::iterator it1 = cellset.begin(); it1 != cellset.end(); it1++)
	{
		for (set<int>::iterator it2 = cellset.begin(); it2 != cellset.end(); it2++)
		{
			vector<int> spstmp;
			spstmp = grid->c_GetCellsSharedPoints(*it1, *it2);
			if (spstmp.size() == 2)
				crossed.insert(make_int2(min(spstmp[0],spstmp[1]), max(spstmp[0],spstmp[1])));
		}
	}
	
	//mexPrintf("cellset size %d\n", cellset.size());
	//mexPrintf("crossed size %d\n", crossed.size());

	// print the set if you want
	for(vector<int>::iterator itr = acells.begin(); itr != acells.end(); itr++)
		// mexPrintf("%d ", *itr);
	// mexPrintf("===\n");
	
	// now check the exits on each cell edge
	cell1 = cell2 = -1;
	for (int k = 0; k < acells.size(); k++)
	{
		cell1 = cell2;
		cell2 = acells[k];
		int i = a2n[k];
		//mexPrintf("cell idxs %d %d\n", i, k);
		if (cell2 != cell1)
		{
			// edges of cell2
			if (crossed.find(make_int2(min(cell2, cell2+grid->c_height),max(cell2, cell2+grid->c_height))) == crossed.end())
			{
				if (IsExitEdge(cell2, cell2+grid->c_height, shearline[i] - shearline[i-1], ispos, cellset, boundaries))
					return true;
			}
			if (crossed.find(make_int2(min(cell2+grid->c_height, cell2+grid->c_height+1),max(cell2+grid->c_height, cell2+grid->c_height+1))) == crossed.end())
			{
				if (IsExitEdge(cell2+grid->c_height, cell2+grid->c_height+1, shearline[i] - shearline[i-1], ispos, cellset, boundaries))
					return true;
			}
			if (crossed.find(make_int2(min(cell2+1, cell2+grid->c_height+1),max(cell2+1, cell2+grid->c_height+1))) == crossed.end())
			{
				if (IsExitEdge(cell2+1, cell2+grid->c_height+1, shearline[i] - shearline[i-1], ispos, cellset, boundaries))
					return true;
			}
			if (crossed.find(make_int2(min(cell2+1, cell2),max(cell2+1, cell2))) == crossed.end())
			{
				if (IsExitEdge(cell2+1, cell2, shearline[i] - shearline[i-1], ispos, cellset, boundaries))
					return true;
			}
		}
	}
	
	// iterate
	return false;
	
}

bool FindMapConvergence(bool ispos, bool isforward, vector<float2>& s_shearline, vector<int>& s_cells, vector<float2>& o_shearline, vector<int>& o_cells)
{
	// mexPrintf("Check convergence\n");
	
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
			// mexPrintf("ERRRROR\n");
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
		// mexPrintf("Convergence map start point is wrong\n");
		return false;
	}
	pt = dm1;
	
	// iterate on items
	vector<float2> shearline;
	vector<int> cells;
	while (1)
	{
		shearline.clear();
		cells.clear();
		int checks = 2;
		int results = 0;
		
		// integrate from the point until a cycle is found
		// mexPrintf("before %f %f\n", pt.x, pt.y);
		IntegrateDP2(true, &pt, 0, intg_lng, error_b, ispos, isforward, shearline, cells, checks, results, e1, e2);
		// mexPrintf("after %f %f %d\n", pt.x, pt.y, results);
		
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
			
		// mexPrintf("intersections %f %f %f %f\n", dm1.x, dm1.y, dm2.x, dm2.y);
			
		odist = length(e1 - e2);
		if (length(dm1 - e1) < length(dm1 - e2))
			e2 = 0.5 * (e1 + e2);
		else
			e1 = 0.5 * (e1 + e2);
		dist = length(e1 - e2);
		
		pt = shearline[shearline.size() - 1];
		
		// mexPrintf("dist %lf\n", dist);
		if (dist > odist)
		{
			// mexPrintf("What the heq!!\n");
			return false;
		}
		if (dist < 0.5 * grid->spc.x)
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

void CheckShearlines(bool ispos, bool isforward, int& streamlineidx, vector<vector<float2> >& shearlines)
{
	int count = 0;
	for (int x = 0; x < grid->c_width; x++)
	{
		//#pragma omp parallel for schedule(dynamic,25)
		for (int y = 0; y < grid->c_height; y++)
		{
			//if ((x != 16) || (y != 13))
				//continue;
		
		
			float2 spt = grid->c_Grid2Space(make_float2(x,y));
			
			vector<float2> shearline;
			vector<int> cells;
			int checks = 1;
			int results = 0;
			
			// coarse resolution check
			IntegrateDP(true, &spt, 0, intg_lng, error_b, ispos, isforward, shearline, cells, checks, results, make_float2(0.0));
			if (results == 0)
				continue;
			
			// get the closed part 
			int c2 = cells.size() - 1;
			int c1 = c2 - 1;
			while (cells[c1] == cells[c2])
				c1--;
			while (cells[c1] != cells[c2])
				c1--;
			shearline.erase(shearline.begin(), shearline.begin() + c1);
			cells.erase(cells.begin(), cells.begin() + c1);
			
			// set for the cells
			set<int> cellset;
			for (vector<int>::iterator it = cells.begin(); it != cells.end(); it++)
				cellset.insert(*it);
			if (cellset.size() <= 3)
				continue;
			
			// check convergence
			/*vector<float2> o_shearline;
			vector<int> o_cells;
			int conv = FindMapConvergence(ispos, isforward, shearline, cells, o_shearline, o_cells);
			mexPrintf("Convergence %d\n", conv);*/
			//if (!conv)
				//continue;
				
			//if (o_shearline.size() > 3000)
			//	continue;
				
			count++;
			
			// now check potential exits
			vector<vector<float2> > boundaries;
			if (RealExits(shearline, cells, ispos, isforward, boundaries))
				continue;

			
			//
			//if (count == 6)
			{
				shearlines.push_back(shearline);
				streamlineidx++;
				/*for (vector<vector<float2>>::iterator it = boundaries.begin(); it != boundaries.end(); it++)
				{
					shearlines.push_back(*it);
					streamlineidx++;
				}*/
				
				//return;
			}
			
			// add the shearline
			// mexPrintf("\nAt %d %d => Shearline %d has size %d and %d distinct cells\n", x, y, count, shearline.size(), cellset.size());
			// mexPrintf("=======================================\n");
			//return;
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
	int xDim, yDim;
	double *outArray;
	int colLen = 2, rowLen = 2;
	int i, j;
	
	
	// get dimensions of any array
	//xDim = (int) mxGetM(prhs[2]);
	//yDim = (int) mxGetN(prhs[2]);
	//mexPrintf("x Dimensions = %d.\n",xDim);
	//mexPrintf("y Dimensions = %d.\n",yDim);
	
	// test type of array data
	//mexPrintf("test = %d.\n", mxIsDouble(prhs[4]));
   
	// get number of parameters
	//mexPrintf("Number of parameters is %d\n", nrhs);
	
	grid = new RegularGrid(eigvals, eigvecs, resolution.x, resolution.y, (domain.y - domain.x) / (resolution.x - 1), (domain.w - domain.z) / (resolution.y - 1));
	//grid->c_SetCoarseRes(resolution.x / 10, resolution.y / 10);
	grid->c_SetCoarseRes(10.0, 10.0);
	// mexPrintf("coarse dimensions are %d %d\n", grid->c_width, grid->c_height);
	
	// integrate
	int streamlineidx = 0;
	vector<vector<float2> > shearlines;
	CheckShearlines(true, true, streamlineidx, shearlines);
	//CheckShearlines(true, false, streamlineidx, shearlines);
	
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
	//CheckShearlines(false, true, streamlineidx, shearlines);
	//CheckShearlines(false, false, streamlineidx, shearlines);
	
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
