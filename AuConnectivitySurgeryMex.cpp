/*
 * AuConnectivitySurgeryMex.cpp
 * Usage: [Skel.E, history] = AuConnectivitySurgery(X, Y, Z, TRIV)
 * Performs the connectivity surgery step by Au et al. 2008 to convert a
 * collapsed mesh into a 1D skeleton.
 * Variables:
 * X - X coordinate for all mesh vertices.
 * Y - Y coordinate for all mesh vertices.
 * Z - Z coordinate for all mesh vertices.
 * TRIV - mesh triangles.
 *
 * David Pickup 2013
 */

#include "mex.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>

// Define how to access arrays with two dimentions in matlab.
#define POS(x,y,M) ((int)(x + (y*M)))

// Compute distance between a vertex and all its neighbours.
double distances(const double *X, const double *Y, const double *Z,
    const double x, const double y, const double z,
        const std::vector<int>& vs)
{
    int i;
    double D, a, b, c;
    
    D = 0;
    for (i = 0; i < vs.size(); i++)
    {
        a = x - X[vs[i]];
        b = y - Y[vs[i]];
        c = z - Z[vs[i]];
        D += sqrt((a*a) + (b*b) + (c*c));
    }
    
    return D;
}

// Compute the collapse cost for a half-edge.
double collapseCost(const double *X, const double *Y, const double *Z,
    const int idx1, const int idx2, const double D, double ***Q)
{
    double cost, tmp, x, y, z;
    
    x = X[idx1];
    y = Y[idx1];
    z = Z[idx1];
    
    // Compute the shape cost.
    cost = x * ((x*Q[idx1][0][0]) + (y*Q[idx1][1][0]) + (z*Q[idx1][2][0]) + (Q[idx1][3][0]));
    cost += y * ((x*Q[idx1][0][1]) + (y*Q[idx1][1][1]) + (z*Q[idx1][2][1]) + (Q[idx1][3][1]));
    cost += z * ((x*Q[idx1][0][2]) + (y*Q[idx1][1][2]) + (z*Q[idx1][2][2]) + (Q[idx1][3][2]));
    cost += ((x*Q[idx1][0][3]) + (y*Q[idx1][1][3]) + (z*Q[idx1][2][3]) + (Q[idx1][3][3]));
    
    cost += x * ((x*Q[idx2][0][0]) + (y*Q[idx2][1][0]) + (z*Q[idx2][2][0]) + (Q[idx2][3][0]));
    cost += y * ((x*Q[idx2][0][1]) + (y*Q[idx2][1][1]) + (z*Q[idx2][2][1]) + (Q[idx2][3][1]));
    cost += z * ((x*Q[idx2][0][2]) + (y*Q[idx2][1][2]) + (z*Q[idx2][2][2]) + (Q[idx2][3][2]));
    cost += ((x*Q[idx2][0][3]) + (y*Q[idx2][1][3]) + (z*Q[idx2][2][3]) + (Q[idx2][3][3]));
    
    
    // Compute the sampling cost.
    x = x - X[idx2];
    y = y - Y[idx2];
    z = z - Z[idx2];
    tmp = sqrt((x*x) + (y*y) + (z*z));
    cost += tmp * (D-tmp);
    
    return cost;
}

/* The matlab gateway function. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Varify that the correct number of function arguments have been given.
    if (nrhs != 6)
        mexErrMsgTxt("Six input arguments required.");
//     if (nlhs != 2)
//         mexErrMsgTxt("Two output argument required.");
    
    // Initialise variables.
    double *X;                              // Mesh vertices X coordinates.
    double *Y;                              // Mesh vertices Y coordinates.
    double *Z;                              // Mesh vertices Z coordinates.
    double *TRIV;                           // Mesh triangles.
    double *I;                              // i index of non-zeros adjacency matrix entries.
    double *J;                              // j index of non-zeros adjacency matrix entries.
    double ***Q;                            // Error matrices.
    std::vector<std::vector<int> > Edges;   // Edges (adjacency list).
    std::vector<int> EdgeIdx;               // Edge indices.
    std::vector<std::vector<int> > TVAL;    // Triangle-vertex adjacency list.
    std::vector<std::vector<int> > VTAL;    // Vertex-triangle adjacency list.
    std::vector<std::vector<int> >::iterator TVAL_Iterator;
    std::vector<std::vector<double> > Costs;// Cost of each edge.
    std::vector<std::vector<int> > history; // History of edge collapses.
    double a[3];                            // Normalised half-edge vector.
    double b[3];                            // a x v.
    double *output;
    double norm, cost, minCost, x, y, z, tmp, D;
    int minCostEdge[2];
    int nPnts, nTris, nEdges, flag;
    int i, j, k, l, idx, idx2, prev_j, v1, v2, count;
    double *debug;
    time_t  timer;
    mxArray *matlabArray;
    
    // Get triangles and vertex coordinates.
    X = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    Z = mxGetPr(prhs[2]);
    nPnts = mxGetM(prhs[0]);
    TRIV = mxGetPr(prhs[3]);
    nTris = mxGetM(prhs[3]);
    I = mxGetPr(prhs[4]);
    J = mxGetPr(prhs[5]);
    nEdges = std::max(mxGetM(prhs[4]),mxGetN(prhs[4]));
        
    // Initialise triangle-vertex adjacency list.
    TVAL.resize(nTris);
    for (i = 0; i < nTris; i++)
        for (j = 0; j < 3; j++)
            TVAL[i].push_back(TRIV[POS(i,j,nTris)]-1);
    
    // Initialise vertex-triangle adjacency list.
    VTAL.resize(nPnts);
    for (i = 0; i < nTris; i++)
    {
        for (j = 0; j < 3; j++)
        {
            idx = TRIV[POS(i,j,nTris)]-1;
            VTAL[idx].push_back(i);
        }
    }
    
    // Initialise adjacency list of edges.
    Edges.resize(nPnts);
    Costs.resize(nPnts);
    history.resize(nPnts);
    for (i = 0; i < nPnts; i++)
        EdgeIdx.push_back(i);
    
    // Compute adjacency list of edges.
    for (i = 0; i < nEdges; i++)
    {
        Edges[I[i]-1].push_back(J[i]-1);
        Costs[I[i]-1].push_back(0);
    }
    
    // Allocate memory for error matrices.
    Q = (double***)calloc(nPnts,sizeof(double**));
    for (i = 0; i < nPnts; i++)
    {
        Q[i] = (double**)calloc(4,sizeof(double*));
        for (j = 0; j < 4; j++)
            Q[i][j] = (double*)calloc(4,sizeof(double));
    }
    
    // Compute error matrix for each half edge.
    for (i = 0; i < Edges.size(); i++)
    {
        for (j = 0; j < Edges[i].size(); j++)
        {
            // Compute normalised half-edge vector.
            a[0] = X[Edges[i][j]] - X[i];
            a[1] = Y[Edges[i][j]] - Y[i];
            a[2] = Z[Edges[i][j]] - Z[i];
            norm = sqrt((a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]));
            a[0] = a[0] / norm;
            a[1] = a[1] / norm;
            a[2] = a[2] / norm;
            
            // Compute cross product between vector and vertex.
            b[0] = (a[1]*Z[i]) - (a[2]*Y[i]);
            b[1] = (a[2]*X[i]) - (a[0]*Z[i]);
            b[2] = (a[0]*Y[i]) - (a[1]*X[i]);
            
            // Update error matrix.
            Q[i][0][0] += (a[2]*a[2]) + (a[1]*a[1]);
            Q[i][0][1] += -(a[0]*a[1]);
            Q[i][0][2] += -(a[0]*a[2]);
            Q[i][0][3] += (b[2]*a[1]) - (b[1]*a[2]);
            Q[i][1][0] += -(a[0]*a[1]);
            Q[i][1][1] += (a[2]*a[2]) + (a[0]*a[0]);
            Q[i][1][2] += -(a[1]*a[2]);
            Q[i][1][3] += (b[0]*a[2]) - (b[2]*a[0]);
            Q[i][2][0] += -(a[0]*a[2]);
            Q[i][2][1] += -(a[1]*a[2]);
            Q[i][2][2] += (a[1]*a[1]) + (a[0]*a[0]);
            Q[i][2][3] += (b[1]*a[0]) - (b[0]*a[1]);
            Q[i][3][0] += (b[2]*a[1]) - (b[1]*a[2]);
            Q[i][3][1] += (b[0]*a[2]) - (b[2]*a[0]);
            Q[i][3][2] += (b[1]*a[0]) - (b[0]*a[1]);
            Q[i][3][3] += (b[0]*b[0]) + (b[1]*b[1]) + (b[2]*b[2]);
        }
    }
    
    prev_j = -1;
    
    // Collapse edges until all faces have been removed.
    while (nTris > 0)//(!TVAL.empty())
    {
        if (prev_j == -1)
        {
            // Iterate through all vertices.
            for (i = 0; i < Edges.size(); i++)
            {
                // Compute length of all edges attached to current vertex.
                D = distances(X,Y,Z,X[i],Y[i],Z[i],Edges[i]);
                
                // Iterate through all half edges computing their collpase cost.
                for (j = 0; j < Edges[i].size(); j++)
                {
                    idx = Edges[i][j];
                    Costs[i][j] = collapseCost(X,Y,Z,idx,i,D,Q);
                }
            }
        }
        else
        {
            // Compute length of all edges attached to current vertex.
            D = distances(X,Y,Z,X[prev_j],Y[prev_j],Z[prev_j],Edges[prev_j]);
            
            // Iterate through all half edges computing their collpase cost.
            for (i = 0; i < Edges[prev_j].size(); i++)
            {
                idx = Edges[prev_j][i];
                Costs[prev_j][i] = collapseCost(X,Y,Z,idx,prev_j,D,Q);
            }
                
            // Iterate though all vertices connected to prev_j.
            for (i = 0; i < Edges[prev_j].size(); i++)
            {
                idx = Edges[prev_j][i];
                
                // Compute length of all edges attached to current vertex.
                D = distances(X,Y,Z,X[idx],Y[idx],Z[idx],Edges[idx]);
                
                // Iterate through all half edges computing their collpase cost.
                for (j = 0; j < Edges[idx].size(); j++)
                {
                    idx2 = Edges[idx][j];
                    Costs[idx][j] = collapseCost(X,Y,Z,idx2,idx,D,Q);
                }
            }
        }
        
        // Find lowest cost.
        minCost = pow(10,10);
        minCostEdge[0] = -1;
        minCostEdge[1] = -1;
        for (i = 0; i < Costs.size(); i++)
        {
            if (EdgeIdx[i] == -1)
                continue;
            for (j = 0; j < Costs[i].size(); j++)
            {
                if (Costs[i][j] < minCost)
                {
                    
                    flag = 0;
                    //Check for triangles connected to current edge.
                    for (k = 0; k < VTAL[i].size(); k++)
                        for (l = 0; l < VTAL[Edges[i][j]].size(); l++)
                            if (VTAL[i][k] == VTAL[Edges[i][j]][l])
                                flag = 1;
                    if (flag == 0)
                        continue;
                    
                    minCost = Costs[i][j];
                    minCostEdge[0] = i;
                    minCostEdge[1] = j;
                }
            }
        }
        
        // Get vertex indices of current edge.
        v1 = minCostEdge[0];
        v2 = Edges[minCostEdge[0]][minCostEdge[1]];
        prev_j = v2;
        
        // Update edge collapse history.
        history[v2].push_back(v1);
        for (i = 0; i < history[v1].size(); i++)
            history[v2].push_back(history[v1][i]);
        std::vector<int>().swap(history[v1]);
        

        
        // Delete all triangles connected to current edge.
        for (i = 0; i < VTAL[v1].size(); i++)
        {
            flag = 0;
            idx = VTAL[v1][i];
            for (j = 0; j < VTAL[v2].size(); j++)
            {
                // If true then triangle is connected to edge to be
                // collapsed.
                if (VTAL[v1][i] == VTAL[v2][j])
                {
                    // Iterate through all the vertices that are part of
                    // the triangle.
                    for (k = 0; k < TVAL[idx].size(); k++)
                    {
                        // Get the index of the current vertex connected to
                        // the triangle.
                        idx2 = TVAL[idx][k];
                        
                        // If the vertex is not part of the edge to be
                        // colapsed.
                        if ((idx2 != v1) && (idx2 != v2))
                        {
                            // Iterate through all the triangles connected
                            // to the vertex.
                            for (l = 0; l < VTAL[idx2].size(); l++)
                            {
                                // Delete the reference to the triangle if
                                // it matches.
                                if (VTAL[idx2][l] == idx)
                                {
                                    VTAL[idx2].erase(VTAL[idx2].begin()+l);
                                    break;
                                }
                            }
                        }
                    }
                    
                    // Delete the reference to the triangle.
                    VTAL[v1].erase(VTAL[v1].begin()+i);
                    i --;
                    VTAL[v2].erase(VTAL[v2].begin()+j);
                    j --;
                    
                    nTris --;
                    flag = 1;   // Remember that the triangle was deleted.
                    break;
                }
            }
            
            // If the triangle was not deleted.
            if (flag == 0)
            {
                // Substiture reference in TVAL to i for j.
                for (j = 0; j < TVAL[idx].size(); j++)
                {
                    if (TVAL[idx][j] == v1)
                    {
                        TVAL[idx][j] = v2;
                        break;
                    }
                }

                // Add reference from j to triangle.
                VTAL[v2].push_back(idx);
            }
        }
        
        
        
        // Iterate through all vertices connected to vertex to be deleted (i).
        for (i = 0; i < Edges[v1].size(); i++)
        {
            // Get vertex index for current connected vertex.
            idx = Edges[v1][i];
            
            // Connect vertex to j, if not equal to j.
            if (idx != v2)
            {
                flag = 1;
                for (j = 0; j < Edges[v2].size(); j++)
                {
                    if (Edges[v2][j] == idx)
                        flag = 0;
                }
                if (flag == 1)
                {
                    Edges[v2].push_back(idx);
                    Costs[v2].push_back(pow(10,9));
//                     mexPrintf("Pushing %d onto %d\n",idx,v2);
                }
            }
            
            // Iterate through all its connections.
            for (j = 0; j < Edges[idx].size(); j++)
            {
                // Subtitute connections to i for connections to j.
                if (Edges[idx][j] == v1)
                {
                    // If current vertex is j, then delete connection to i.
                    if (idx == v2)
                    {
                        Edges[idx].erase(Edges[idx].begin()+j);
                        Costs[idx].erase(Costs[idx].begin()+j);
                        break;
                    }
                    
                    // If already connected to j, then delete connection
                    // to i.
                    flag = 1;
                    for (k = 0; k < Edges[idx].size(); k++)
                    {
                        if (Edges[idx][k] == v2)
                        {
                            Edges[idx].erase(Edges[idx].begin()+j);
                            Costs[idx].erase(Costs[idx].begin()+j);
                            flag = 0;
                            break;
                        }
                    }
                    if (flag == 0)
                        break;
                    
                    // Make the substitution.
//                     mexPrintf("Subtituting %d for %d\n",Edges[idx][j],v2);
                    Edges[idx][j] = v2;
                    break;
                }
            }
        }
        
        // Mark collapsed edges as deleted.
        EdgeIdx[v1] = -1;
        

    }
    
    plhs[0] = mxCreateDoubleMatrix(nEdges, 2, mxREAL);
    output = mxGetPr(plhs[0]);
    idx = 0;
    for (i = 0; i < Edges.size(); i++)
    {
        if (EdgeIdx[i] != -1)
        {
            for (j = 0; j < Edges[i].size(); j++)
            {
                output[POS(idx,0,nEdges)] = i+1;
                output[POS(idx,1,nEdges)] = Edges[i][j]+1;
                idx ++;
            }
        }
    }
    
    // Initialise cell array for output images.
    plhs[1] = mxCreateCellMatrix((mwSize)history.size(), 1);
    for (i = 0; i < history.size(); i++)
    {
        matlabArray = mxCreateDoubleMatrix(history[i].size(),1,mxREAL);
        output = mxGetPr(matlabArray);
        for (j = 0; j < history[i].size(); j++)
            output[j] = history[i][j] + 1;
        mxSetCell(plhs[1], (mwIndex)i, matlabArray);
    }
}