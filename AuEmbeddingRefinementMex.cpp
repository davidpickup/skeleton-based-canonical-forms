/*
 * AuEmbeddingRefinementMex.cpp
 * Usage: TODO
 * Refines the embedding of a skeleton using the method by
 * Au et al. SIGGRAPH 2008.
 * Variables:
 * TODO
 *
 * David Pickup 2014
 */

#include "mex.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

// Define how to access arrays with two dimentions in matlab.
#define POS(x,y,M) ((int)(x + (y*M)))

std::ofstream debug; int isdebug = 0;

/* Computes the centeredness of a joint as the standard deviation of the */
/* distances between the joint position and its assigned vertices. */
double centeredness(const double NX, const double NY, const double NZ,
        const double *MX, const double *MY, const double *MZ,
        const std::vector<int> &H)
{
    double D[H.size()];
    double mu = 0;
    
    for (int i = 0; i < H.size(); i++)
    {
        D[i] = sqrt(pow((NX-MX[H[i]]),2) + pow((NY-MY[H[i]]),2) +
                pow((NZ-MZ[H[i]]),2));
        mu += D[i];
    }
    mu = mu / H.size();
    
    double sigma = 0;
    for (int i = 0; i < H.size(); i++)
        sigma += pow((D[i] - mu),2);
    return sqrt(sigma/(double)H.size());
}

/* The matlab gateway function. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Varify that the correct number of function arguments have been given.
    if (nrhs != 11)
        mexErrMsgTxt("Eleven input arguments required.");
    if (nlhs != 3)
        mexErrMsgTxt("Three output argument required.");
    
    // Initialise variables.
    double *MX;                             // Mesh vertices X coordinates.
    double *MY;                             // Mesh vertices Y coordinates.
    double *MZ;                             // Mesh vertices Z coordinates.
    double *SX;                             // Skeleton nodes X coordinates.
    double *SY;                             // Skeleton nodes Y coordinates.
    double *SZ;                             // Skeleton nodes Z coordinates.
    double *MI;                             // i index of non-zeros mesh adjacency matrix entries.
    double *MJ;                             // j index of non-zeros mesh adjacency matrix entries.
    double *SI;                             // i index of non-zeros skeleton adjacency matrix entries.
    double *SJ;                             // j index of non-zeros skeleton adjacency matrix entries.
    std::vector<std::vector<int> > MEdges;  // Edges of mesh (adjacency list).
    std::vector<std::vector<int> > SEdges;  // Edges of skeleton (adjacency list).
    int nNodes;                             // Number of skeleton nodes.
    int nPnts;                              // Number of mesh vertices.
    int nMEdges;                            // Number of mesh edges.
    int nSEdges;                            // Number of skeleton edges (bones).
    std::vector<std::vector<int> > H;       // For each node, all the assigned mesh vertices.
    int *VNA;                               // For each vertex, the assigned skeleton node.
    double *OX;                             // Output X coordinates.
    double *OY;                             // Output Y coordinates.
    double *OZ;                             // Output Z coordinates.
    
    if (isdebug) debug.open("debug.txt");
    
    if (isdebug) debug << "81" << std::endl;
    
    // Get inputs.
    MX = mxGetPr(prhs[0]);
    MY = mxGetPr(prhs[1]);
    MZ = mxGetPr(prhs[2]);
    SX = mxGetPr(prhs[3]);
    SY = mxGetPr(prhs[4]);
    SZ = mxGetPr(prhs[5]);
    MI = mxGetPr(prhs[6]);
    MJ = mxGetPr(prhs[7]);
    SI = mxGetPr(prhs[8]);
    SJ = mxGetPr(prhs[9]);
    nPnts = mxGetM(prhs[0]);
    nNodes = mxGetM(prhs[3]);
    nMEdges = mxGetM(prhs[6]);
    nSEdges = mxGetM(prhs[8]);
    
    if (isdebug) debug << "99" << std::endl;
    
    // Initialise node-vertex assignment storage.
    H.resize(nNodes);
    VNA = (int*)calloc(nPnts,sizeof(int));
    
    // Get node-vertex assignments from input.
    for (int i = 0; i < nNodes; i++)
    {
        mxArray *cell = mxGetCell(prhs[10],i);
        double *data = mxGetPr(cell);
        int n = mxGetM(cell);
        for (int j = 0; j < n; j++)
        {
            H[i].push_back((int)data[j]-1);
            VNA[(int)data[j]-1] = i;
        }
    }
    
    if (isdebug) debug << "118" << std::endl;
    
    // Initialise adjacency list of mesh edges.
    MEdges.resize(nPnts);
    
    // Compute adjacency list of mesh edges.
    for (int i = 0; i < nMEdges; i++)
        MEdges[MI[i]-1].push_back(MJ[i]-1);
    
    if (isdebug) debug << "127" << std::endl;
    
    // Initialise adjacency list of skeleton edges.
    SEdges.resize(nNodes);
    
    // Compute adjacency list of skeleton edges.
    for (int i = 0; i < nSEdges; i++)
        SEdges[SI[i]-1].push_back(SJ[i]-1);
    
    if (isdebug) debug << "136" << std::endl;
    
    // Initialise output.
    plhs[0] = mxCreateDoubleMatrix(nNodes, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nNodes, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nNodes, 1, mxREAL);
    OX = mxGetPr(plhs[0]);
    OY = mxGetPr(plhs[1]);
    OZ = mxGetPr(plhs[2]);
    
    if (isdebug) debug << "146" << std::endl;
    
    // Iterate through all skeleton nodes.
    for (int i = 0; i < nNodes; i++)
    {
        if (isdebug) debug << "Node = " << i+1 << std::endl;
        // Initialise node displacement.
        double d1=0, d2=0, d3=0;
        
        // If a leaf node, then just take the mean displacement of all
        // vertices assigned to the joint.
        if (SEdges[i].size() == 1)
        {
            for (int j = 0; j < H[i].size(); j++)
            {
                d1 += SX[i] - MX[H[i][j]];
                d2 += SY[i] - MY[H[i][j]];
                d3 += SZ[i] - MZ[H[i][j]];
            }
            OX[i] = SX[i] - (d1 / H[i].size());
            OY[i] = SY[i] - (d2 / H[i].size());
            OZ[i] = SZ[i] - (d3 / H[i].size());
            continue;
        }
        
        // Iterate though all neighbours.
        for (int j = 0; j < SEdges[i].size(); j++)
        {
            std::vector<int> boundary;
            double dj1=0, dj2=0, dj3=0, w=0;
            // Iterate through all vertices assigned to node.
            for (int k = 0; k < H[i].size(); k++)
            {
                int v = H[i][k];
                // If the vertex is on the boundary, add it to a list of
                // boundary vertices.
                for (int l = 0; l < MEdges[v].size(); l++)
                {
                    int n = MEdges[v][l];
                    
                    if (VNA[n] == SEdges[i][j])
                    {
                        boundary.push_back(v);
                        break;
                    }
                }
            }
            if (isdebug) debug << "Boundary size = " << boundary.size() << std::endl << std::endl;
            for (int k = 0; k < boundary.size(); k++)
            {
                int v = boundary[k];
                double wCurr = 0;
                for (int l = 0; l < MEdges[v].size(); l++)
                {
                    for (int m = 0; m < boundary.size(); m++)
                    {
                        if (MEdges[v][l] == boundary[m])
                        {
                            wCurr += sqrt(pow((MX[v] - MX[MEdges[v][l]]),2) +
                                pow((MY[v] - MY[MEdges[v][l]]),2) +
                                pow((MZ[v] - MZ[MEdges[v][l]]),2));
                            break;
                        }
                    }
                }
                w += wCurr;
                // Add displacement of current boundary vertex.
                dj1 += (SX[i] - MX[v]) * wCurr;
                dj2 += (SY[i] - MY[v]) * wCurr;
                dj3 += (SZ[i] - MZ[v]) * wCurr;
                            
            }
            d1 += dj1 / w;
            d2 += dj2 / w;
            d3 += dj3 / w;
        }
        OX[i] = SX[i] - (d1 / SEdges[i].size());
        OY[i] = SY[i] - (d2 / SEdges[i].size());
        OZ[i] = SZ[i] - (d3 / SEdges[i].size());
    }
    
    if (isdebug) debug.close();
}
