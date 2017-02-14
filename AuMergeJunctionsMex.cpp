/*
 * AuMergeJunctionsMex.cpp
 * Usage: AuMergeJunctionsMex(mesh.X,mesh.Y,mesh.Z,Skel.X,Skel.Y,Skel.Z,SI,SJ,Skel.H);
 * Refines the embedding of a skeleton using the method by merging
 * junctions as described by Au et al. SIGGRAPH 2008.
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
    if (nrhs != 9)
        mexErrMsgTxt("Nine input arguments required.");
    if (nlhs != 2)
        mexErrMsgTxt("Two output argument required.");
    
    // Initialise variables.
    double *MX;                             // Mesh vertices X coordinates.
    double *MY;                             // Mesh vertices Y coordinates.
    double *MZ;                             // Mesh vertices Z coordinates.
    double *SX;                             // Skeleton nodes X coordinates.
    double *SY;                             // Skeleton nodes Y coordinates.
    double *SZ;                             // Skeleton nodes Z coordinates.
    double *SI;                             // i index of non-zeros skeleton adjacency matrix entries.
    double *SJ;                             // j index of non-zeros skeleton adjacency matrix entries.
    std::vector<std::vector<int> > SEdges;  // Edges of skeleton (adjacency list).
    int nNodes;                             // Number of skeleton nodes.
    int nPnts;                              // Number of mesh vertices.
    int nSEdges;                            // Number of skeleton edges (bones).
    std::vector<std::vector<int> > H;       // For each node, all the assigned mesh vertices.
    double *OSEdges;                        // Output skeleton edges.
    std::vector<int> Jncs;                  // Junction nodes.
    double *NodeEnabled;                    // For each skeleton node, is it enabled (does it exist?)?
    
    if (isdebug) debug.open("debug.txt");
    
    // Get inputs.
    MX = mxGetPr(prhs[0]);
    MY = mxGetPr(prhs[1]);
    MZ = mxGetPr(prhs[2]);
    SX = mxGetPr(prhs[3]);
    SY = mxGetPr(prhs[4]);
    SZ = mxGetPr(prhs[5]);
    SI = mxGetPr(prhs[6]);
    SJ = mxGetPr(prhs[7]);
    nPnts = mxGetM(prhs[0]);
    nNodes = mxGetM(prhs[3]);
    nSEdges = mxGetM(prhs[6]);
    
    // Initialise node-vertex assignment storage.
    H.resize(nNodes);
    
    // Get node-vertex assignments from input.
    for (int i = 0; i < nNodes; i++)
    {
        mxArray *cell = mxGetCell(prhs[8],i);
        double *data = mxGetPr(cell);
        int n = mxGetM(cell);
        for (int j = 0; j < n; j++)
        {
            H[i].push_back((int)data[j]-1);
        }
    }
    
    // Initialise adjacency list of skeleton edges.
    SEdges.resize(nNodes);
    
    // Compute adjacency list of skeleton edges.
    for (int i = 0; i < nSEdges; i++)
        SEdges[SI[i]-1].push_back(SJ[i]-1);
    
    // Create list of junction nodes and enable all nodes.
    NodeEnabled = (double*)calloc(nNodes,sizeof(double));
    for (int i = 0; i < SEdges.size(); i++)
    {
        NodeEnabled[i] = 1;
        if (SEdges[i].size() > 2)
        {
            if (isdebug) debug << "Node " << i+1 << " is a junction" << std::endl;
            Jncs.push_back(i);
        }
    }
    if (isdebug) debug << std::endl;
    
    // Iterate through all junction nodes.
    for (int i = 0; i < Jncs.size(); i++)
    {
        // Compute the centeredness of the current junction node.
        int v = Jncs[i];
        double sigma = centeredness(SX[v], SY[v], SZ[v], MX, MY, MZ, H[v]);
        
        double sigmaj[SEdges[v].size()];
        
        // Iterate through all the neighbours of the junction node.
        for (int j = 0; j < SEdges[v].size(); j++)
        {
            int n = SEdges[v][j];
            
            std::vector<int> Htmp;
            for (int k = 0; k < H[v].size(); k++)
                Htmp.push_back(H[v][k]);
            for (int k = 0; k < H[n].size(); k++)
                Htmp.push_back(H[n][k]);
            
            sigmaj[j] = centeredness(SX[n], SY[n], SZ[n], MX, MY, MZ, Htmp);
        }
        
        // Find neighbour with the potential merger with the best
        // centeredness.
        double smallSigma = pow(10,10);
        int n = -1;
        for (int j = 0; j < SEdges[v].size(); j++)
        {
            if (sigmaj[j] < smallSigma)
            {
                smallSigma = sigmaj[j];
                n = SEdges[v][j];
            }
        }
        
        // If the best potential merger is better centered than the
        // original junction, then merge.
        if (smallSigma < (0.9*sigma))
        {
            // Assign all vertices assigned to current joint to new merged
            // joint.
            for (int j = 0; j < H[v].size(); j++)
                H[n].push_back(H[v][j]);
            
            // Update the neighbours of the merged junction.
            for (int j = 0; j < SEdges[v].size(); j++)
            {
                int nj = SEdges[v][j];
                if (isdebug) debug << "Node " << nj+1 << " is a neighbour of a deleted node" << std::endl;
                
                // If neighbour is equal to merged junction then skip.
                if (nj == n)
                {
                    for (int k = 0; k < SEdges[nj].size(); k++)
                        if (SEdges[nj][k] == v)
                            SEdges[nj].erase(SEdges[nj].begin()+k);
                    continue;
                }
                
                // Add neighbour from original junction to merged junction.
                SEdges[n].push_back(nj);
                if (isdebug) debug << "Node " << n+1 << " is now a neightbour of node " << nj+1 << std::endl;
                
                // Change references to original junction to references to
                // merged junction.
                for (int k = 0; k < SEdges[nj].size(); k++)
                    if (SEdges[nj][k] == v)
                        SEdges[nj][k] = n;
                if (isdebug) debug << "Node " << nj+1 << " is now a neightbour of node " << n+1 << std::endl;
            }
            
            // Set merged junction as disabled.
            NodeEnabled[v] = 0;
            if (isdebug) debug << "Node " << v+1 << " merged with node " << n+1 << std::endl << std::endl;
            
            Jncs[i] = n;
            i--;
        }
    }
    
    // Initialise output.
    nSEdges = 0;
    for (int i = 0; i < SEdges.size(); i++)
        if (NodeEnabled[i] == 1)
            for (int j = 0; j < SEdges[i].size(); j++)
                nSEdges++;
    plhs[0] = mxCreateDoubleMatrix(nSEdges, 2, mxREAL);
    OSEdges = mxGetPr(plhs[0]);
    
    // Populate output skeleton edges.
    int idx = 0;
    for (int i = 0; i < SEdges.size(); i++)
    {
        if (NodeEnabled[i] == 1)
        {
            for (int j = 0; j < SEdges[i].size(); j++)
            {
                OSEdges[POS(idx,0,nSEdges)] = i+1;
                OSEdges[POS(idx,1,nSEdges)] = SEdges[i][j]+1;
                idx ++;
            }
        }
    }
    
    // Initialise cell array for output history.
    plhs[1] = mxCreateCellMatrix((mwSize)H.size(), 1);
    for (int i = 0; i < H.size(); i++)
    {
        mxArray *matlabArray = mxCreateDoubleMatrix(H[i].size(),1,mxREAL);
        double *output = mxGetPr(matlabArray);
        for (int j = 0; j < H[i].size(); j++)
            output[j] = H[i][j]+1;
        mxSetCell(plhs[1], (mwIndex)i, matlabArray);
    }
    if (isdebug) debug.close();
}
