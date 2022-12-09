#include<vector>
#include<iostream>
#include<cmath>
using namespace std;

void pincell_trim(const int gridind_in[2][3], 
                  std::vector< std::vector< std::vector<int> > > pincell,
                  int ngridx_trim[3]);

// function to recurisvely bisect pincell grid to obtain subdomain boundaries
void bisect(const int& dim, 
            const int& npoints_in,
            const int& npartitions,
            const int gridind_in[2][3],
            std::vector< std::vector< std::vector<int> > > pincell,
            const int& nprocs_in,
            const int& node_in,
            std::vector< std::vector< std::vector<int> > >& bounds,
            int& nleafreached)
{
    
    int ngridx_trim[3];
    pincell_trim(gridind_in, pincell, ngridx_trim);
    int A[3];
    A[0] = ngridx_trim[1]*ngridx_trim[2];
    A[1] = ngridx_trim[0]*ngridx_trim[2];
    A[2] = ngridx_trim[0]*ngridx_trim[1];
    int cax;
    if ( (A[0] <= A[1]) && (A[0] <= A[2]) )
    {
        cax = 0; // cutting plane orthogonal to x-axis
    } else if ( (A[1] <= A[0]) && (A[1] <= A[2]) )
    {
        cax = 1; // cutting plane orthogonal to y-axis
    } else 
    {
        cax = 2; // cutting plane orthogonal to z-axis
    }
    
    int i = gridind_in[0][cax], n_p = 0, pincol = 0;
    float nprocs_in_float = static_cast<float>(nprocs_in);
    int np_per_node = static_cast<int>(std::ceil(nprocs_in_float/2.)/nprocs_in_float*npoints_in);
    while (n_p < np_per_node)
    {

        pincol = 0;
        if (cax==0) // cutting plane orthogonal to x-axis
        {
            for (int j = gridind_in[0][1]; j < gridind_in[1][1]; ++j)
            {
                for (int k = gridind_in[0][2]; k < gridind_in[1][2]; ++k)
                {
                    pincol += pincell[i][j][k];
                }
            }
        } else if (cax==1) // cutting plane orthogonal to y-axis
        {
            for (int j = gridind_in[0][0]; j < gridind_in[1][0]; ++j)
            {
                for (int k = gridind_in[0][2]; k < gridind_in[1][2]; ++k)
                {
                    pincol += pincell[j][i][k];
                }
            }
        } else // cutting plane orthogonal to z-axis
        {
            for (int j = gridind_in[0][0]; j < gridind_in[1][0]; ++j)
            {
                for (int k = gridind_in[0][1]; k < gridind_in[1][1]; ++k)
                {
                    pincol += pincell[j][k][i];
                }
            }
        }
        i++;
        n_p += pincol;
    }
    if ((np_per_node - (n_p - pincol)) < (n_p - np_per_node))
    {
        i--;
        n_p -= pincol;
    }
    
    int gridind_out[2][3], node_out, npoints_out, nprocs_out;
    for (int i = 0; i < 2; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            gridind_out[i][d] = gridind_in[i][d];
        }
    }
    
    // lower first
    node_out = 2*node_in;
    npoints_out = n_p;
    gridind_out[0][cax] = gridind_in[0][cax];
    gridind_out[1][cax] = i;
    nprocs_out = std::ceil(static_cast<float>(nprocs_in)/2.);
    if (nprocs_out==1)
    {
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < dim; ++d)
            {
                bounds[nleafreached][i][d] = gridind_out[i][d];
            }
        }
        nleafreached++;
    } else
    {
        bisect(dim, npoints_out, npartitions, gridind_out, pincell, nprocs_out, node_out, bounds, nleafreached);
    }

    // upper next
    node_out = 2*node_in + 1;
    npoints_out = npoints_in - n_p;
    gridind_out[0][cax] = i;
    gridind_out[1][cax] = gridind_in[1][cax];
    nprocs_out = nprocs_in - nprocs_out; 
    if (nprocs_out==1)
    {
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < dim; ++d)
            {
                bounds[nleafreached][i][d] = gridind_out[i][d];
            }
        }
        nleafreached++;
    } else
    {
        bisect(dim, npoints_out, npartitions, gridind_out, pincell, nprocs_out, node_out, bounds, nleafreached);
    }
}

// helper function to trim the pincell grid to remove border 0 layers
// forwards declared because it's pretty verbose.
void pincell_trim(const int gridind_in[2][3], 
                  std::vector< std::vector< std::vector<int> > > pincell,
                  int ngridx_trim[3])
{
    int oldi[2], oldj[2], oldk[2], newi[2], newj[2], newk[2];
    for (int i = 0; i < 2; ++i)
    {
        oldi[i] = gridind_in[i][0];
        newi[i] = oldi[i];
        oldj[i] = gridind_in[i][1];
        newj[i] = oldj[i];
        oldk[i] = gridind_in[i][2];
        newk[i] = oldk[i];
    }

    // finding new start for i
    for (int i = oldi[0]; i < oldi[1]; ++i)
    {
        for (int j = oldj[0]; j < oldj[1]; ++j)
        {
            for (int k = oldk[0]; k < oldk[1]; ++k)
            {
                if (pincell[i][j][k] != 0) 
                {
                    newi[0] = i;
                    goto newendi;
                }
            }
        }
    }
    // finding new end for i
    newendi:
    for (int i = oldi[1]-1; i >= newi[0]; --i)
    {
        for (int j = oldj[0]; j < oldj[1]; ++j)
        {
            for (int k = oldk[0]; k < oldk[1]; ++k)
            {
                if (pincell[i][j][k] != 0) 
                {
                    newi[1] = i+1;
                    goto newstartj;
                }
            }
        }
    }
    // finding new start for j
    newstartj:
    for (int j = oldj[0]; j < oldj[1]; ++j)
    {
        for (int i = newi[0]; i < newi[1]; ++i)
        {
            for (int k = oldk[0]; k < oldk[1]; ++k)
            {
                if (pincell[i][j][k] != 0)
                {
                    newj[0] = j;
                    goto newendj;
                }
            }
        }
    }
    //finding new end for j
    newendj:
    for (int j = oldj[1]-1; j >= newj[0]; --j)
    {
        for (int i = newi[0]; i < newi[1]; ++i)
        {
            for (int k = oldk[0]; k < oldk[1]; ++k)
            {
                if (pincell[i][j][k] != 0)
                {
                    newj[1] = j+1;
                    goto newstartk;
                }
            }
        }
    }
    newstartk:
    for (int k = oldk[0]; k < oldk[1]; ++k)
    {
        for (int i = newi[0]; i < newi[1]; ++i)
        {
            for (int j = newj[0]; j < newj[1]; ++j)
            {
                if (pincell[i][j][k] != 0)
                {
                    newk[0] = k;
                    goto newendk;
                }
            }
        }
    }
    newendk:
    for (int k = oldk[1]-1; k >= newk[0]; --k)
    {
        for (int i = newi[0]; i < newi[1]; ++i)
        {
            for (int j = newj[0]; j < newj[1]; ++j)
            {
                if (pincell[i][j][k] != 0)
                {
                    newk[1] = k+1;
                    goto finishup;
                }
            }
        }
    }
    // calculating new number of grid cells in each dimension
    finishup:
    ngridx_trim[0] = newi[1] - newi[0];
    ngridx_trim[1] = newj[1] - newj[0];
    ngridx_trim[2] = newk[1] - newk[0];
}