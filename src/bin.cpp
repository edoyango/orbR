#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

std::vector< std::vector< std::vector<int> > > bin(Rcpp::NumericMatrix& x, 
                                                   const double& dcell,
                                                   int ngridx[3],
                                                   double mingridx[3])
{

    int dim = x.nrow(), npoints = x.ncol();

    // Finding min-max extents
    // Initialization
    double maxgridx[3];
    for (int d = 0; d < dim; ++d)
    {
        mingridx[d] = x(d+1, 1);
        maxgridx[d] = x(d+1, 1);
    }
    if (dim == 2) 
    {
        mingridx[2] = 0.;
        maxgridx[2] = 0.;
    }
    // finding min-max
    for (int i = 0; i < npoints; ++i)
    {
        for (int d = 0; d < dim; ++d)
        {
            int i1 = i+1;
            int d1 = d+1;
            if (x(d1, i1) < mingridx[d]) mingridx[d] = x(d1, i1);
            if (x(d1, i1) > maxgridx[d]) maxgridx[d] = x(d1, i1);
        }
    }

    // number of grid cells
    for (int d = 0; d < 3; ++d)
    {
        ngridx[d] = static_cast<int>((maxgridx[d]-mingridx[d])/dcell) + 1;
    }

    // initializing and zeroing binned grid
    std::vector< std::vector< std::vector<int> > > \
    pincell(ngridx[0], std::vector< std::vector<int> >(ngridx[1], std::vector<int>(ngridx[2])));
    for (int i = 0; i < ngridx[0]; ++i)
    {
        for (int j = 0; j < ngridx[1]; ++j)
        {
            for (int k = 0; k < ngridx[2]; ++k)
            {
                pincell[i][j][k] = 0;
            }
        }
    }

    // counting points in array
    int icell[3];
    for (int i = 0; i < npoints; ++i)
    {
        for (int d = 0; d < dim; ++d)
        {
            icell[d] = static_cast<int>((x(d+1,i+1) - mingridx[d])/dcell);
        }
        if (dim==2) icell[2] = 0;
        pincell[icell[0]][icell[1]][icell[2]] ++;
    }

    return(pincell);
}