#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

std::vector< std::vector< std::vector<int> > > bin(Rcpp::NumericMatrix& x, 
                                                   const double& dcell,
                                                   int ngridx[3],
                                                   double mingridx[3]);
void bisect(const int& dim, 
            const int& npoints_in,
            const int& npartitions,
            const int gridind_in[2][3],
            std::vector< std::vector< std::vector<int> > > pincell,
            const int& nprocs_in,
            const int& node_in,
            std::vector< std::vector< std::vector<int> > >& bounds,
            int& nleafreached);

//' Find pairs of points within fixed cutoff
//'
//' @param x Matrix of positions where a column corrresponds to a point
//' @param npartitions Number of partitions to generate
//' @param dcell Side-length if square grid cell used to bin points
//'
//' @return List of bounadries
//' @export
//'
//' @examples
//' x <- matrix(rnorm(5000), nrow=2)
//' orb(x, npartitions=8, dcell=0.15)
// [[Rcpp::export]]
std::vector< std::vector< std::vector<double> > > orb(Rcpp::NumericMatrix& x,
                                                   const int& npartitions,
                                                   const double& dcell)
{
    // Retrieving number of geometric dimensions (rows) and number of points (columns)
    int dim = x.nrow(), npoints = x.ncol();

    // Throwing error message if not in 2 or 3 spatial dimensions
    if ( (dim != 2) && (dim != 3)) Rcpp::stop("Input matrix rows should be either 2 or 3!");

    int ngridx[3];
    double mingridx[3];
    std::vector< std::vector< std::vector<int> > > pincell = bin(x, dcell, ngridx, mingridx);

    int gridind0[2][3], nleafreached = 0;
    for (int i = 0; i < 3; ++i)
    {
        gridind0[0][i] = 0;
        gridind0[1][i] = ngridx[i];
    }
    std::vector< std::vector< std::vector<int> > > \
    bounds_ind(npartitions,std::vector< std::vector<int> >(2, vector<int>(dim)));
    bisect(dim, npoints, npartitions, gridind0, pincell, npartitions, 1, bounds_ind, nleafreached);

    std::vector< std::vector< std::vector<double> > > \
    bounds(npartitions,std::vector< std::vector<double> >(2, vector<double>(dim)));
    for (int i = 0; i < npartitions; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int d = 0; d < dim; ++d)
            {
                bounds[i][j][d] = mingridx[d] + (bounds_ind[i][j][d]+j)*dcell;
            }
        }
    }

    return(bounds);
}