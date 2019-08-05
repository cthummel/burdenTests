#include <vector>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>




class skato
{
    public:
    skato(gsl_matrix* genoData, gsl_matrix* covariates, gsl_vector* variantFrequency, gsl_vector* phenotype);
    ~skato();
    double getPvalue(){return finalPvalue;}

    private:
    double generateQrho(double rho);
    void setupMatricies();
    double generateTestStatistic();

    gsl_matrix* data;
    gsl_matrix* cov;
    gsl_matrix* Vinverse;
    gsl_matrix* sqrtP;
    gsl_matrix* W;
    gsl_matrix* M;
    gsl_matrix* Z;

    gsl_vector* rhoGrid;
    gsl_vector* rhoStatistics;
    gsl_vector* rhoPvalues;
    gsl_vector* rhoQmin;
    gsl_vector* pheno;

    double finalPvalue;
};