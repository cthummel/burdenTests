
#include "skato.hpp"
#include "davies.cpp"
#include "cdflib.cpp"
#include <chrono>
#include <gsl/gsl_integration.h>

struct QrhoInfo
{
    double muQ;
    double sigmaQ;
    double l;
    double delta;
    double a;
    double rho;
    gsl_vector *lambdas;

    ~QrhoInfo()
    {
        if (lambdas)
        {
            gsl_vector_free(lambdas);
        }
    }
};

struct pvalueInfo
{
    gsl_vector *Qmins;
    gsl_vector *tauRhos;
    gsl_vector *rhos;
    gsl_vector *lambdak;
    double muQ;
    double sigmaQ;
    double sigmaXi;

    ~pvalueInfo()
    {
        if (Qmins)
        {
            gsl_vector_free(Qmins);
        }
        if (tauRhos)
        {
            gsl_vector_free(tauRhos);
        }
        if (rhos)
        {
            gsl_vector_free(rhos);
        }
        if (lambdak)
        {
            gsl_vector_free(lambdak);
        }
    }
};

double pvalueIntegration(double x, void *p)
{
    struct pvalueInfo *params = (pvalueInfo *)p;
    double minValue = DBL_MAX;

    //Find minimum
    for (int i = 0; i < params->Qmins->size; i++)
    {
        double temp = (gsl_vector_get(params->Qmins, i) - gsl_vector_get(params->tauRhos, i) * x) / (1.0 - gsl_vector_get(params->rhos, i));
        if (temp < minValue)
        {
            minValue = temp;
        }
    }

    double tempStatistic = (minValue - params->muQ) * sqrt((params->sigmaQ * params->sigmaQ) - (params->sigmaXi * params->sigmaXi)) / params->sigmaQ;
    tempStatistic += params->muQ;

    int eigenCount = 0;
    int fault;
    double pvalue = 0.0;
    double sigma = 0.0;
    double lim = 10000;  //Max iterations
    double acc = .00001; //Target accuracy
    double trace[7];     //Error cache

    gsl_vector *noncentral = gsl_vector_calloc(params->lambdak->size);
    gsl_vector_int *df = gsl_vector_int_calloc(params->lambdak->size);
    gsl_vector_int_set_all(df, 1);

    for (int j = 0; j < params->lambdak->size; j++)
    {
        if (gsl_vector_get(params->lambdak, j) != 0)
        {
            eigenCount++;
        }
    }

    pvalue = 1 - qf(params->lambdak->data, noncentral->data, df->data, eigenCount,
                    sigma, tempStatistic, lim, acc, trace, &fault);

    //indefinite integral of F(delta(x)|lambdak)f(x|chisq(1))
    return 1.0 - (pvalue * gsl_ran_chisq_pdf(x, 1));
}

gsl_vector *logisticRegression(gsl_vector *pheno, gsl_matrix *cov)
{
    gsl_vector *result = gsl_vector_calloc(pheno->size);

    double mean = gsl_stats_mean(pheno->data, 1, pheno->size);
    gsl_vector_set_all(result, 1.0 / (1 + exp(-mean)));
    //Return the logit^-1
    return result;
}

gsl_matrix *matrixInverse(gsl_matrix *m)
{
    int signum;
    gsl_matrix *inverse = gsl_matrix_alloc(m->size1, m->size1);
    gsl_permutation *perm = gsl_permutation_alloc(m->size1);

    gsl_linalg_LU_decomp(m, perm, &signum);
    gsl_linalg_LU_invert(m, perm, inverse);

    gsl_permutation_free(perm);
    return inverse;
}

int sqrtMatrix(gsl_matrix *m)
{
    int errorCode = 0;
    gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(m->size1);
    gsl_matrix *result = gsl_matrix_alloc(m->size1, m->size1);
    gsl_vector *eval = gsl_vector_alloc(m->size1);
    gsl_matrix *evec = gsl_matrix_alloc(m->size1, m->size1);
    gsl_matrix *D = gsl_matrix_calloc(m->size1, m->size1);

    gsl_eigen_symmv(m, eval, evec, work);
    for (int i = 0; i < eval->size; i++)
    {
        gsl_matrix_set(D, i, i, sqrt(gsl_vector_get(eval, i)));
    }
    gsl_matrix *temp = gsl_matrix_alloc(m->size1, m->size1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, D, 0, temp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, matrixInverse(evec), 0, result);

    //m = result;
    //If that doesnt work use below.
    gsl_matrix_memcpy(m, result);

    gsl_eigen_symmv_free(work);
    gsl_matrix_free(temp);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(D);
    gsl_matrix_free(result);
    return errorCode;
}

skato::skato(gsl_matrix *genoData, gsl_matrix *covariates, gsl_vector *variantFrequency, gsl_vector *phenotype)
{
    //Save parameters.
    data = genoData;
    cov = covariates;
    pheno = phenotype;

    //Generate Weight Matrix.
    cout << "Generate Weight Matrix" << endl;
    W = gsl_matrix_calloc(variantFrequency->size, variantFrequency->size);
    cout << W->size1 << ", " << W->size2 << endl;
    for (int i = 0; i < variantFrequency->size; i++)
    {
        gsl_matrix_set(W, i, i, gsl_ran_beta_pdf(gsl_vector_get(variantFrequency, i), 1, 25));
    }

    //Generate rho grid
    cout << "Generate rho grid" << endl;
    rhoGrid = gsl_vector_calloc(11);
    rhoStatistics = gsl_vector_calloc(rhoGrid->size);
    rhoPvalues = gsl_vector_calloc(rhoGrid->size);
    rhoQmin = gsl_vector_calloc(rhoGrid->size);
    vector<QrhoInfo> rhoQminInfo(rhoGrid->size);
    //gsl_matrix* rhoQminInfo = gsl_matrix_calloc(rhoGrid->size, 3);
    for (int i = 0; i < rhoGrid->size; i++)
    {
        cout << ((rhoGrid->size - 1) / 100.0) * i << " ";
        gsl_vector_set(rhoGrid, i, ((rhoGrid->size - 1) / 100.0) * i);
    }
    cout << endl;

    //Assuming Binary Response Variable.
    cout << "Assuming Binary Response Variable" << endl;
    gsl_vector *uhat = logisticRegression(pheno, cov);
    Vinverse = gsl_matrix_calloc(uhat->size, uhat->size);
    for (int i = 0; i < uhat->size; i++)
    {
        gsl_matrix_set(Vinverse, i, i, gsl_vector_get(uhat, i) * (1 - gsl_vector_get(uhat, i)));
    }

    //Generate sqrtP
    cout << "Generate sqrtP" << endl;
    sqrtP = gsl_matrix_calloc(Vinverse->size1, Vinverse->size2);
    gsl_matrix_memcpy(sqrtP, Vinverse);
    gsl_matrix *VX = gsl_matrix_alloc(Vinverse->size1, cov->size2);           //VX
    gsl_matrix *XtransverseV = gsl_matrix_alloc(cov->size2, Vinverse->size2); //X'V
    gsl_matrix *inverse = gsl_matrix_alloc(cov->size2, cov->size2);           //(X'VX)^-1
    gsl_matrix *partRight = gsl_matrix_alloc(cov->size2, Vinverse->size2);    //(X'VX)^-1 * X'V
    gsl_matrix *sqrtPright = gsl_matrix_alloc(cov->size1, Vinverse->size2);   //VX * ((X'VX)^-1) * X'V

    //cout << "(" << Vinverse->size1 << ", " << Vinverse->size2 << ") * (" << cov->size1 << ", " << cov->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Vinverse, cov, 0.0, VX);
    //cout << "(" << cov->size2 << ", " << cov->size1 << ") * (" << Vinverse->size1 << ", " << Vinverse->size2 << ")" << endl;
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, cov, Vinverse, 0.0, XtransverseV);
    //cout << "(" << XtransverseV->size1 << ", " << XtransverseV->size2 << ") * (" << cov->size1 << ", " << cov->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, XtransverseV, cov, 0.0, inverse);
    //cout << gsl_matrix_get(inverse, 0, 0) << endl;
    inverse = matrixInverse(inverse);
    //cout << gsl_matrix_get(inverse, 0, 0) << endl;
    //cout << "(" << inverse->size1 << ", " << inverse->size2 << ") * (" << XtransverseV->size1 << ", " << XtransverseV->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverse, XtransverseV, 0.0, partRight);
    //cout << "(" << VX->size1 << ", " << VX->size2 << ") * (" << partRight->size1 << ", " << partRight->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, VX, partRight, 0.0, sqrtPright);
    //cout << "(" << sqrtPright->size1 << ", " << sqrtPright->size2 << ")" << endl;
    gsl_matrix_sub(sqrtP, sqrtPright);
    //sqrtMatrix(sqrtP);

    gsl_matrix_free(VX);
    gsl_matrix_free(XtransverseV);
    gsl_matrix_free(inverse);
    gsl_matrix_free(partRight);
    gsl_matrix_free(sqrtPright);

    //Setup for Qrho
    gsl_vector *diff = gsl_vector_calloc(uhat->size);
    //gsl_vector* left = gsl_vector_calloc(genoData->size2);
    //gsl_vector* right = gsl_vector_calloc(genoData->size2);

    //Generate phenotype difference vector.
    cout << "Generate phenotype difference vector" << endl;
    for (int i = 0; i < uhat->size; i++)
    {
        gsl_vector_set(diff, i, gsl_vector_get(pheno, i) - gsl_vector_get(uhat, i));
    }

    /*
    //Left side setup to minimize calculations during generation of Qrho. (y-uhat)'GW
    gsl_vector* tempLeft = gsl_vector_calloc(genoData->size2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, genoData, diff, 0.0, tempLeft); //xG = (G'x')'
    gsl_blas_dgemv(CblasNoTrans, 1.0, W, tempLeft, 0.0, left);
    gsl_vector_free(tempLeft);

    //Right side setup to minimize calculations during generation of Qrho. WG(y-uhat)
    gsl_vector* tempRight = gsl_vector_calloc(genoData->size2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, genoData, diff, 0.0, tempRight);
    gsl_blas_dgemv(CblasNoTrans, 1.0, W, tempRight, 0.0, right);
    gsl_vector_free(tempRight);
    */

    //G is coded as a "v by n" matrix so we need to transpose it to the right shape.
    gsl_matrix *GW = gsl_matrix_calloc(genoData->size2, W->size2);
    gsl_matrix *WG = gsl_matrix_calloc(W->size1, genoData->size2);
    //cout << "(" << genoData->size2 << ", " << genoData->size1 << ") * (" << W->size1 << ", " << W->size2 << ")" << endl;
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, genoData, W, 0.0, GW);
    //cout << "(" << W->size1 << ", " << W->size2 << ") * (" << genoData->size1 << ", " << genoData->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W, genoData, 0.0, WG);

    //Generate Qrho, eigenvalues, pvalues, and Qmin
    gsl_matrix *R = gsl_matrix_calloc(genoData->size1, genoData->size1);
    
    for (int i = 0; i < rhoGrid->size; i++)
    {
        auto runStartTime = std::chrono::high_resolution_clock::now();
        //Generate Krho
        cout << "Generate K_" << gsl_vector_get(rhoGrid, i) << endl;
        gsl_matrix_set_all(R, gsl_vector_get(rhoGrid, i));
        gsl_matrix_add_diagonal(R, 1.0 - gsl_vector_get(rhoGrid, i));
        gsl_matrix *leftKernel = gsl_matrix_calloc(genoData->size2, genoData->size1);
        gsl_matrix *Krho = gsl_matrix_calloc(genoData->size2, genoData->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, GW, R, 0.0, leftKernel);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, leftKernel, WG, 0.0, Krho);

        //Generate sqrtP * Krho * sqrtP
        cout << "Generate sqrtP * K_" << gsl_vector_get(rhoGrid, i) << " * sqrtP" << endl;
        gsl_matrix *leftmatrix = gsl_matrix_calloc(genoData->size2, genoData->size2);
        gsl_matrix *pvalueMatrix = gsl_matrix_calloc(genoData->size2, genoData->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sqrtP, Krho, 0.0, leftmatrix);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, leftmatrix, sqrtP, 0.0, pvalueMatrix);

        //Generate Qrho
        cout << "Generate Q_" << gsl_vector_get(rhoGrid, i) << endl;
        double result;
        gsl_vector *tempQrho = gsl_vector_calloc(genoData->size2);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Krho, diff, 0.0, tempQrho);
        gsl_blas_ddot(diff, tempQrho, &result);
        gsl_vector_set(rhoStatistics, i, result);

        //Generate pvalue
        cout << "Generate pvalue" << endl;
        int eigenCount = 0;
        int fault;
        double pvalue = 0.0;
        double sigma = 0.0;
        double lim = 10000;  //Max iterations
        double acc = .00001; //Target accuracy
        double trace[7];     //Error cache

        gsl_vector *eigenValues = gsl_vector_calloc(pvalueMatrix->size2);
        gsl_vector *noncentral = gsl_vector_calloc(eigenValues->size);
        gsl_vector_int *df = gsl_vector_int_calloc(eigenValues->size);
        gsl_eigen_symm_workspace *work = gsl_eigen_symm_alloc(pvalueMatrix->size2);

        gsl_vector_int_set_all(df, 1);
        gsl_eigen_symm(pvalueMatrix, eigenValues, work);

        for (int j = 0; j < eigenValues->size; j++)
        {
            if (gsl_vector_get(eigenValues, j) != 0)
            {
                eigenCount++;
            }
        }
        cout << "Eigen Count: " << eigenCount << endl;
        cout << "Q: " << result << endl;

        pvalue = 1 - qf(eigenValues->data, noncentral->data, df->data, eigenCount,
                        sigma, result, lim, acc, trace, &fault);
        cout << "pvalue: " << pvalue << endl;
        gsl_vector_set(rhoPvalues, i, pvalue);

        //Generate rhoQminInfo
        cout << "Generate QminInfo" << endl;
        vector<double> c(4, 0);
        for (int j = 0; j < eigenValues->size; j++)
        {
            double lambda = gsl_vector_get(eigenValues, j);
            double lambda2 = lambda * lambda;
            c[0] += lambda;
            c[1] += lambda2;
            c[2] += lambda2 * lambda;
            c[3] += lambda2 * lambda2;
        }
        double muQ = c[0];
        double sigmaQ = sqrt(2 * c[1]);
        double s1 = c[2] / (c[1] * sqrt(c[1]));
        double s2 = c[3] / (c[1] * c[1]);
        double a, delta, l;
        if (s1 * s1 > s2)
        {
            a = 1.0 / (s1 - sqrt((s1 * s1) - s2));
            delta = a * a * ((s1 * a) - 1);
            l = (a * a) - (2 * delta);
        }
        else
        {
            a = 1.0 / s1;
            delta = 0;
            l = a * a;
        }
        double muX = l + delta;
        double sigmaX = sqrt(2) * a;

        //Save relevent info for later.
        rhoQminInfo[i].a = a;
        rhoQminInfo[i].l = l;
        rhoQminInfo[i].delta = delta;
        rhoQminInfo[i].muQ = muQ;
        rhoQminInfo[i].sigmaQ = sigmaQ;
        rhoQminInfo[i].rho = gsl_vector_get(rhoGrid, i);
        rhoQminInfo[i].lambdas = eigenValues;

        //Cleanup
        auto runEndTime = std::chrono::high_resolution_clock::now();
        cout << "rho iteration took " << std::chrono::duration_cast<std::chrono::milliseconds>(runEndTime - runStartTime).count() / 60000.0 << " minutes." << endl;
        cout << endl;
        gsl_matrix_free(leftKernel);
        gsl_matrix_free(Krho);
        gsl_matrix_free(leftmatrix);
        gsl_matrix_free(pvalueMatrix);
        //gsl_vector_free(eigenValues);
        gsl_vector_free(noncentral);
        gsl_vector_int_free(df);
        gsl_vector_free(tempQrho);
        gsl_eigen_symm_free(work);
    }
    
    //gsl_vector_free(left);
    //gsl_vector_free(right);
    gsl_vector_free(diff);
    gsl_vector_free(uhat);
    gsl_matrix_free(R);
    gsl_matrix_free(WG);

    //Generate Z
    cout << "Generate Z"<< endl;
    Z = gsl_matrix_calloc(genoData->size2, genoData->size1);
    cout << "(" << genoData->size2 << ", " << genoData->size2 << ") * (" << genoData->size2 << ", " << W->size2 << ") -> (" << Z->size1 << ", " << Z->size2 << ")" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sqrtP, GW, 0.0, Z);
    gsl_matrix_free(GW);

    //Generate M
    cout << "Generate M"<< endl;
    M = gsl_matrix_calloc(genoData->size2, genoData->size2);
    gsl_vector *zbar = gsl_vector_calloc(genoData->size2);
    for (int j = 0; j < Z->size1; j++)
    {
        gsl_vector_set(zbar, j, gsl_stats_mean(gsl_matrix_row(Z, j).vector.data, 1, Z->size2));
    }
    double zbarPrimeZbar;
    gsl_blas_ddot(zbar, zbar, &zbarPrimeZbar); //z'z
    gsl_vector *tempzbarleft = gsl_vector_calloc(zbar->size);
    gsl_vector_memcpy(tempzbarleft, zbar);
    gsl_vector_scale(tempzbarleft, 1.0 / zbarPrimeZbar);
    gsl_matrix_view zbarMatrix = gsl_matrix_view_vector(zbar, 1, zbar->size);
    gsl_matrix_view zbarleftMatrix = gsl_matrix_view_vector(tempzbarleft, zbar->size, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &zbarleftMatrix.matrix, &zbarMatrix.matrix, 0.0, M);
    gsl_vector_free(tempzbarleft);

    //Generate tau(rho)
    cout << "Generate taus" << endl;
    gsl_vector *tauRho = gsl_vector_calloc(rhoGrid->size);
    double tauSumAccumulator = 0;
    for (int i = 0; i < genoData->size1; i++)
    {
        double temp;
        gsl_vector_view tempView = gsl_matrix_column(Z, i);
        gsl_blas_ddot(zbar, &tempView.vector, &temp);
        tauSumAccumulator += temp * temp;
    }
    tauSumAccumulator = tauSumAccumulator / zbarPrimeZbar;
    double leftTau = genoData->size1 * genoData->size1 * zbarPrimeZbar;
    for (int i = 0; i < rhoGrid->size; i++)
    {
        double right = (1 - gsl_vector_get(rhoGrid, i)) * tauSumAccumulator;
        gsl_vector_set(tauRho, i, (leftTau * gsl_vector_get(rhoGrid, i)) + right);
    }

    //Generate eigenvalues
    cout << "Generate lambdas "<< endl;
    gsl_matrix *center = gsl_matrix_calloc(M->size1, M->size2);
    gsl_matrix *eigenRight = gsl_matrix_calloc(Z->size1, Z->size2);
    gsl_matrix *eigenReady = gsl_matrix_calloc(Z->size2, Z->size2);
    gsl_vector *eigenValues = gsl_vector_calloc(Z->size2);
    gsl_matrix_set_identity(center);
    cout << "Building (I - M)" << endl;
    gsl_matrix_sub(center, M);
    cout << "Building (I - M)Z" << endl;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, center, Z, 0.0, eigenRight);
    cout << "Building Z'(I - M)Z" << endl;
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Z, eigenRight, 0.0, eigenReady);
    gsl_eigen_symm_workspace *work = gsl_eigen_symm_alloc(Z->size2);
    gsl_eigen_symm(eigenReady, eigenValues, work);
    gsl_matrix_free(center);
    gsl_matrix_free(eigenRight);
    //gsl_matrix_free(eigenReady); we need this for later
    gsl_eigen_symm_free(work);

    //Generate mu_Q
    cout << "Generate mu_Q"<< endl;
    double meanQ = 0;
    for (int i = 0; i < eigenValues->size; i++)
    {
        meanQ += gsl_vector_get(eigenValues, i);
    }

    //Generate sigmaXi
    cout << "Generate sigma Xi"<< endl;
    gsl_matrix *ZM = gsl_matrix_calloc(Z->size2, Z->size1);
    gsl_matrix *ZMZ = gsl_matrix_calloc(Z->size2, Z->size2);
    gsl_matrix *traceOn = gsl_matrix_calloc(Z->size2, Z->size2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Z, M, 0.0, ZM);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ZM, Z, 0.0, ZMZ);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ZMZ, eigenReady, 0.0, traceOn);
    double sigmaXi = 0;
    for (int i = 0; i < traceOn->size1; i++)
    {
        sigmaXi += gsl_matrix_get(traceOn, i, i);
    }
    sigmaXi = 2.0 * sqrt(sigmaXi);
    gsl_matrix_free(ZM);
    gsl_matrix_free(ZMZ);
    gsl_matrix_free(eigenReady);
    gsl_matrix_free(traceOn);

    //Generate sigmaQ
    double sigmaQ = 0;
    for (int i = 0; i < eigenValues->size; i++)
    {
        sigmaQ += gsl_vector_get(eigenValues, i) * gsl_vector_get(eigenValues, i);
    }
    sigmaQ *= 2;
    sigmaQ += sigmaXi * sigmaXi;
    sigmaQ = sqrt(sigmaQ);

    //Generate T
    double T = 1.0;
    int index = 0;
    for (int i = 0; i < rhoPvalues->size; i++)
    {
        if (gsl_vector_get(rhoPvalues, i) < T)
        {
            T = gsl_vector_get(rhoPvalues, i);
            index = i;
        }
    }

    //Generate Qmins
    for (int i = 0; i < rhoGrid->size; i++)
    {
        double P, Q, bound, tstar;
        int status;
        int which = 1;
        tstar = (T - rhoQminInfo[i].muQ) / rhoQminInfo[i].sigmaQ * sqrt(2) * rhoQminInfo[i].a;
        tstar += rhoQminInfo[i].l + rhoQminInfo[i].delta;
        cdfchn(&which, &P, &Q, &tstar, &rhoQminInfo[i].l, &rhoQminInfo[i].delta, &status, &bound);
        if (status == 0)
        {
            //Use Q here because we want upper tail based off paper.
            gsl_vector_set(rhoQmin, i, Q);
        }
        else if (status < 0)
        {
            cout << "When generating qmin for rho " << gsl_vector_get(rhoGrid, i) << " the parameter " << -status << " is out of range." << endl;
            if(-status == 1){cout << "Parameter 1 = " << &P << " when range must be: [0, 1.0-1.0D-16)" << endl;}
            if(-status == 2){cout << "Parameter 2 = " << &Q << " when range must be: [0, 1.0-1.0D-16)" << endl;}
            if(-status == 3){cout << "Parameter 3 = " << &tstar << " when range must be: [0,1.0D+300]" << endl;}
            if(-status == 4){cout << "Parameter 4 = " << &rhoQminInfo[i].l << " when range must be: [ 1.0D-300, 1.0D+300]" << endl;}
            if(-status == 5){cout << "Parameter 5 = " << &rhoQminInfo[i].delta << " when range must be: [0,1.0D+4]" << endl;}
        }
        else if (status == 1)
        {
            cout << "When generating qmin for rho " << gsl_vector_get(rhoGrid, i) << "\"The answer appears to be lower than the lowest search bound\"" << endl;
        }
        else if (status == 2)
        {
            cout << "When generating qmin for rho " << gsl_vector_get(rhoGrid, i) << "\"The answer appears to be higher than the greatest search bound.\"" << endl;
        }
    }

    gsl_function integrationFunction;
    pvalueInfo parameters;
    parameters.Qmins = rhoQmin;
    parameters.rhos = rhoGrid;
    parameters.tauRhos = tauRho;
    parameters.lambdak = eigenValues;
    parameters.muQ = meanQ;
    parameters.sigmaQ = sigmaQ;
    parameters.sigmaXi = sigmaXi;

    integrationFunction.function = &pvalueIntegration;
    integrationFunction.params = &parameters;

    double pvalue, abserr;
    double epsabs = 0.0;
    double epsrel = 0.000001;
    size_t limit = 1000;

    gsl_integration_workspace *integrationWork = gsl_integration_workspace_alloc(limit);
    gsl_integration_qagi(&integrationFunction, epsabs, epsrel, limit, integrationWork, &pvalue, &abserr);

    finalPvalue = pvalue;
    cout << "Pvalue is " << pvalue << endl;

    //Final Cleanup of remaining gsl data structures.
    gsl_integration_workspace_free(integrationWork);
}

skato::~skato()
{
    if (Vinverse)
    {
        gsl_matrix_free(Vinverse);
    }
    if (sqrtP)
    {
        gsl_matrix_free(sqrtP);
    }
    if (W)
    {
        gsl_matrix_free(W);
    }
    if (M)
    {
        gsl_matrix_free(M);
    }
    if (Z)
    {
        gsl_matrix_free(Z);
    }
    if (rhoGrid)
    {
        gsl_vector_free(rhoGrid);
    }
    if (rhoPvalues)
    {
        gsl_vector_free(rhoPvalues);
    }
    if (rhoQmin)
    {
        gsl_vector_free(rhoQmin);
    }
    if (rhoStatistics)
    {
        gsl_vector_free(rhoStatistics);
    }
}