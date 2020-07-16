#include <iostream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <list>
#include <cmath>
#include <limits>
#include <set>
// #include <matlab>
// Lib IGL includes
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/parallel_for.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <Eigen/SVD>


// Los que tienen el tipo de objeto definido son input, el resto output??

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // std::cout << "entered mex file" << std::endl;
    using namespace igl;
    using namespace igl::matlab;
   // using namespace igl::copyleft::cgal;
    using namespace Eigen;
    VectorXd AZ,U,X,C,weights;
    VectorXi II;
    double rho;
    
    
    igl::matlab::MexStream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    
    
    mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
    parse_rhs_double(prhs,AZ); // Aquí se pasa del prhs a la matriz V de Eigen
    parse_rhs_double(prhs+1,U); // Aquí se pasa del prhs a la matriz F de Eigen
    parse_rhs_index(prhs+2,II);
    rho = *mxGetPr(prhs[3]);
    parse_rhs_index(prhs+4,weights);
    C = AZ-U;
    
    

    X.resize(C.size());
//        MatrixXd MM(2,2);
//        MatrixXd Hval(2,2);
//        MatrixXd HH(2,2);
//        VectorXd S(2);
//        double s1,h1,s2,h2;
//    HH(1,0) = 0;
//    HH(0,1) = 0;
//
//    for (int i = 0; i < II.size(); i++) {
//        MM(0,0) = C(4*i);
//        MM(0,1) = C(4*i+1);
//        MM(1,0) = C(4*i+2);
//        MM(1,1) = C(4*i+3);
////
//        JacobiSVD<MatrixXd> svd( MM, ComputeFullV | ComputeFullU );
//        S = svd.singularValues();
//
//        if ((S(0)-(1/rho))>0) {
//            h1 = S(0)-(1/rho);
//                  }else{
//                      if((S(0)+(1/rho))<0){
//                          h1 = S(0)+(1/rho);
//                      }else{
//                          h1 = 0;
//                      }
//                  }
//        if ((S(1)-(1/rho))>0) {
//            h2 = S(1)-(1/rho);
//                  }else{
//                      if((S(1)+(1/rho))<0){
//                          h2 = S(1)+(1/rho);
//                                }else{
//                                    h2 = 0;
//                                }
//                        }
//        HH(0,0) = h1;
//        HH(1,1) = h2;
//
//        Hval = svd.matrixU()*HH*svd.matrixV().transpose();
//
//        X(4*i) = Hval(0,0);
//        X(4*i+1) = Hval(0,1);
//        X(4*i+2) = Hval(1,0);
//        X(4*i+3) = Hval(1,1);
//    }
//
    igl::parallel_for(C.size()/4,[&] (const int i){
        MatrixXd MM(2,2);
        MatrixXd Hval(2,2);
        MatrixXd HH(2,2);
        VectorXd S(2);
        double rho_with_weight;
        rho_with_weight = rho;
        //rho = rho_with_weight;
        double s1,h1,s2,h2;
        HH(1,0) = 0;
        HH(0,1) = 0;

        MM(0,0) = C(4*i);
        MM(0,1) = C(4*i+1);
        MM(1,0) = C(4*i+2);
        MM(1,1) = C(4*i+3);
        //
        JacobiSVD<MatrixXd> svd( MM, ComputeFullV | ComputeFullU );
        S = svd.singularValues();

        if ((S(0)-(1/rho_with_weight))>0) {
            h1 = S(0)-(1/rho_with_weight);
        }else{
            if((S(0)+(1/rho_with_weight))<0){
                h1 = S(0)+(1/rho_with_weight);
            }else{
                h1 = 0;
            }
        }
        if ((S(1)-(1/rho_with_weight))>0) {
            h2 = S(1)-(1/rho_with_weight);
        }else{
            if((S(1)+(1/rho_with_weight))<0){
                h2 = S(1)+(1/rho_with_weight);
            }else{
                h2 = 0;
            }
        }
        HH(0,0) = h1;
        HH(1,1) = h2;

        Hval = svd.matrixU()*HH*svd.matrixV().transpose();

        X(4*i) = Hval(0,0);
        X(4*i+1) = Hval(0,1);
        X(4*i+2) = Hval(1,0);
        X(4*i+3) = Hval(1,1);
    },0);
    

   // X = AZ;

    
    switch(nlhs)
    {
        case 1:
            prepare_lhs_double(X,plhs+0); // Sustituir W con las curvaturas?
        default:break;
    }
    
    // Restore the std stream buffer Important!
    std::cout.rdbuf(outbuf);
    return;
}
