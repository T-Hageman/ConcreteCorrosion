#include "IPScheme.h"

/// @brief Generates a standard Gauss integration scheme
/// @param wout output: Weights of integration points
/// @param coordsIPout output: Location of integration points, x = coords[ip](x/y)
/// @param ipCount input: Number of integration points per dimension
/// @param dim input: Dimension to generate scheme for
/// @return Total amount of integration points used in this scheme
int getIPScheme(std::vector<double> *wout, std::vector<Eigen::VectorXd> *coordsIPout, int ipCount, int dim){
    std::vector<double> w1D, C1D;
    getIPScheme1D(&w1D, &C1D, ipCount);

    int n_ip = std::pow(ipCount, dim);
    wout->resize(n_ip);
    coordsIPout->resize(n_ip);

    int cnt = -1;
    switch (dim){
        case 1:
            for (int i = 0; i < ipCount; i++){
                cnt=cnt+1;
                (*wout)[cnt] = w1D[i];
                (*coordsIPout)[cnt].resize(1);
                (*coordsIPout)[cnt](0) = C1D[i];
            }
            break;
        case 2:
            for (int j = 0; j < ipCount; j++){
                for (int i = 0; i < ipCount; i++){
                    cnt=cnt+1;
                    (*wout)[cnt] = w1D[i]*w1D[j];
                    (*coordsIPout)[cnt].resize(2);
                    (*coordsIPout)[cnt](0) = C1D[i];
                    (*coordsIPout)[cnt](1) = C1D[j];
                }
            }
            break;
        case 3:
            for (int k = 0; k < ipCount; k++){
                for (int j = 0; j < ipCount; j++){
                    for (int i = 0; i < ipCount; i++){
                        cnt=cnt+1;
                        (*wout)[cnt] = w1D[i]*w1D[j]*w1D[k];
                        (*coordsIPout)[cnt].resize(3);
                        (*coordsIPout)[cnt](0) = C1D[i];
                        (*coordsIPout)[cnt](1) = C1D[j];
                        (*coordsIPout)[cnt](2) = C1D[k];
                    }
                }
            }
            break;
    }

    return n_ip;
}


/// @brief Generates a one-dimensional Gauss integration scheme
/// @param wout output: integration weights
/// @param coordsIPout output: integration point coordinates
/// @param ipCount input: number of integration points
void getIPScheme1D(std::vector<double> *wout, std::vector<double> *coordsIPout, int ipCount){
    wout->resize(ipCount);
    coordsIPout->resize(ipCount);

    switch (ipCount) {
        case 1:
            coordsIPout->assign({0});
            wout->assign({2.0});
            break;
        case 2:
            coordsIPout->assign(-1/std::sqrt(3.0), 1/std::sqrt(3.0));
            wout->assign({1.0, 1.0});
            break;
        case 3:
            coordsIPout->assign({-std::sqrt(3.0/5.0), 0, std::sqrt(3.0/5.0)});
            wout->assign({5.0/9.0, 8.0/9.0, 5.0/9.0});
            break;
        case 4:
            coordsIPout->assign({-std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0)), 
                                 -std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)), 
                                 std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)), 
                                 std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0))});
            wout->assign({(18.0-std::sqrt(30.0))/36.0,
                          (18.0+std::sqrt(30.0))/36.0, 
                          (18.0+std::sqrt(30.0))/36.0, 
                          (18.0-std::sqrt(30.0))/36.0});
            break;
        case 5:
            coordsIPout->assign({-1.0/3.0*std::sqrt(5.0+2.0*std::sqrt(10.0/7.0)),
                                 -1.0/3.0*std::sqrt(5.0-2.0*std::sqrt(10.0/7.0)),
                                 0.0,
                                 1.0/3.0*std::sqrt(5.0-2.0*std::sqrt(10.0/7.0)),
                                 1.0/3.0*std::sqrt(5.0+2.0*std::sqrt(10.0/7.0))});
            wout->assign({(322.0-13.0*std::sqrt(70.0))/900.0,
                          (322.0+13.0*std::sqrt(70.0))/900.0,
                          128.0/225.0,
                          (322.0+13.0*std::sqrt(70.0))/900.0,
                          (322.0-13.0*std::sqrt(70.0))/900.0});
            break;
        default:
            throw std::invalid_argument("Higer order ip schemes not implemented in IPScheme.h");
            break;
    }
    for (size_t i = 0; i < wout->size(); i++){
        (*coordsIPout)[i] = ((*coordsIPout)[i]+1)/2;
        (*wout)[i] = (*wout)[i]/2;
    }
}

int getIPSchemeTriangle(std::vector<double> *wout, std::vector<Eigen::VectorXd> *coordsIPout, int ipCount, int dim){
    int n_ip = 0;

    if (dim==2){
        switch (ipCount) {
            case 1:
                n_ip = 1;
                wout->resize(n_ip);
                    (*wout)[0] = 0.5;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(2);
                    (*coordsIPout)[0](0)=1.0/3.0; (*coordsIPout)[0](1)=1.0/3.0; 
                break;
            case 2:
                n_ip = 3;
                wout->resize(n_ip);
                    (*wout)[0] = 1.0/6.0;
                    (*wout)[1] = 1.0/6.0;
                    (*wout)[2] = 1.0/6.0;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(2);
                    (*coordsIPout)[0](0)=1.0/6.0; (*coordsIPout)[0](1)=1.0/6.0; 

                    (*coordsIPout)[1].resize(2);
                    (*coordsIPout)[1](0)=4.0/6.0; (*coordsIPout)[1](1)=1.0/6.0; 

                    (*coordsIPout)[2].resize(2);
                    (*coordsIPout)[2](0)=1.0/6.0; (*coordsIPout)[2](1)=4.0/6.0;             
                break;
            case 3:
                n_ip = 4;
                wout->resize(n_ip);
                    (*wout)[0] = -27.0/96.0;
                    (*wout)[1] = 25.0/96.0;
                    (*wout)[2] = 25.0/96.0;
                    (*wout)[3] = 25.0/96.0;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(2);
                    (*coordsIPout)[0](0)=1.0/3.0; (*coordsIPout)[0](1)=1.0/3.0; 

                    (*coordsIPout)[1].resize(2);
                    (*coordsIPout)[1](0)=1.0/5.0; (*coordsIPout)[1](1)=1.0/5.0; 

                    (*coordsIPout)[2].resize(2);
                    (*coordsIPout)[2](0)=3.0/5.0; (*coordsIPout)[2](1)=1.0/5.0;      

                    (*coordsIPout)[3].resize(2);
                    (*coordsIPout)[3](0)=1.0/5.0; (*coordsIPout)[3](1)=3.0/5.0;    

                break;
            default:
                throw std::invalid_argument("Higer order ip schemes not implemented in IPScheme.h");
                break;
        }
    } else if (dim==3){ 
        if (true){ //https://www.cfd-online.com/Wiki/Code:_Quadrature_on_Tetrahedra
        switch (ipCount) {
            case 1:
                n_ip = 1;
                wout->resize(n_ip);
                    (*wout)[0] = 1;
                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=1./3.; (*coordsIPout)[0](1)=1./3.; (*coordsIPout)[0](2)=1./3.; 
                break;
            case 2:
                n_ip = 4;
                wout->resize(n_ip);
                    (*wout)[0] = 0.25;
                    (*wout)[1] = 0.25;
                    (*wout)[2] = 0.25;
                    (*wout)[3] = 0.25;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=0.5854101966249685; (*coordsIPout)[0](1)=0.1381966011250105; (*coordsIPout)[0](2)=0.1381966011250105; 

                (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=0.1381966011250105; (*coordsIPout)[1](1)=0.1381966011250105; (*coordsIPout)[1](2)=0.1381966011250105; 

                (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=0.1381966011250105; (*coordsIPout)[2](1)=0.1381966011250105; (*coordsIPout)[2](2)=0.5854101966249685; 

                (*coordsIPout)[3].resize(3);
                    (*coordsIPout)[3](0)=0.1381966011250105; (*coordsIPout)[3](1)=0.5854101966249685; (*coordsIPout)[3](2)=0.1381966011250105; 
                break;
            case 3:
                n_ip = 11;
                wout->resize(n_ip);
                    (*wout)[0] = -0.0789333333333333;
                    (*wout)[1] = 0.0457333333333333;
                    (*wout)[2] = 0.0457333333333333;
                    (*wout)[3] = 0.0457333333333333;
                    (*wout)[4] = 0.0457333333333333;
                    (*wout)[5] = 0.1493333333333333;
                    (*wout)[6] = 0.1493333333333333;
                    (*wout)[7] = 0.1493333333333333;
                    (*wout)[8] = 0.1493333333333333;
                    (*wout)[9] = 0.1493333333333333;
                    (*wout)[10] = 0.1493333333333333;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=0.2500000000000000; (*coordsIPout)[0](1)=0.2500000000000000; (*coordsIPout)[0](2)=0.2500000000000000; 

                (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=0.7857142857142857; (*coordsIPout)[1](1)=0.0714285714285714; (*coordsIPout)[1](2)=0.0714285714285714; 

                (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=0.0714285714285714; (*coordsIPout)[2](1)=0.0714285714285714; (*coordsIPout)[2](2)=0.0714285714285714; 

                (*coordsIPout)[3].resize(3);
                    (*coordsIPout)[3](0)=0.0714285714285714; (*coordsIPout)[3](1)=0.0714285714285714; (*coordsIPout)[3](2)=0.7857142857142857; 

                (*coordsIPout)[4].resize(3);
                    (*coordsIPout)[4](0)=0.0714285714285714; (*coordsIPout)[4](1)=0.7857142857142857; (*coordsIPout)[4](2)=0.0714285714285714; 

                (*coordsIPout)[5].resize(3);
                    (*coordsIPout)[5](0)=0.1005964238332008; (*coordsIPout)[5](1)=0.3994035761667992; (*coordsIPout)[5](2)=0.3994035761667992; 

                (*coordsIPout)[6].resize(3);
                    (*coordsIPout)[6](0)=0.3994035761667992; (*coordsIPout)[6](1)=0.1005964238332008; (*coordsIPout)[6](2)=0.3994035761667992; 

                (*coordsIPout)[7].resize(3);
                    (*coordsIPout)[7](0)=0.3994035761667992; (*coordsIPout)[7](1)=0.3994035761667992; (*coordsIPout)[7](2)=0.1005964238332008; 

                (*coordsIPout)[8].resize(3);
                    (*coordsIPout)[8](0)=0.3994035761667992; (*coordsIPout)[8](1)=0.1005964238332008; (*coordsIPout)[8](2)=0.1005964238332008; 

                (*coordsIPout)[9].resize(3);
                    (*coordsIPout)[9](0)=0.1005964238332008; (*coordsIPout)[9](1)=0.3994035761667992; (*coordsIPout)[9](2)=0.1005964238332008; 

                (*coordsIPout)[10].resize(3);
                    (*coordsIPout)[10](0)=0.1005964238332008; (*coordsIPout)[10](1)=0.1005964238332008; (*coordsIPout)[10](2)=0.3994035761667992; 

                break;
            case 4:
                n_ip = 15;
                wout->resize(n_ip);
                    (*wout)[0] = 0.1817020685825351;
                    (*wout)[1] = 0.0361607142857143;
                    (*wout)[2] = 0.0361607142857143;
                    (*wout)[3] = 0.0361607142857143;
                    (*wout)[4] = 0.0361607142857143;
                    (*wout)[5] = 0.0698714945161738;
                    (*wout)[6] = 0.0698714945161738;
                    (*wout)[7] = 0.0698714945161738;
                    (*wout)[8] = 0.0698714945161738;
                    (*wout)[9] = 0.0656948493683187;
                    (*wout)[10] = 0.0656948493683187;
                    (*wout)[11] = 0.0656948493683187;
                    (*wout)[12] = 0.0656948493683187;
                    (*wout)[13] = 0.0656948493683187;
                    (*wout)[14] = 0.0656948493683187;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=0.2500000000000000; (*coordsIPout)[0](1)=0.2500000000000000; (*coordsIPout)[0](2)=0.2500000000000000; 

                (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=0.0000000000000000; (*coordsIPout)[1](1)=0.3333333333333333; (*coordsIPout)[1](2)=0.3333333333333333; 

                (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=0.3333333333333333; (*coordsIPout)[2](1)=0.3333333333333333; (*coordsIPout)[2](2)=0.3333333333333333; 

                (*coordsIPout)[3].resize(3);
                    (*coordsIPout)[3](0)=0.3333333333333333; (*coordsIPout)[3](1)=0.3333333333333333; (*coordsIPout)[3](2)=0.0000000000000000; 

                (*coordsIPout)[4].resize(3);
                    (*coordsIPout)[4](0)=0.3333333333333333; (*coordsIPout)[4](1)=0.0000000000000000; (*coordsIPout)[4](2)=0.3333333333333333; 

                (*coordsIPout)[5].resize(3);
                    (*coordsIPout)[5](0)=0.7272727272727273; (*coordsIPout)[5](1)=0.0909090909090909; (*coordsIPout)[5](2)=0.0909090909090909; 

                (*coordsIPout)[6].resize(3);
                    (*coordsIPout)[6](0)=0.0909090909090909; (*coordsIPout)[6](1)=0.0909090909090909; (*coordsIPout)[6](2)=0.0909090909090909; 

                (*coordsIPout)[7].resize(3);
                    (*coordsIPout)[7](0)=0.0909090909090909; (*coordsIPout)[7](1)=0.0909090909090909; (*coordsIPout)[7](2)=0.7272727272727273; 

                (*coordsIPout)[8].resize(3);
                    (*coordsIPout)[8](0)=0.0909090909090909; (*coordsIPout)[8](1)=0.7272727272727273; (*coordsIPout)[8](2)=0.0909090909090909; 

                (*coordsIPout)[9].resize(3);
                    (*coordsIPout)[9](0)=0.4334498464263357; (*coordsIPout)[9](1)=0.0665501535736643; (*coordsIPout)[9](2)=0.0665501535736643; 

                (*coordsIPout)[10].resize(3);
                    (*coordsIPout)[10](0)=0.0665501535736643; (*coordsIPout)[10](1)=0.4334498464263357; (*coordsIPout)[10](2)=0.0665501535736643; 

                (*coordsIPout)[11].resize(3);
                    (*coordsIPout)[11](0)=0.0665501535736643; (*coordsIPout)[11](1)=0.0665501535736643; (*coordsIPout)[11](2)=0.4334498464263357; 

                (*coordsIPout)[12].resize(3);
                    (*coordsIPout)[12](0)=0.0665501535736643; (*coordsIPout)[12](1)=0.4334498464263357; (*coordsIPout)[12](2)=0.4334498464263357; 

                (*coordsIPout)[13].resize(3);
                    (*coordsIPout)[13](0)=0.4334498464263357; (*coordsIPout)[13](1)=0.0665501535736643; (*coordsIPout)[13](2)=0.4334498464263357; 

                (*coordsIPout)[14].resize(3);
                    (*coordsIPout)[14](0)=0.4334498464263357; (*coordsIPout)[14](1)=0.4334498464263357; (*coordsIPout)[14](2)=0.0665501535736643; 

                break;
            default:
                throw std::invalid_argument("Higer order ip schemes not implemented in IPScheme.h");
                break;
        }
        } else { //http://aero-comlab.stanford.edu/Papers/williams_jcam_2014.pdf
        switch (ipCount) {
            case 1:  
                n_ip = 1;
                wout->resize(n_ip);
                    (*wout)[0] = 1.;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=1.0/3.0; (*coordsIPout)[0](1)=1.0/3.0; (*coordsIPout)[0](2)=1.0/3.0; 
                break;
            case 2:
                n_ip = 3;
                wout->resize(n_ip);
                    (*wout)[0] = 1./3.;
                    (*wout)[1] = 1./3.;
                    (*wout)[2] = 1./3.;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=2.0/3.0; (*coordsIPout)[0](1)=1.0/3.0; (*coordsIPout)[0](2)=1.0/3.0; 

                (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=1.0/3.0; (*coordsIPout)[1](1)=2.0/3.0; (*coordsIPout)[1](2)=1.0/3.0; 

                (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=1.0/3.0; (*coordsIPout)[2](1)=1.0/3.0; (*coordsIPout)[2](2)=2.0/3.0; 
                break;
            case 3:
                n_ip = 6;
                wout->resize(n_ip);
                    (*wout)[0] = 0.109951743655333;
                    (*wout)[1] = 0.109951743655333;
                    (*wout)[2] = 0.109951743655333;
                    (*wout)[3] = 0.223381589678000;
                    (*wout)[4] = 0.223381589678000;
                    (*wout)[5] = 0.223381589678000;

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=0.816847572980440; (*coordsIPout)[0](1)=0.091576213509780; (*coordsIPout)[0](2)=0.091576213509780; 

                (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=0.091576213509780; (*coordsIPout)[1](1)=0.816847572980440; (*coordsIPout)[1](2)=0.091576213509780; 

                (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=0.091576213509780; (*coordsIPout)[2](1)=0.091576213509780; (*coordsIPout)[2](2)=0.816847572980440; 

                (*coordsIPout)[3].resize(3);
                    (*coordsIPout)[3](0)=0.108103018168071; (*coordsIPout)[3](1)=0.445948490915964; (*coordsIPout)[3](2)=0.445948490915964; 

                (*coordsIPout)[4].resize(3);
                    (*coordsIPout)[4](0)=0.445948490915964; (*coordsIPout)[4](1)=0.108103018168071; (*coordsIPout)[4](2)=0.445948490915964; 

                (*coordsIPout)[5].resize(3);
                    (*coordsIPout)[5](0)=0.445948490915964; (*coordsIPout)[5](1)=0.445948490915964; (*coordsIPout)[5](2)=0.108103018168071; 
                break;
            case 4:
                n_ip = 10;
                wout->resize(n_ip);
                coordsIPout->resize(n_ip);

                wout->resize(n_ip);
                    (*wout)[0] = 0.041955512996649;
                    (*wout)[1] = 0.041955512996649;
                    (*wout)[2] = 0.041955512996649;
                    (*wout)[3] = 0.112098412070887;
                    (*wout)[4] = 0.112098412070887;
                    (*wout)[5] = 0.112098412070887;       
                    (*wout)[6] = 0.112098412070887;
                    (*wout)[7] = 0.112098412070887;
                    (*wout)[8] = 0.112098412070887;
                    (*wout)[9] = 0.201542988584730;         

                coordsIPout->resize(n_ip);
                    (*coordsIPout)[0].resize(3);
                    (*coordsIPout)[0](0)=0.055564052669793; (*coordsIPout)[0](1)=0.055564052669793; (*coordsIPout)[0](2)=0.888871894660413; 

                    (*coordsIPout)[1].resize(3);
                    (*coordsIPout)[1](0)=0.055564052669793; (*coordsIPout)[1](1)=0.888871894660413; (*coordsIPout)[1](2)=0.055564052669793; 

                    (*coordsIPout)[2].resize(3);
                    (*coordsIPout)[2](0)=0.888871894660413; (*coordsIPout)[2](1)=0.055564052669793; (*coordsIPout)[2](2)=0.055564052669793; 

                    (*coordsIPout)[3].resize(3);
                    (*coordsIPout)[3](0)=0.295533711735893; (*coordsIPout)[3](1)=0.634210747745723; (*coordsIPout)[3](2)=0.070255540518384; 

                    (*coordsIPout)[4].resize(3);
                    (*coordsIPout)[4](0)=0.295533711735893; (*coordsIPout)[4](1)=0.070255540518384; (*coordsIPout)[4](2)=0.634210747745723; 

                    (*coordsIPout)[5].resize(3);
                    (*coordsIPout)[5](0)=0.070255540518384; (*coordsIPout)[5](1)=0.295533711735893; (*coordsIPout)[5](2)=0.634210747745723; 

                    (*coordsIPout)[6].resize(3);
                    (*coordsIPout)[6](0)=0.634210747745723; (*coordsIPout)[6](1)=0.295533711735893; (*coordsIPout)[6](2)=0.070255540518384; 

                    (*coordsIPout)[7].resize(3);
                    (*coordsIPout)[7](0)=0.634210747745723; (*coordsIPout)[7](1)=0.070255540518384; (*coordsIPout)[7](2)=0.295533711735893; 

                    (*coordsIPout)[8].resize(3);
                    (*coordsIPout)[8](0)=0.070255540518384; (*coordsIPout)[8](1)=0.634210747745723; (*coordsIPout)[8](2)=0.295533711735893; 

                    (*coordsIPout)[9].resize(3);
                    (*coordsIPout)[9](0)=0.333333333333333; (*coordsIPout)[9](1)=0.333333333333333; (*coordsIPout)[9](2)=0.333333333333333; 
                break;
            default:
                throw std::invalid_argument("Higer order ip schemes not implemented in IPScheme.h");
                break;
        }
        }
        for (int i = 0; i < n_ip; i++){
            (*wout)[i] = (*wout)[i]/6.0;
        }
    }

    return n_ip;
}