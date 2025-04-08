#ifndef UTILITY_H
#define UTILITY_H

#define EIGEN_NO_AUTOMATIC_RESIZING

#include <string>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <petscvec.h>
#include <petsc.h>
#include <numeric>
#include <cmath>
#include <map>
#include <stdint.h>
#include <limits.h>
#include <Eigen/Dense>

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 7, 7> Matrix7d;
typedef Eigen::Matrix<double, 7, 1> Vector7d;

typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 18, 18> Matrix18d;
typedef Eigen::Matrix<double, 18, 1> Vector18d;
typedef Eigen::Matrix<double, 20, 20> Matrix20d;
typedef Eigen::Matrix<double, 20, 1> Vector20d;

void KronProd(Eigen::VectorXd* result, Eigen::VectorXd* in1, Eigen::VectorXd* in2);
void KronProd(Eigen::VectorXd& result, Eigen::VectorXd& in1, Eigen::VectorXd& in2);
void KronProd(Eigen::VectorXd& result, Eigen::VectorXd& in1, Eigen::VectorXd& in2, Eigen::VectorXd& in3);
double min(std::vector<std::vector<double>> &Matrix);
double max(std::vector<std::vector<double>> &Matrix);

void MatAdd(Mat& K, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols, Eigen::MatrixXd &K_el);
void MatAdd(Mat& K, PetscInt &rows, std::vector<PetscInt> &cols, Eigen::MatrixXd &K_el);
void MatAdd(Mat& K, std::vector<PetscInt> &rows, PetscInt &cols, Eigen::MatrixXd &K_el);
void MatAdd(Mat& K, PetscInt &rows, PetscInt &cols, double &K_el);

void VecAdd(Vec& V, std::vector<PetscInt> &rows, Eigen::VectorXd &F_el);
void VecAdd(Vec& V, PetscInt &rows, double &F_el);

double sgn(double x);
double sum(std::vector<double>& ToSum);
size_t sum(std::vector<size_t>& ToSum);
void GetLog10(std::vector<std::vector<double>> &NodeData);
size_t Find(std::vector<std::string> array, std::string ToFind, bool& found);

void ConstructBMat(Eigen::MatrixXd& B, Eigen::MatrixXd& G, int dim);

class inputData;

/// @brief Class for printing info to cmd (still to add, printing to log file)
class MyLogger {
   public:
      MyLogger();
      ~MyLogger();

      void Init(inputData* inputs);

      void PrintEvery(std::string outString, uint infoLevel);
      void PrintSingle(std::string outString, uint infoLevel);
      void PrintEvery(std::ostringstream& outString, uint infoLevel);
      void PrintSingle(std::ostringstream& outString, uint infoLevel);

      uint info_Level;

   private:
      std::string logName;
};

extern MyLogger Logs;

/// @brief GhostVector which allows the use of global indices
class GhostVector {
    public:
        GhostVector();
        ~GhostVector();

        void SetSizeSequential(size_t nGlobal);
        void SetSize(size_t nLocal, size_t nGlobal, std::vector<PetscInt> &GhostLocs);
        void Set(std::vector<PetscInt> &locs, std::vector<double> &vals, InsertMode operation);
        void Set(PetscInt locs, double vals, InsertMode operation);
        void Set(GhostVector& ToCopy);
        void GetValues(std::vector<size_t> &locs, Eigen::VectorXd &vals_out);
        void GetValues(std::vector<size_t> &locs, std::vector<double> &vals_out);
        void GetValues(std::vector<PetscInt> &locs, Eigen::VectorXd &vals_out);
        double GetValue(PetscInt locs);
        void Zero();

        void AssemblyStart();
        void AssemblyEnd();
        void SyncStart(InsertMode op);
        void SyncEnd(InsertMode op);

        Vec DataVector; //Vector which contains the actual data for this ghostvector, local values first then appended with ghost values from other cores
        std::map<size_t, size_t> GlobalToLocal; //Mapping from glabal numbering to local numbering
        std::map<size_t, size_t> LocalToGlobal; //Mapping from local numbering to global numbering
        std::vector<PetscInt> Ghosts;  //Locations of ghost points that are owned by other cores
    protected:

    private:

};

#endif