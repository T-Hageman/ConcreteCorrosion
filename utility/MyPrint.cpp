#include "utility.h"

#include "../InputsOutputs/inputData.h"
#include <petsc.h>

MyLogger Logs;

MyLogger::MyLogger(){
    info_Level = 99;
}

MyLogger::~MyLogger(){

}

void MyLogger::Init(inputData* inputs){
    bool infoset = inputs->GetOptional(info_Level,{"Logs","InfoLevel"});
    if (infoset==false){
        info_Level = 2;
    }
}

/// @brief Appends index of cpu core and prints string from every core
/// @param outString 
/// @param infoLevel 
void MyLogger::PrintEvery(std::string outString, uint infoLevel){
    if (infoLevel<=info_Level){
        PetscMPIInt    rank;
        PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
        for (size_t i = 0; i < infoLevel; i++){
            outString = "\t" + outString;
        }
        outString = std::to_string(rank) +"\t"+ outString;
        PetscCallThrow(PetscPrintf(PETSC_COMM_SELF,"%s",outString.c_str()));
    }
}

/// @brief Appends index of cpu core and prints string from core 0
/// @param outString 
/// @param infoLevel 
void MyLogger::PrintSingle(std::string outString, uint infoLevel){
    if (infoLevel<=info_Level){
        PetscMPIInt    rank;
        PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
        if (rank == 0){
            for (size_t i = 0; i < infoLevel; i++){
                outString = "\t" + outString;
            }
            outString = std::to_string(rank) +"\t"+ outString;
            PetscCallThrow(PetscPrintf(PETSC_COMM_SELF,"%s",outString.c_str()));
        }
    }
}

/// @brief Appends index of cpu core and prints string from every core
/// @param outStringStream 
/// @param infoLevel 
void MyLogger::PrintEvery(std::ostringstream& outStringStream, uint infoLevel){
    std::string outString = outStringStream.str();
    PrintEvery(outString, infoLevel);
}

/// @brief Appends index of cpu core and prints string from core 0
/// @param outStringStream 
/// @param infoLevel 
void MyLogger::PrintSingle(std::ostringstream& outStringStream, uint infoLevel){
    std::string outString = outStringStream.str();
    PrintSingle(outString, infoLevel);
}