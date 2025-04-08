#include "inputData.h"
#include "../utility/utility.h"

inputData::inputData(){
    
}

/// @brief parses input data provided in json formatted file filename
/// @param filename input file name (string)
inputData::inputData(char const* filename){
    parse_file( filename );
    print();
}

inputData::~inputData(){

}

/// @brief Reads the input file
/// @param filename 
void inputData::parse_file( char const* filename ){
    std::ifstream ifs(filename);
    rapidjson::IStreamWrapper isw(ifs);
    
    inputs.ParseStream(isw);
} 

/// @brief Prints the input file contents to the terminal/output log
void inputData::print(){
    PetscMPIInt    rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (rank == 0){
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        inputs.Accept(writer);
        std::string inputPrint(buffer.GetString());
        Logs.PrintSingle("input file:\n"+inputPrint + "\n",1);        
    }
}

void inputData::KeyErrorMessage(std::vector<std::string> inputKey){
    std::string msg = "Error: InputFile does not have a key for ";
    for (size_t i = 0; i < inputKey.size(); i++){
        msg += inputKey[i] + ".";
    }
    Logs.PrintEvery(msg, 0);
    throw std::invalid_argument(msg);
}

void inputData::TypeErrorMessage(std::vector<std::string> inputKey, std::string typestring){
    std::string msg = "Error: InputFile key ";
    for (size_t i = 0; i < inputKey.size(); i++){
        msg += inputKey[i] + ".";
    }
    msg += "needs to be of type ";
    msg += typestring;
    Logs.PrintEvery(msg, 0);
    throw std::invalid_argument(msg);    
}

std::string inputData::GetType(std::vector<std::string> inputKey){
    bool exists;
    rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
    std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
    inputIter[0].CopyFrom(inputs, a);

    for (size_t i = 0; i < inputKey.size(); i++){
        exists = inputIter[i].HasMember(inputKey[i].c_str());
        if (exists){
            inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
        } else {
            return "None";
        }
    }

    rapidjson::Type dataType = inputIter[inputKey.size()].GetType();
    if (dataType == rapidjson::kStringType){
        return "String";
    } else if (dataType == rapidjson::kNumberType){
        return "Number";
    } else {
        return "Other";
    }

}

void inputData::GetRequired(std::string &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::string &data, std::vector<std::string> inputKey){
    bool exists;
    try {
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetString();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(double &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(double &data, std::vector<std::string> inputKey){
    bool exists;
    try {
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetDouble();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(bool &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(bool &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetBool();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<bool> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<bool> &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        size_t ArraySize = inputIter[inputKey.size()].Size();
        data.resize(ArraySize);
        for (size_t i = 0; i < ArraySize; i++){
            data[i] = inputIter[inputKey.size()][i].GetBool();
        }
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(size_t &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(size_t &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetUint64();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(uint &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(uint &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetUint();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(int &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(int &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        data = inputIter[inputKey.size()].GetInt();
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<std::string> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<std::string> &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        size_t ArraySize = inputIter[inputKey.size()].Size();
        data.resize(ArraySize);
        for (size_t i = 0; i < ArraySize; i++){
            data[i] = inputIter[inputKey.size()][i].GetString();
        }
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<int> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<int> &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        size_t ArraySize = inputIter[inputKey.size()].Size();
        data.resize(ArraySize);
        for (size_t i = 0; i < ArraySize; i++){
            data[i] = inputIter[inputKey.size()][i].GetInt();
        }
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<double> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<double> &data, std::vector<std::string> inputKey){
    bool exists;
    try{
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        size_t ArraySize = inputIter[inputKey.size()].Size();
        data.resize(ArraySize);
        for (size_t i = 0; i < ArraySize; i++){
            data[i] = inputIter[inputKey.size()][i].GetDouble();
        }
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<size_t> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<size_t> &data, std::vector<std::string> inputKey){
    bool exists;
    try{
    rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
    std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
    inputIter[0].CopyFrom(inputs, a);

    for (size_t i = 0; i < inputKey.size(); i++){
        exists = inputIter[i].HasMember(inputKey[i].c_str());
        if (exists){
            inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
        } else {
            return false;
        }
    }
    size_t ArraySize = inputIter[inputKey.size()].Size();
    data.resize(ArraySize);
    for (size_t i = 0; i < ArraySize; i++){
        data[i] = inputIter[inputKey.size()][i].GetUint64();
    }
    return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

void inputData::GetRequired(std::vector<std::vector<std::string>> &data, std::vector<std::string> inputKey){
    bool exists = GetOptional(data, inputKey);
    if (exists == false){
        KeyErrorMessage(inputKey);
    }
}

bool inputData::GetOptional(std::vector<std::vector<std::string>> &data, std::vector<std::string> inputKey){
    bool exists;
    try {
        rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
        std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
        inputIter[0].CopyFrom(inputs, a);

        for (size_t i = 0; i < inputKey.size(); i++){
            exists = inputIter[i].HasMember(inputKey[i].c_str());
            if (exists){
                inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
            } else {
                return false;
            }
        }
        size_t ArraySize = inputIter[inputKey.size()].Size();
        data.resize(ArraySize);
        for (size_t i = 0; i < ArraySize; i++){
            rapidjson::Value& dVec = inputIter[inputKey.size()][i].GetArray();
            data[i].resize(dVec.Size());
            for (size_t j = 0; j < dVec.Size(); j++){
                data[i][j] = dVec[j].GetString();
            }
        }
        return true;
    } catch (...){
        TypeErrorMessage(inputKey, "string");
        return false;
    }
}

bool inputData::HasKey(std::vector<std::string> inputKey){
    bool exists;
    rapidjson::Document::AllocatorType& a = inputs.GetAllocator();
    std::vector<rapidjson::Value> inputIter(inputKey.size()+1);
    inputIter[0].CopyFrom(inputs, a);

    for (size_t i = 0; i < inputKey.size(); i++){
        exists = inputIter[i].HasMember(inputKey[i].c_str());
        if (exists){
            inputIter[i+1] = inputIter[i][inputKey[i].c_str()];
        } else {
            return false;
        }
    }
    return true;
}