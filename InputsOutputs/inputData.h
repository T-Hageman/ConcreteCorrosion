#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include <petsc.h>
#include <fstream>
#include <iostream>
#include <vector>

/// @brief Class to handle input data
class inputData {
    public:
        inputData();
        inputData(char const* filename);
        ~inputData();

        void print();

        std::string GetType(std::vector<std::string> inputKey);

        void GetRequired(std::string &data, std::vector<std::string> inputKey);
        bool GetOptional(std::string &data, std::vector<std::string> inputKey);

        void GetRequired(double      &data, std::vector<std::string> inputKey);
        bool GetOptional(double      &data, std::vector<std::string> inputKey);

        void GetRequired(bool        &data, std::vector<std::string> inputKey);
        bool GetOptional(bool        &data, std::vector<std::string> inputKey);

        void GetRequired(size_t        &data, std::vector<std::string> inputKey);
        bool GetOptional(size_t        &data, std::vector<std::string> inputKey);

        void GetRequired(uint        &data, std::vector<std::string> inputKey);
        bool GetOptional(uint        &data, std::vector<std::string> inputKey);

        void GetRequired(int        &data, std::vector<std::string> inputKey);
        bool GetOptional(int        &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<std::string> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<std::string> &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<int> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<int> &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<bool> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<bool> &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<double> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<double> &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<size_t> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<size_t> &data, std::vector<std::string> inputKey);

        void GetRequired(std::vector<std::vector<std::string>> &data, std::vector<std::string> inputKey);
        bool GetOptional(std::vector<std::vector<std::string>> &data, std::vector<std::string> inputKey);

        bool HasKey(std::vector<std::string> inputKey);
    private:
        void parse_file( char const* filename );
        void KeyErrorMessage(std::vector<std::string> inputKey);
        void TypeErrorMessage(std::vector<std::string> inputKey, std::string typestring);

        rapidjson::Document inputs; // Object holding actual input data, accessible through ["string"].getDouble etc.    
};

#endif