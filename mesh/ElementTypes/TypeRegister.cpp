#include "TypeRegister.h"
#include <iostream>


enum ElementTypes{ELEM_L9, ELEM_L3, ELEM_L4, ELEM_L2, 
                  ELEM_B2_1, ELEM_B1_1, ELEM_B2_2, ELEM_B2_3, ELEM_B2_4, ELEM_B1_2, ELEM_B1_3, ELEM_B1_4, 
                  ELEM_L3B, ELEM_T6B, ELEM_L2B, ELEM_T3B,
                  ELEM_T2_6, ELEM_T3_10, ELEM_T2_6B, ELEM_T3_10B,
                  ELEM_NURBS3_CUBE2, ELEM_NURBS3_PLANE2};

std::unordered_map<std::string, ElementTypes> ElementTypesNames = {
    {"L9", ELEM_L9},
    {"L3", ELEM_L3},
    {"L4", ELEM_L4},
    {"L2", ELEM_L2},
    {"B2_1", ELEM_B2_1},
    {"B1_1", ELEM_B1_1},
    {"B2_2", ELEM_B2_2},
    {"B2_3", ELEM_B2_3},
    {"B2_4", ELEM_B2_4},
    {"B1_2", ELEM_B1_2},
    {"B1_3", ELEM_B1_3},
    {"B1_4", ELEM_B1_4},
    {"L3B", ELEM_L3B},
    {"T6B", ELEM_T6B},
    {"L2B", ELEM_L2B},
    {"T3B", ELEM_T3B},
    {"T2_6", ELEM_T2_6},
    {"T3_10", ELEM_T3_10},
    {"T2_6B", ELEM_T2_6B},
    {"T3_10B", ELEM_T3_10B},
    {"NURBS3_Cube2",ELEM_NURBS3_CUBE2},
    {"NURBS3_Plane2",ELEM_NURBS3_PLANE2}
};


/// @brief Factory which initializes and returns references to the correct base element types
/// @param elemName Type of parametric element to initialize
/// @param ipcount1D Number of integration points per dimension
/// @return Pointer to parametric element type (remember to add to destructor)
BaseElemType* CreateElem(std::string elemName, int ipcount1D){

    ElementTypes EType;                                          
    BaseElemType* Elem;

    if (ElementTypesNames.count(elemName)){
        EType = ElementTypesNames[elemName];
    } else {
        std::stringstream err;
        err << "Element of type "+elemName+" not defined, valid options are: ";
        for(const auto& i : ElementTypesNames){
            err << i.first << ", ";
        }
        err << "\n";
        throw std::invalid_argument(err.str());
    }

    switch (EType) {
        case ELEM_L9:{
            Elem = new L9(ipcount1D);
        } break;
        case ELEM_L3:{
            Elem = new L3(ipcount1D);
        } break;
        case ELEM_L4:{
            Elem = new L4(ipcount1D);
        } break;
        case ELEM_B2_1:{
            Elem = new L4(ipcount1D);
        } break;
        case ELEM_L2:{
            Elem = new L2(ipcount1D);
        } break;
        case ELEM_B1_1:{
            Elem = new L2(ipcount1D);
        } break;
        case ELEM_B2_2:{
            Elem = new B2_2(ipcount1D);
        } break;
        case ELEM_B2_3:{
            Elem = new B2_3(ipcount1D);
        } break;
        case ELEM_B2_4:{
            Elem = new B2_4(ipcount1D);
        } break;
        case ELEM_B1_2:{
            Elem = new B1_2(ipcount1D);
        } break;
        case ELEM_B1_3:{
            Elem = new B1_3(ipcount1D);
        } break;
        case ELEM_B1_4:{
            Elem = new B1_4(ipcount1D);
        } break;
        case ELEM_L3B:{
            Elem = new L3B(ipcount1D);
        } break;
        case ELEM_T6B:{
            Elem = new T6B(ipcount1D);
        } break;
        case ELEM_L2B:{
            Elem = new L2B(ipcount1D);
        } break;
        case ELEM_T3B:{
            Elem = new T3B(ipcount1D);
        } break;
        case ELEM_T2_6:{
            Elem = new T2_6(ipcount1D);
        } break;
        case ELEM_T3_10:{
            Elem = new T3_10(ipcount1D);
        } break;
        case ELEM_T2_6B:{
            Elem = new T2_6B(ipcount1D);
        } break;
        case ELEM_T3_10B:{
            Elem = new T3_10B(ipcount1D);
        } break;
        case ELEM_NURBS3_CUBE2:{
            Elem = new NURBS3_Cube2(ipcount1D);
        } break;
        case ELEM_NURBS3_PLANE2:{
            Elem = new NURBS3_Plane2(ipcount1D);
        } break;
        default:{
            std::stringstream err;
            err << "Element of type "+elemName+" not defined, valid options are: ";
            for(const auto& i : ElementTypesNames){
                err << i.first << ", ";
            }
            err << "\n";
            throw std::invalid_argument(err.str());
        }
    }

    return Elem;
}