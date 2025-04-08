#ifndef TYPEREGISTER_H
#define TYPEREGISTER_H

#include "BaseElemType.h"
#include "Lagrange/L9.h"
#include "Lagrange/L3.h"
#include "Lagrange/L4.h"
#include "Lagrange/L2.h"
#include "BernStein/T6B.h"
#include "BernStein/L3B.h"
#include "BernStein/T3B.h"
#include "BernStein/L2B.h"
#include "Bezier/B2_2.h"
#include "Bezier/B2_3.h"
#include "Bezier/B2_4.h"
#include "Bezier/B1_2.h"
#include "Bezier/B1_3.h"
#include "Bezier/B1_4.h"
#include "3D/T2_6.h"
#include "3D/T3_10.h"
#include "3D/T2_6B.h"
#include "3D/T3_10B.h"
#include "NURBS/NURBS3_Cube2.h"
#include "NURBS/NURBS3_Plane2.h"
#include <string>

BaseElemType* CreateElem(std::string elemName, int ipcount1D);

#endif