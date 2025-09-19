//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// geometry_msgs_PointStruct.cpp
//
// Code generation for function 'geometry_msgs_PointStruct'
//

// Include files
#include "geometry_msgs_PointStruct.h"
#include "pratik_2_myNode_types.h"

// Function Definitions
geometry_msgs_PointStruct_T geometry_msgs_PointStruct()
{
  static const geometry_msgs_PointStruct_T b_msg{
      {'g', 'e', 'o', 'm', 'e', 't', 'r', 'y', '_', 'm', 's', 'g', 's', '/',
       'P', 'o', 'i', 'n', 't'}, // MessageType
      0.0,                       // X
      0.0,                       // Y
      0.0                        // Z
  };
  geometry_msgs_PointStruct_T msg;
  msg = b_msg;
  //  Message struct definition for geometry_msgs/Point
  //(&b_msg);
  return msg;
}

// End of code generation (geometry_msgs_PointStruct.cpp)
