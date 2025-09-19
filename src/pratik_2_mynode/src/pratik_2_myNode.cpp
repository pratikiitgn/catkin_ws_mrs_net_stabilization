//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// pratik_2_myNode.cpp
//
// Code generation for function 'pratik_2_myNode'
//

// Include files
#include "pratik_2_myNode.h"
#include "Subscriber.h"
#include "pause.h"
#include "pratik_2_myNode_data.h"
#include "pratik_2_myNode_initialize.h"
#include "validate_print_arguments.h"
#include <cstdio>

// Function Definitions
void pratik_2_callback(double msg_X, double msg_Y, double msg_Z)
{
  double validatedHoleFilling[3];
  //  Subscriber callback function
  coder::internal::validate_print_arguments(msg_X, msg_Y, msg_Z,
                                            validatedHoleFilling);
  std::printf("(X,Y,Z): (%f,%f,%f)\n", validatedHoleFilling[0],
              validatedHoleFilling[1], validatedHoleFilling[2]);
  std::fflush(stdout);
}

void pratik_2_myNode()
{
  coder::ros::Subscriber sub;
  char varargin_1[7];
  if (!isInitialized_pratik_2_myNode) {
    pratik_2_myNode_initialize();
  }
  sub.init();
  for (int i{0}; i < 6; i++) {
    varargin_1[i] = sub.TopicName[i];
  }
  varargin_1[6] = '\x00';
  std::printf("Created %s subscriber", &varargin_1[0]);
  std::fflush(stdout);
  while (1) {
    std::printf("Node is alive .. \n");
    std::fflush(stdout);
    coder::pause();
  }
}

// End of code generation (pratik_2_myNode.cpp)
