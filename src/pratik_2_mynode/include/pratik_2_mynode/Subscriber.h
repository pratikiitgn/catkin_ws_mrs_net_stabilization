//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// Subscriber.h
//
// Code generation for function 'Subscriber'
//

#ifndef SUBSCRIBER_H
#define SUBSCRIBER_H

// Include files
#include "pratik_2_myNode_types.h"
#include "rtwtypes.h"
#include "mlroscpp_sub.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace ros {
class Subscriber {
public:
  Subscriber *init();
  void callback();
  double get_MessageCount() const;
  char TopicName[6];
  double BufferSize;
  double MessageCount;

private:
  std::unique_ptr<
      MATLABSubscriber<geometry_msgs::Point, geometry_msgs_PointStruct_T>>
      SubscriberHelper;
  geometry_msgs_PointStruct_T MsgStruct;
  bool IsInitialized;
};

} // namespace ros
} // namespace coder

#endif
// End of code generation (Subscriber.h)
