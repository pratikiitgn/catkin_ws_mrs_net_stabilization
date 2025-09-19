//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// Subscriber.cpp
//
// Code generation for function 'Subscriber'
//

// Include files
#include "Subscriber.h"
#include "geometry_msgs_PointStruct.h"
#include "pratik_2_myNode.h"
#include "pratik_2_myNode_types.h"
#include "mlroscpp_sub.h"

// Function Definitions
namespace coder {
namespace ros {
void Subscriber::callback()
{
  MessageCount = get_MessageCount() + 1.0;
  if (IsInitialized) {
    pratik_2_callback(MsgStruct.X, MsgStruct.Y, MsgStruct.Z);
  }
}

double Subscriber::get_MessageCount() const
{
  return MessageCount;
}

Subscriber *Subscriber::init()
{
  static const char topic[6]{'/', 'p', 'o', 'i', 'n', 't'};
  Subscriber *obj;
  obj = this;
  obj->IsInitialized = false;
  for (int i{0}; i < 6; i++) {
    obj->TopicName[i] = topic[i];
  }
  obj->BufferSize = 1.0;
  obj->MessageCount = 0.0;
  obj->MsgStruct = geometry_msgs_PointStruct();
  auto structPtr = (&obj->MsgStruct);
  obj->SubscriberHelper = std::unique_ptr<
      MATLABSubscriber<geometry_msgs::Point, geometry_msgs_PointStruct_T>>(
      new MATLABSubscriber<geometry_msgs::Point, geometry_msgs_PointStruct_T>(
          structPtr, [this] { this->callback(); })); //();
  MATLABSUBSCRIBER_createSubscriber(obj->SubscriberHelper, &obj->TopicName[0],
                                    6.0, obj->BufferSize);
  obj->callback();
  obj->IsInitialized = true;
  return obj;
}

} // namespace ros
} // namespace coder

// End of code generation (Subscriber.cpp)
