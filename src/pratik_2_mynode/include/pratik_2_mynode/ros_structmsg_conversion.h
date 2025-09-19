#ifndef _ROS_STRUCTMSG_CONVERSION_H_
#define _ROS_STRUCTMSG_CONVERSION_H_

#include "ros_structmsg_declare.h"
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include "pratik_2_myNode_types.h"
#include "mlroscpp_msgconvert_utils.h"


void struct2msg(geometry_msgs::Point* msgPtr, geometry_msgs_PointStruct_T const* structPtr);
void msg2struct(geometry_msgs_PointStruct_T* structPtr, geometry_msgs::Point const* msgPtr);


#endif
