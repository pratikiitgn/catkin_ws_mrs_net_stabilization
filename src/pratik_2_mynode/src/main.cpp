#include "ros/ros.h"
#include <thread>
#include "pratik_2_myNode.h"

bool threadTerminating = false;

void threadFunction(void)
{
   try
   {
       pratik_2_myNode();
   }
   catch (std::runtime_error & e)
   {
       std::cout << "Caught exception: " << e.what() << std::endl;
   }
   catch (...)
   {
       std::cout << "Caught unknown exception, terminating the program." << std::endl;
   }
    threadTerminating = true;
    ros::shutdown();
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "pratik_2_myNode");
    ros::NodeHandlePtr MLROSNodePtr = ros::NodeHandlePtr(new ros::NodeHandle);
    std::thread threadObj(threadFunction);

    ros::spin();
    if (threadTerminating) {
    threadObj.join();
    }

    return 0;
}
