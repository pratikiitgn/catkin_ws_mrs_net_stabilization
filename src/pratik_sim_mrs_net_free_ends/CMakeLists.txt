cmake_minimum_required(VERSION 3.5)
project(pratik_sim_mrs_net_free_ends)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-Wall)
add_definitions(-Wextra)

set(CATKIN_DEPENDENCIES
  cmake_modules
  dynamic_reconfigure
  geometry_msgs
  mrs_lib
  mrs_msgs
  mrs_uav_hw_api
  nav_msgs
  nodelet
  roscpp
  rosgraph_msgs
  sensor_msgs
  std_msgs
  std_srvs
  )

find_package(catkin REQUIRED COMPONENTS
  visualization_msgs
  ${CATKIN_DEPENDENCIES}
  )

generate_dynamic_reconfigure_options(
  cfg/sim_mrs_net_free_ends.cfg
  )

set(LIBRARIES
  PratikSimMrsNetFreeEnds_Simulator
  PratikSimMrsNetFreeEnds_SimulatorSystemRos
  PratikSimMrsNetFreeEnds_HwApiPlugin
  )

catkin_package(
  LIBRARIES ${LIBRARIES}
  INCLUDE_DIRS include 
  include/pratik_sim_mrs_net_free_ends 
  include/pratik_sim_mrs_net_free_ends/uav_system 
  include/pratik_sim_mrs_net_free_ends/uav_system/controllers
  CATKIN_DEPENDS ${CATKIN_DEPENDENCIES} visualization_msgs
  )

include_directories(
  include
  include/pratik_sim_mrs_net_free_ends
  include/pratik_sim_mrs_net_free_ends/uav_system
  include/pratik_sim_mrs_net_free_ends/uav_system/controllers
  ${catkin_INCLUDE_DIRS}
  ${dynamic_reconfigure_PACKAGE_PATH}/cmake/cfgbuild.cmake
  )

# ROS Wrapper for the Multirotor UAV System + controllers

add_library(PratikSimMrsNetFreeEnds_SimulatorSystemRos
  src/uav_system_ros.cpp
  )

add_dependencies(PratikSimMrsNetFreeEnds_SimulatorSystemRos
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimMrsNetFreeEnds_SimulatorSystemRos
  ${catkin_LIBRARIES}
  ${Eigen_LIBRARIES}
  )

# Multirotor simulator
# - ROS node with the MRS Multirotor Simulator

add_library(PratikSimMrsNetFreeEnds_Simulator
  src/sim_mrs_net_free_ends.cpp
  )

add_dependencies(PratikSimMrsNetFreeEnds_Simulator
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimMrsNetFreeEnds_Simulator
  PratikSimMrsNetFreeEnds_SimulatorSystemRos
  ${catkin_LIBRARIES}
  ${Eigen_LIBRARIES}
  )

# Hw API Plugin
# - pluginlib library for interfacing the Simulator with the MRS UAV System

add_library(PratikSimMrsNetFreeEnds_HwApiPlugin
  src/hw_api_plugin.cpp
  )

add_dependencies(PratikSimMrsNetFreeEnds_HwApiPlugin
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimMrsNetFreeEnds_HwApiPlugin
  ${catkin_LIBRARIES}
  ${Eigen_LIBRARIES}
  )

## --------------------------------------------------------------
## |                           Install                          |
## --------------------------------------------------------------

install(TARGETS ${LIBRARIES}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
  )

install(DIRECTORY launch config
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
  )

install(DIRECTORY tmux/
  USE_SOURCE_PERMISSIONS
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/tmux
  )

install(DIRECTORY include/pratik_sim_mrs_net_free_ends/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  )

install(DIRECTORY scripts/
  USE_SOURCE_PERMISSIONS
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )

install(DIRECTORY scripts
  USE_SOURCE_PERMISSIONS
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
  )

install(FILES nodelets.xml plugins.xml
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
  )
