cmake_minimum_required(VERSION 3.5)
project(pratik_sim_n_link_constraint)

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
  ${CATKIN_DEPENDENCIES}
  )

generate_dynamic_reconfigure_options(
  cfg/multirotor_simulator.cfg
  )

set(LIBRARIES
  PratikSimNLinkConstraint_Simulator
  PratikSimNLinkConstraint_SimulatorSystemRos
  PratikSimNLinkConstraint_HwApiPlugin
  )

catkin_package(
  LIBRARIES ${LIBRARIES}
  INCLUDE_DIRS include include/mrs_multirotor_simulator include/mrs_multirotor_simulator/uav_system include/mrs_multirotor_simulator/uav_system/controllers
  CATKIN_DEPENDS ${CATKIN_DEPENDENCIES}
  )

include_directories(
  include
  include/mrs_multirotor_simulator
  include/mrs_multirotor_simulator/uav_system
  include/mrs_multirotor_simulator/uav_system/controllers
  ${catkin_INCLUDE_DIRS}
  ${dynamic_reconfigure_PACKAGE_PATH}/cmake/cfgbuild.cmake
  )

# ROS Wrapper for the Multirotor UAV System + controllers

add_library(PratikSimNLinkConstraint_SimulatorSystemRos
  src/uav_system_ros.cpp
  )

add_dependencies(PratikSimNLinkConstraint_SimulatorSystemRos
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimNLinkConstraint_SimulatorSystemRos
  ${catkin_LIBRARIES}
  ${Eigen_LIBRARIES}
  )

# Multirotor simulator
# - ROS node with the MRS Multirotor Simulator

add_library(PratikSimNLinkConstraint_Simulator
  src/multirotor_simulator.cpp
  )

add_dependencies(PratikSimNLinkConstraint_Simulator
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimNLinkConstraint_Simulator
  PratikSimNLinkConstraint_SimulatorSystemRos
  ${catkin_LIBRARIES}
  ${Eigen_LIBRARIES}
  )

# Hw API Plugin
# - pluginlib library for interfacing the Simulator with the MRS UAV System

add_library(PratikSimNLinkConstraint_HwApiPlugin
  src/hw_api_plugin.cpp
  )

add_dependencies(PratikSimNLinkConstraint_HwApiPlugin
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(PratikSimNLinkConstraint_HwApiPlugin
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

install(DIRECTORY include/mrs_multirotor_simulator/
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
