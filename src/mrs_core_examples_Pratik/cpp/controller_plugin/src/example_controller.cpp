/* includes //{ */

#include <ros/ros.h>
#include <ros/package.h>
#include <nodelet/nodelet.h>

// #include <pid.hpp>

/* for loading dynamic parameters while the nodelet is running */
#include <mrs_uav_managers/controller.h>
#include <dynamic_reconfigure/server.h>

#include <example_controller_plugin/example_controllerConfig.h>

#include <mrs_lib/param_loader.h>
#include <mrs_lib/subscribe_handler.h>
#include <mrs_lib/mutex.h>
#include <mrs_lib/attitude_converter.h>
#include <mrs_lib/msg_extractor.h>
#include <mrs_lib/geometry/misc.h>

#include <mrs_msgs/ControlManagerDiagnostics.h>

// | ----------------- Calling required libraries ----------------- |
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// | ----------------- Calling required libraries from gazebo ----------------- |
// #include <gazebo_msgs/LinkStates.h>

/* for storing information about the state of the uav (position) */
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PointStamped.h>

// #include "mrs_multirotor_simulator/include/mrs_multirotor_simulator/uav_system/multirotor_model_cable_suspended_load.hpp"
#include <nav_msgs/Odometry.h>

// Cable-suspended load

// #include <rate.h>

// | ----------------- Basic math variables ----------------- |
Eigen::Vector3d e1(1.0,0.0,0.0);
Eigen::Vector3d e2(0.0,1.0,0.0);
Eigen::Vector3d e3(0.0,0.0,1.0);

// | ----------------- System parameters     ---------------- |

float mq                = 2.0;      // in kgs mass of the quadcopter
float mp                = 0.1;      // in kgs mass of the payload
float g_acceleration    = 9.81;     // in m/s^2
float PI_value          = 3.1415926535;

// | ----------------- Time related variables ----------------- |
float MRS_text_start_time = 0.0;
double initial_ros_time_custom_controller = 0.0;

float t1_MRS_traj = 0.0;
float t2_MRS_traj = 0.0;

float t1_straight_traj = 0.0;
float t2_straight_traj = 0.0;

// | ----------------- Position related variables for MRS text ----------------- |
float sty_MRS_traj = 0.0;
float stz_MRS_traj = 0.0;
float eny_MRS_traj = 0.0;
float enz_MRS_traj = 0.0;

// | ----------------- Position related variables for Striaght Line ----------------- |
float sty_straight_traj = 0.0;
float stz_straight_traj = 0.0;
float eny_straight_traj = 0.0;
float enz_straight_traj = 0.0;

// | ----------------- Quadcopter State ----------------- |

Eigen::Vector3d   pos_of_quad(0.0,0.0,0.0);
Eigen::Vector3d   vel_of_quad(0.0,0.0,0.0);

// | ----------------- Desired quadcopter State ----------------- |

float des_quad_x      = 0.0;
float des_quad_y      = 0.0;
float des_quad_z      = 0.0;

float des_quad_x_dot  = 0.0;
float des_quad_y_dot  = 0.0;
float des_quad_z_dot  = 0.0;

float des_quad_x_dot_dot  = 0.0;
float des_quad_y_dot_dot  = 0.0;
float des_quad_z_dot_dot  = 0.0;

Eigen::Vector3d des_pos_of_quad(0.0,0.0,0.0);
Eigen::Vector3d des_vel_of_quad(0.0,0.0,0.0);
Eigen::Vector3d des_acc_of_quad(0.0,0.0,0.0);

Eigen::Vector3d b_1_des(1.0,0.0,0.0);
Eigen::Vector3d b_2_des(0.0,1.0,0.0);
Eigen::Vector3d b_3_des(0.0,0.0,1.0);

float commanded_yaw_angle = 0.0;

Eigen::Vector3d b_1_c(1.0,0.0,0.0);

Eigen::Matrix3d R_des;
Eigen::Vector3d des_rpy;
// | ----------------- Cable attitude State ----------------- |

Eigen::Vector3d     q (0.0,0.0,-1.0);
Eigen::Vector3d q_dot (0.0,0.0, 0.0);
Eigen::Vector3d q_old (0.0,0.0, 0.0);

// | ----------------- Desired cable attitude State ----------------- |

Eigen::Vector3d     q_d (0.0,0.0,-1.0);
Eigen::Vector3d q_d_dot (0.0,0.0, 0.0);

// | ----------------- Payload position State ----------------- |

Eigen::Vector3d pos_of_payload(0.0,0.0,0.0);

// | ----------------- Error definition ----------------- |
Eigen::Vector3d     e_x_q     (0.0,0.0,0.0);
Eigen::Vector3d     e_x_q_dot (0.0,0.0,0.0);

Eigen::Vector3d     e_q       (0.0,0.0,0.0);
Eigen::Vector3d     e_q_dot   (0.0,0.0,0.0);

// | ----------------- Custom Gains ----------------- |

float kx_1        = 0.0;
float kx_2        = 0.0;
float kx_3        = 0.0;

float kx_1_dot    = 0.0;
float kx_2_dot    = 0.0;
float kx_3_dot    = 0.0;

float kq_1        = 0.0;
float kq_2        = 0.0;
float kq_3        = 0.0;

float kq_1_dot    = 0.0;
float kq_2_dot    = 0.0;
float kq_3_dot    = 0.0;

Eigen::Array3d kx(0.0,0.0,0.0);
Eigen::Array3d kx_dot(0.0,0.0,0.0);

Eigen::Array3d kq(0.0,0.0,0.0);
Eigen::Array3d kq_dot(0.0,0.0,0.0);

// | ----------------- High level commands ----------------- |
double desired_thrust_force = 0.0;

// | ----------------- Thrust force ----------------- |
Eigen::Vector3d u_control_input (0.0,0.0,0.0);
Eigen::Vector3d u_cable_input   (0.0,0.0,0.0);
Eigen::Vector3d u_quad_input    (0.0,0.0,0.0);

//}

namespace example_controller_plugin
{

namespace example_controller
{

/* //{ class ExampleController */

class ExampleController : public mrs_uav_managers::Controller {

public:
  bool initialize(const ros::NodeHandle& nh, std::shared_ptr<mrs_uav_managers::control_manager::CommonHandlers_t> common_handlers,
                  std::shared_ptr<mrs_uav_managers::control_manager::PrivateHandlers_t> private_handlers);

  bool activate(const ControlOutput& last_control_output);

  void deactivate(void);

  void updateInactive(const mrs_msgs::UavState& uav_state, const std::optional<mrs_msgs::TrackerCommand>& tracker_command);

  ////////////////////////////////////////////////
  //// for custom controller
  ////////////////////////////////////////////////
  float distance_bt_two_pts(Eigen::Vector3d A, Eigen::Vector3d B);
  float min_acc_first_coefficient(float t1, float t2, float st, float en);
  float min_acc_second_coefficient(float t1, float t2, float st, float en);
  float min_acc_third_coefficient(float t1, float t2, float st, float en);
  float min_acc_fourth_coefficient(float t1, float t2, float st, float en);
  float clipping_angle(float max_value, float current_angle);
  Eigen::Vector3d Matrix_vector_mul(Eigen::Matrix3d R, Eigen::Vector3d v);
  float clipping_net_thrust_force(float max_value, float current_thrust);
  Eigen::Vector3d clipping_e_x_q(Eigen::Vector3d e_x_q_vector);
  Eigen::Vector3d clipping_e_x_q_dot(Eigen::Vector3d e_x_q_dot_vector);
  Eigen::Vector3d Rotation_matrix_to_Euler_angle(Eigen::Matrix3d R);

  ////////////////////////////////////////////////
  //// for custom controller
  ////////////////////////////////////////////////

  ControlOutput updateActive(const mrs_msgs::UavState& uav_state, const mrs_msgs::TrackerCommand& tracker_command);

  const mrs_msgs::ControllerStatus getStatus();

  void switchOdometrySource(const mrs_msgs::UavState& new_uav_state);

  void resetDisturbanceEstimators(void);

  const mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr setConstraints(const mrs_msgs::DynamicsConstraintsSrvRequest::ConstPtr& cmd);

  // 

private:
  ros::NodeHandle nh_;

  bool is_initialized_ = false;
  bool is_active_      = false;

  std::shared_ptr<mrs_uav_managers::control_manager::CommonHandlers_t>  common_handlers_;
  std::shared_ptr<mrs_uav_managers::control_manager::PrivateHandlers_t> private_handlers_;

  // | ------------------------ uav state ----------------------- |

  mrs_msgs::UavState uav_state_;
  std::mutex         mutex_uav_state_;

  // | --------------- dynamic reconfigure server --------------- |

  boost::recursive_mutex                                      mutex_drs_;
  typedef example_controller_plugin::example_controllerConfig DrsConfig_t;
  typedef dynamic_reconfigure::Server<DrsConfig_t>            Drs_t;
  boost::shared_ptr<Drs_t>                                    drs_;
  void                                                        callbackDrs(example_controller_plugin::example_controllerConfig& config, uint32_t level);
  DrsConfig_t                                                 drs_params_;
  std::mutex                                                  mutex_drs_params_;

  // | ----------------------- constraints ---------------------- |

  mrs_msgs::DynamicsConstraints constraints_;
  std::mutex                    mutex_constraints_;

  // | --------- throttle generation and mass estimation -------- |

  double _uav_mass_;
  double uav_mass_difference_;

  // | ------------------ activation and output ----------------- |

  ControlOutput last_control_output_;
  ControlOutput activation_control_output_;

  ros::Time         last_update_time_;
  std::atomic<bool> first_iteration_ = true;

  // | ------------------------ integrals ----------------------- |

  Eigen::Vector2d Ib_b_;  // body error integral in the body frame
  Eigen::Vector2d Iw_w_;  // world error integral in the world_frame

  // | ------------------------- rampup ------------------------- |

  bool   _rampup_enabled_ = false;
  double _rampup_speed_;

  bool      rampup_active_ = false;
  double    rampup_throttle_;
  int       rampup_direction_;
  double    rampup_duration_;
  ros::Time rampup_start_time_;
  ros::Time rampup_last_time_;

  // | --------------------- timer callbacks -------------------- |

  // | ---------------------- msg callbacks --------------------- |
  mrs_lib::SubscribeHandler<nav_msgs::Odometry>                  sh_cable_states;
  void              callback_cable_states(const nav_msgs::Odometry::ConstPtr msg);

};

//}

// --------------------------------------------------------------
// |                   controller's interface                   |
// --------------------------------------------------------------

/* //{ initialize() */

bool ExampleController::initialize(const ros::NodeHandle& nh, std::shared_ptr<mrs_uav_managers::control_manager::CommonHandlers_t> common_handlers,
                                   std::shared_ptr<mrs_uav_managers::control_manager::PrivateHandlers_t> private_handlers) {

  nh_ = nh;

  common_handlers_  = common_handlers;
  private_handlers_ = private_handlers;

  _uav_mass_ = common_handlers->getMass();

  last_update_time_ = ros::Time(0);
  initial_ros_time_custom_controller =ros::Time::now().toSec();

  ros::Time::waitForValid();

  // | ------------------- loading parameters ------------------- |

  bool success = true;

  // FYI
  // This method will load the file using `rosparam get`
  //   Pros: you can the full power of the official param loading
  //   Cons: it is slower
  //
  // Alternatives:
  //   You can load the file directly into the ParamLoader as shown below.

  success *= private_handlers->loadConfigFile(ros::package::getPath("example_controller_plugin") + "/config/example_controller.yaml");

  if (!success) {
    return false;
  }

  mrs_lib::ParamLoader param_loader(nh_, "ExampleController");

  // This is the alternaive way of loading the config file.
  //
  // Files loaded using this method are prioritized over ROS params.
  //
  // param_loader.addYamlFile(ros::package::getPath("example_tracker_plugin") + "/config/example_tracker.yaml");

  param_loader.loadParam("desired_roll",          drs_params_.roll);
  param_loader.loadParam("desired_pitch",         drs_params_.pitch);
  param_loader.loadParam("desired_yaw",           drs_params_.yaw);
  param_loader.loadParam("desired_thrust_force",  drs_params_.force);
  param_loader.loadParam("kx_1_value",      kx_1);
  param_loader.loadParam("kx_2_value",      kx_2);
  param_loader.loadParam("kx_3_value",      kx_3);
  param_loader.loadParam("kx_1_dot_value",  kx_1_dot); 
  param_loader.loadParam("kx_2_dot_value",  kx_2_dot); 
  param_loader.loadParam("kx_3_dot_value",  kx_3_dot);

  param_loader.loadParam("kq_1_value",      kq_1);
  param_loader.loadParam("kq_2_value",      kq_2);
  param_loader.loadParam("kq_3_value",      kq_3);
  param_loader.loadParam("kq_1_dot_value",  kq_1_dot); 
  param_loader.loadParam("kq_2_dot_value",  kq_2_dot); 
  param_loader.loadParam("kq_3_dot_value",  kq_3_dot);

  // | -------- initialize a publisher -------- |

  // | ------------------ finish loading params ----------------- |

  if (!param_loader.loadedSuccessfully()) {
    ROS_ERROR("[ExampleController]: could not load all parameters!");
    return false;
  }

  // | --------------- dynamic reconfigure server --------------- |

  drs_.reset(new Drs_t(mutex_drs_, nh_));
  drs_->updateConfig(drs_params_);
  Drs_t::CallbackType f = boost::bind(&ExampleController::callbackDrs, this, _1, _2);
  drs_->setCallback(f);

  // | ------------------ initialize subscribers ----------------- |

  mrs_lib::SubscribeHandlerOptions shopts;
  shopts.nh                 = nh;
  shopts.node_name          = "ExampleController";
  shopts.no_message_timeout = ros::Duration(1.0);
  shopts.threadsafe         = true;
  shopts.autostart          = true;
  shopts.queue_size         = 10;
  shopts.transport_hints    = ros::TransportHints().tcpNoDelay();

  sh_cable_states           = mrs_lib::SubscribeHandler<nav_msgs::Odometry>(shopts, "/multirotor_simulator/uav1/cable_state",
                                                                                            &ExampleController::callback_cable_states, this);

  // initialize the integrals
  uav_mass_difference_ = 0;

  // | ----------------------- finish init ---------------------- |

  ROS_INFO("[ExampleController]: initialized");

  is_initialized_ = true;

  return true;
}

//}

/* //{ activate() */

bool ExampleController::activate(const ControlOutput& last_control_output) {

  activation_control_output_ = last_control_output;

  double activation_mass = _uav_mass_;

  if (activation_control_output_.diagnostics.mass_estimator) {
    uav_mass_difference_ = activation_control_output_.diagnostics.mass_difference;
    activation_mass += uav_mass_difference_;
    ROS_INFO("[ExampleController]: setting mass difference from the last control output: %.2f kg", uav_mass_difference_);
  }

  last_control_output_.diagnostics.controller_enforcing_constraints = false;

  first_iteration_ = true;

  is_active_ = true;

  ROS_INFO("[ExampleController]: activated");

  return true;
}

//}

/* //{ deactivate() */

void ExampleController::deactivate(void) {

  is_active_            = false;
  first_iteration_      = false;
  uav_mass_difference_  = 0;

  ROS_INFO("[ExampleController]: deactivated");
}

//}

/* updateInactive() //{ */

void ExampleController::updateInactive(const mrs_msgs::UavState& uav_state, [[maybe_unused]] const std::optional<mrs_msgs::TrackerCommand>& tracker_command) {

  mrs_lib::set_mutexed(mutex_uav_state_, uav_state, uav_state_);

  last_update_time_ = uav_state.header.stamp;

  first_iteration_ = false;
}

//}

/* //{ updateActive() */

ExampleController::ControlOutput ExampleController::updateActive(const mrs_msgs::UavState& uav_state, const mrs_msgs::TrackerCommand& tracker_command) {

  auto drs_params = mrs_lib::get_mutexed(mutex_drs_params_, drs_params_);

  mrs_lib::set_mutexed(mutex_uav_state_, uav_state, uav_state_);

  // clear all the optional parts of the result
  last_control_output_.desired_heading_rate          = {};
  last_control_output_.desired_orientation           = {};
  last_control_output_.desired_unbiased_acceleration = {};
  last_control_output_.control_output                = {};

  if (!is_active_) {
    return last_control_output_;
  }

  // | ---------- calculate dt from the last iteration ---------- |
  double dt;

  if (first_iteration_) {
    // dt               = 0.01;
    dt               = 0.004;
    first_iteration_ = false;
  } else {
    dt = (uav_state.header.stamp - last_update_time_).toSec();
  }

  last_update_time_ = uav_state.header.stamp;

  if (fabs(dt) < 0.001) {

    ROS_DEBUG("[ExampleController]: the last odometry message came too close (%.2f s)!", dt);
    // dt = 0.01;
    dt               = 0.004;
  }

  // | -------- check for the available output modalities ------- |

  // you can decide what to return, but it needs to be available
  if (common_handlers_->control_output_modalities.attitude) {
    ROS_INFO_THROTTLE(1.0, "[ExampleController]: desired attitude output modality is available");
  }

  // | ---------- extract the detailed model parameters --------- |

  if (common_handlers_->detailed_model_params) {

    mrs_uav_managers::control_manager::DetailedModelParams_t detailed_model_params = common_handlers_->detailed_model_params.value();

    // ROS_INFO_STREAM_THROTTLE(1.0, "[ExampleController]: UAV inertia is: " << detailed_model_params.inertia);
  }

////////////////////////////////////////////////
//////         Custom controller starts
////////////////////////////////////////////////

  // | -------------- prepare the control reference ------------- |

  geometry_msgs::PoseStamped position_reference;

  position_reference.header           = tracker_command.header;
  position_reference.pose.position    = tracker_command.position;
  position_reference.pose.orientation = mrs_lib::AttitudeConverter(0, 0, 0).setHeading(tracker_command.heading);

  des_pos_of_quad[0] = tracker_command.position.x;
  des_pos_of_quad[1] = tracker_command.position.y;
  des_pos_of_quad[2] = tracker_command.position.z;

  des_vel_of_quad[0] = tracker_command.velocity.x;
  des_vel_of_quad[1] = tracker_command.velocity.y;
  des_vel_of_quad[2] = tracker_command.velocity.z;

  des_acc_of_quad[0] = tracker_command.acceleration.x;
  des_acc_of_quad[1] = tracker_command.acceleration.y;
  des_acc_of_quad[2] = tracker_command.acceleration.z;

  commanded_yaw_angle  = tracker_command.heading;

  // | ---------------- Custom PD Controller for altitude control --------------- |

  MRS_text_start_time = ros::Time::now().toSec() - initial_ros_time_custom_controller;
  ROS_INFO_STREAM_THROTTLE(1, "[ExampleController]: Current Time: " << MRS_text_start_time);

  // Getting positional state of the drone
  pos_of_quad[0] = uav_state.pose.position.x;
  pos_of_quad[1] = uav_state.pose.position.y;
  pos_of_quad[2] = uav_state.pose.position.z;
  vel_of_quad[0] = uav_state.velocity.linear.x;
  vel_of_quad[1] = uav_state.velocity.linear.y;
  vel_of_quad[2] = uav_state.velocity.linear.z;

  Eigen::Quaterniond quad_rot_in_quat(uav_state.pose.orientation.w, uav_state.pose.orientation.x, uav_state.pose.orientation.y, uav_state.pose.orientation.z);

  Eigen::Matrix3d R_curr = mrs_lib::AttitudeConverter(quad_rot_in_quat);

  // ROS_INFO_STREAM_THROTTLE(0.5, "Just Debugging" << R_curr);

  // // quad_rot_in_quat.normalize();
  // Eigen::Matrix3d quad_R = quad_rot_in_quat.toRotationMatrix();

  // ROS_INFO_STREAM_THROTTLE(0.5, "Just Debugging" << pos_of_quad);
  // ROS_INFO_STREAM_THROTTLE(1, "[ExampleController]: cable Attitude" << uav_state.cable);
  // MultirotorModel::State &state;

  // | ---------------- Get the gains values --------------- |

  kx      << kx_1 ,     kx_2 ,    kx_3;
  kx_dot  << kx_1_dot , kx_2_dot, kx_3_dot;

  kq      << kq_1 ,     kq_2 ,    kq_3;
  kq_dot  << kq_1_dot , kq_2_dot, kq_3_dot;

  // | ---------------- Error computation --------------- |

  e_x_q       = des_pos_of_quad - pos_of_quad;
  e_x_q       = clipping_e_x_q(e_x_q);
  e_x_q_dot   = des_vel_of_quad - vel_of_quad;
  e_x_q_dot   = clipping_e_x_q_dot(e_x_q_dot);

  // ROS_INFO("x des: %2.2f, y des: %2.2f, z des: %2.2f", des_pos_of_quad(0), des_pos_of_quad(1), des_pos_of_quad(2));
  // ROS_INFO("x pos: %2.2f, y pos: %2.2f, z pos: %2.2f", pos_of_quad(0), pos_of_quad(1), pos_of_quad(2));
  // ROS_INFO("x err: %2.2f, y err: %2.2f, z err: %2.2f", e_x_q(0), e_x_q(1), e_x_q(2));
  // ROS_INFO("----------------");

  // ROS_INFO_STREAM_THROTTLE(0.2, "x pos:" << pos_of_quad(0) ", y pos:" << pos_of_quad(1) ", z pos:" << pos_of_quad(2));
  // ROS_INFO_STREAM_THROTTLE(0.2, "x err:" << e_x_q(0) ", y err:" << e_x_q(1) ", z err:" << e_x_q(2));

  // ROS_INFO_STREAM_THROTTLE(0.2, "error   :" << e_x_q);
  // ROS_INFO_STREAM_THROTTLE(0.2, "quad pos:" << pos_of_quad);

  // ROS_INFO_STREAM_THROTTLE(0.5, "Just Debugging" << e_q);

  // | ---------------- prepare the quadcopter control output --------------- |

  // Eigen::Vector3d feed_forward   = (mq + mp) * g_acceleration * e3;
  // Eigen::Vector3d feed_forward   = _uav_mass_ * g_acceleration * e3 + _uav_mass_ * des_acc_of_quad;

  double total_mass                 = _uav_mass_ + uav_mass_difference_;
  Eigen::Vector3d feed_forward      = total_mass * g_acceleration * e3;
  Eigen::Vector3d position_feedback = kx     * e_x_q.array();
  Eigen::Vector3d velocity_feedback = kx_dot * e_x_q_dot.array();

  u_quad_input      = position_feedback + velocity_feedback + feed_forward;

  // ROS_INFO("x err: %2.2f, y err: %2.2f, z err: %2.2f", e_x_q(0), e_x_q(1), e_x_q(2));
  ROS_INFO_THROTTLE(0.5, "Mass Estimator: %2.2f", total_mass);

  // | ---------------- Get Desired Cable Attitude
  // q_d               = -u_quad_input.normalized();

  // | ---------------- prepare the cable attitude control output --------------- |
  
  e_q               = q.cross(q.cross(q_d));
  e_q_dot           = q_dot - (q_d.cross(q_d_dot)).cross(q);
  u_cable_input     = kq * e_q.array()  + kq_dot * e_q_dot.array();

  // | ---------------- prepare the final control output --------------- |
  u_control_input   = u_quad_input;
  // u_control_input  = u_quad_input + u_cable_input;

  if (u_control_input(2) < 0) {
    ROS_WARN_THROTTLE(1.0, "[ExampleController]: the calculated downwards desired force is negative (%.2f) -> mitigating flip", u_control_input(2));
    u_control_input << 0, 0, 1;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Previous computation ----------------

  // Desired quadcopter attitude
  b_3_des[0]          = u_control_input[0] / u_control_input.norm();
  b_3_des[1]          = u_control_input[1] / u_control_input.norm();
  b_3_des[2]          = u_control_input[2] / u_control_input.norm();

  b_1_c[0]            = cosf(commanded_yaw_angle);
  b_1_c[1]            = sinf(commanded_yaw_angle);
  b_1_c[2]            = 0.0;

  b_2_des             = b_3_des.cross(b_1_c);
  b_1_des             = b_2_des.cross(b_3_des);

  R_des <<  b_1_des[0], b_2_des[0], b_3_des[0],
                      b_1_des[1], b_2_des[1], b_3_des[1],
                      b_1_des[2], b_2_des[2], b_3_des[2];

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Just try it ----------------

  // const Eigen::Vector3d fd_norm = u_control_input.normalized();

  // // | ------------------------- body z ------------------------- |
  // R_des.col(2) = fd_norm;

  // // | ------------------------- body x ------------------------- |

  // // construct the oblique projection
  // Eigen::Matrix3d projector_body_z_compl = (Eigen::Matrix3d::Identity(3, 3) - fd_norm * fd_norm.transpose());

  // // create a basis of the body-z complement subspace
  // Eigen::MatrixXd A_Matrix = Eigen::MatrixXd(3, 2);
  // A_Matrix.col(0)          = projector_body_z_compl.col(0);
  // A_Matrix.col(1)          = projector_body_z_compl.col(1);

  // // create the basis of the projection null-space complement
  // Eigen::MatrixXd B_Matrix = Eigen::MatrixXd(3, 2);
  // B_Matrix.col(0)          = Eigen::Vector3d(1, 0, 0);
  // B_Matrix.col(1)          = Eigen::Vector3d(0, 1, 0);

  // // oblique projector to <range_basis>
  // Eigen::MatrixXd Bt_A               = B_Matrix.transpose() * A_Matrix;
  // Eigen::MatrixXd Bt_A_pseudoinverse = ((Bt_A.transpose() * Bt_A).inverse()) * Bt_A.transpose();
  // Eigen::MatrixXd oblique_projector  = A_Matrix * Bt_A_pseudoinverse * B_Matrix.transpose();

  // R_des.col(0) = oblique_projector * b_1_c;
  // R_des.col(0).normalize();

  // // | ------------------------- body y ------------------------- |

  // R_des.col(1) = R_des.col(2).cross(R_des.col(0));
  // R_des.col(1).normalize();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  desired_thrust_force     = u_control_input.dot(R_curr.col(2));

  double throttle          = 0.0;

  if (desired_thrust_force >= 0) {
    throttle = mrs_lib::quadratic_throttle_model::forceToThrottle(common_handlers_->throttle_model, desired_thrust_force);
  } else {
    ROS_WARN_THROTTLE(1.0, "[ExampleController]: just so you know, the desired throttle force is negative (%.2f)", desired_thrust_force);
  }

  Eigen::Vector3d des_roll_pitch_yaw = Rotation_matrix_to_Euler_angle(R_des);

  des_roll_pitch_yaw(0) = clipping_angle(0.7,des_roll_pitch_yaw(0));
  des_roll_pitch_yaw(1) = clipping_angle(0.7,des_roll_pitch_yaw(1));

  // ROS_INFO("Roll des: %2.2f, Pitch des: %2.2f, Yaw des: %2.2f", des_roll_pitch_yaw(0), des_roll_pitch_yaw(1), des_roll_pitch_yaw(2));

  mrs_msgs::HwApiAttitudeCmd attitude_cmd;
  attitude_cmd.stamp       = ros::Time::now();
  attitude_cmd.orientation = mrs_lib::AttitudeConverter(des_roll_pitch_yaw[0],des_roll_pitch_yaw[1],des_roll_pitch_yaw[2]);
  attitude_cmd.throttle    = throttle;

  // | ----------------- set the control output ----------------- |

  last_control_output_.control_output = attitude_cmd;

  // | --------------- fill in the optional parts --------------- |

  //// it is recommended to fill the optinal parts if you know them

  /// this is used for:
  // * plotting the orientation in the control_refence topic (optional)
  // * checking for attitude control error
  // last_control_output_.desired_orientation = mrs_lib::AttitudeConverter(R_des);
  last_control_output_.desired_orientation = mrs_lib::AttitudeConverter(des_roll_pitch_yaw[0],des_roll_pitch_yaw[1],des_roll_pitch_yaw[2]);

  /// IMPORANT
  // The acceleration and heading rate in 3D (expressed in the "fcu" frame of reference) that the UAV will actually undergo due to the control action.
  last_control_output_.desired_unbiased_acceleration = Eigen::Vector3d(0, 0, 0);
  last_control_output_.desired_heading_rate          = 0;

  // | ----------------- fill in the diagnostics ---------------- |

  last_control_output_.diagnostics.controller = "ExampleController";

  // | ----------------- Return last control input ---------- |
  return last_control_output_;
}

//}

/* //{ getStatus() */

const mrs_msgs::ControllerStatus ExampleController::getStatus() {

  mrs_msgs::ControllerStatus controller_status;

  controller_status.active = is_active_;

  return controller_status;
}

//}

/* switchOdometrySource() //{ */

void ExampleController::switchOdometrySource([[maybe_unused]] const mrs_msgs::UavState& new_uav_state) {
}

//}

/* resetDisturbanceEstimators() //{ */

void ExampleController::resetDisturbanceEstimators(void) {
}

//}

/* setConstraints() //{ */

const mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr ExampleController::setConstraints([
    [maybe_unused]] const mrs_msgs::DynamicsConstraintsSrvRequest::ConstPtr& constraints) {

  if (!is_initialized_) {
    return mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr(new mrs_msgs::DynamicsConstraintsSrvResponse());
  }

  mrs_lib::set_mutexed(mutex_constraints_, constraints->constraints, constraints_);

  ROS_INFO("[ExampleController]: updating constraints");

  mrs_msgs::DynamicsConstraintsSrvResponse res;
  res.success = true;
  res.message = "constraints updated";

  return mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr(new mrs_msgs::DynamicsConstraintsSrvResponse(res));
}

//}

// --------------------------------------------------------------
// |                          callbacks                         |
// --------------------------------------------------------------

/* //{ callbackDrs() */

// void ExampleController::callback_quad_state(const nav_msgs::Odometry& msg) {

//   pos_of_quad[0] = msg.pose.pose.position.x;
//   pos_of_quad[1] = msg.pose.pose.position.y;
//   pos_of_quad[2] = msg.pose.pose.position.z;
//   // ROS_INFO_STREAM_THROTTLE(0.5, "From the subscriber" << pos_of_quad);

// }

// void ExampleController::callback_user_reference(const nav_msgs::Odometry& msg) {

//   des_pos_of_quad(0) = msg.pose.pose.position.x;
//   des_pos_of_quad(1) = msg.pose.pose.position.y;
//   des_pos_of_quad(2) = msg.pose.pose.position.z;

//   des_vel_of_quad(0) = msg.twist.twist.linear.x;
//   des_vel_of_quad(1) = msg.twist.twist.linear.y;
//   des_vel_of_quad(2) = msg.twist.twist.linear.z;

//     ROS_INFO("x des: %2.2f, y des: %2.2f, z des: %2.2f", des_pos_of_quad(0), des_pos_of_quad(1), des_pos_of_quad(2));

// }



void ExampleController::callbackDrs(example_controller_plugin::example_controllerConfig& config, [[maybe_unused]] uint32_t level) {

  mrs_lib::set_mutexed(mutex_drs_params_, config, drs_params_);

  ROS_INFO("[ExampleController]: dynamic reconfigure params updated");
}

void ExampleController::callback_cable_states(const nav_msgs::Odometry::ConstPtr msg) {

  q[0] = msg->pose.pose.position.x;
  q[1] = msg->pose.pose.position.y;
  q[2] = msg->pose.pose.position.z;

  q_dot[0] = msg->twist.twist.linear.x;
  q_dot[1] = msg->twist.twist.linear.y;
  q_dot[2] = msg->twist.twist.linear.z;

}

float ExampleController::min_acc_first_coefficient(float t1, float t2, float st, float en){
  float a = (2*(en - st)) /  ((t1 - t2)*(t1 - t2)*(t1 - t2));
  return a;
}

float ExampleController::min_acc_second_coefficient(float t1, float t2, float st, float en){
  float b = -(3*(t1 + t2)*(en - st)) /  ((t1 - t2)*(t1 - t2)*(t1 - t2));
  return b;
}

float ExampleController::min_acc_third_coefficient(float t1, float t2, float st, float en){
  float c = (6*t1*t2*(en - st)) /  ((t1 - t2)*(t1 - t2)*(t1 - t2));
  return c;
}

float ExampleController::min_acc_fourth_coefficient(float t1, float t2, float st, float en){
  float d = (en*t1*t1*t1 - 3*en*t1*t1*t2 + 3*st*t1*t2*t2 - st*t2*t2*t2) /  ((t1 - t2)*(t1 - t2)*(t1 - t2));
  return d;
}

float ExampleController::clipping_angle(float max_value, float current_angle){
  if (current_angle > max_value) // 0.35 rad means 20 deg
  {
    current_angle = max_value;
  }

  if (current_angle < -max_value) // 0.35 rad means 20 deg
  {
    current_angle = -max_value;
  }
  return current_angle;
}

float ExampleController::clipping_net_thrust_force(float max_value, float current_thrust){
  if (current_thrust > max_value) // 24.959870582 * 4
  {
    current_thrust = max_value;
  }

  if (current_thrust < 0.0) // 24.959870582 * 4
  {
    current_thrust = 0.0;
  }
  return current_thrust;
}

Eigen::Vector3d ExampleController::clipping_e_x_q(Eigen::Vector3d e_x_q_vector){
  float max_error_lim = 5.0; // in meter
  for (int i=0;i<=2;i++){
    if (e_x_q_vector(i) > max_error_lim ){
      e_x_q_vector(i) = max_error_lim;
    }
    if (e_x_q_vector(i) < -max_error_lim ){
      e_x_q_vector(i) = -max_error_lim;
    }
  }
return e_x_q_vector;
}

Eigen::Vector3d ExampleController::clipping_e_x_q_dot(Eigen::Vector3d e_x_q_dot_vector){
  float max_error_lim = 5.0; // in meter per second
  for (int i=0;i<=2;i++){
    if (e_x_q_dot_vector(i) > max_error_lim ){
      e_x_q_dot_vector(i) = max_error_lim;
    }
    if (e_x_q_dot_vector(i) < -max_error_lim ){
      e_x_q_dot_vector(i) = -max_error_lim;
    }
  }
return e_x_q_dot_vector;
}

float ExampleController::distance_bt_two_pts(Eigen::Vector3d A, Eigen::Vector3d B){

  float norm__  = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
  norm__ = sqrtf(norm__ );
  return norm__;

}

Eigen::Vector3d ExampleController::Matrix_vector_mul(Eigen::Matrix3d R, Eigen::Vector3d v){
  Eigen::Vector3d mul_vector (R(0,0)*v[0] + R(0,1)*v[1] + R(0,2)*v[2], R(1,0)*v[0] + R(1,1)*v[1] + R(1,2)*v[2],  R(2,0)*v[0] + R(2,1)*v[1] + R(2,2)*v[2]);
  return mul_vector;
}

Eigen::Vector3d ExampleController::Rotation_matrix_to_Euler_angle(Eigen::Matrix3d R){

  float sy      = sqrtf( R(0,0) * R(0,0) +  R(1,0) * R(1,0) );
  float ph_des  = atan2f(R(2,1) , R(2,2));
  float th_des  = atan2f(-R(2,0), sy);
  float ps_des  = atan2f(R(1,0), R(0,0));

  Eigen::Vector3d des_roll_pitch_yaw(ph_des, th_des, ps_des);
  return des_roll_pitch_yaw;
}

// void ExampleController::Publisher_QuadState(const mrs_msgs::UavState& uav_state){

//   if (!is_initialized_) {
//     return;
//   }

//   geometry_msgs::Pose current_pose;
//   current_pose.position         = uav_state.pose.position;
//   current_pose.orientation      = uav_state.pose.orientation;

//   try {
//     pub_quad_state_.publish(current_pose);
//   }
//   catch (...) {
//     ROS_ERROR("Exception caught during publishing topic %s.", pub_quad_state_.getTopic().c_str());
//   }

//   ros::spinOnce();
//   // loop_rate.sleep();

// }

/* timerPublishQuadState() //{ */

// void ExampleController::timerPublishQuadState([[maybe_unused]] const ros::TimerEvent& te, const mrs_msgs::UavState& uav_state) {

//   if (!is_initialized_) {
//     return;
//   }

//   geometry_msgs::Pose current_pose;
//   current_pose.position         = uav_state.pose.position;
//   current_pose.orientation      = uav_state.pose.orientation;

//   try {
//     pub_quad_state_.publish(current_pose);
//   }
//   catch (...) {
//     ROS_ERROR("Exception caught during publishing topic %s.", pub_quad_state_.getTopic().c_str());
//   }
// }

//}

}  // namespace example_controller

}  // namespace example_controller_plugin

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(example_controller_plugin::example_controller::ExampleController, mrs_uav_managers::Controller)
