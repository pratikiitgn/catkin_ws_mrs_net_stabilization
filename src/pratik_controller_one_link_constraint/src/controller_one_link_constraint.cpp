/* includes //{ */

#include <ros/ros.h>
#include <ros/package.h>
#include <nodelet/nodelet.h>

/* for loading dynamic parameters while the nodelet is running */
#include <mrs_uav_managers/controller.h>
#include <dynamic_reconfigure/server.h>

#include <pratik_controller_one_link_constraint/controller_one_link_constraintConfig.h>

#include <mrs_lib/param_loader.h>
#include <mrs_lib/subscribe_handler.h>
#include <mrs_lib/mutex.h>
#include <mrs_lib/attitude_converter.h>
#include <mrs_lib/msg_extractor.h>
#include <mrs_lib/geometry/misc.h>
#include <mrs_lib/profiler.h>
#include <mrs_lib/utils.h>
#include <mrs_lib/geometry/cyclic.h>
#include <mrs_lib/publisher_handler.h>

#include <mrs_msgs/ControlManagerDiagnostics.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

// | ----------------- Calling required libraries ----------------- |
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

/* for storing information about the state of the uav (position) */
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <nav_msgs/Odometry.h>

// | ----------------- Basic math variables ----------------- |
Eigen::Vector3d e1(1.0,0.0,0.0);
Eigen::Vector3d e2(0.0,1.0,0.0);
Eigen::Vector3d e3(0.0,0.0,1.0);

// | ----------------- System parameters     ---------------- |
float g_acceleration    = 9.81;     // in m/s^2
float PI_value          = 3.1415926535;

Eigen::Vector3d b_1_c(1.0,0.0,0.0);

Eigen::Vector3d des_rpy;

// | ----------------- Payload position State ----------------- |

Eigen::Vector2d Iw_w_   = Eigen::Vector2d(0, 0);
Eigen::Vector2d Ib_w    = Eigen::Vector2d(0, 0);
Eigen::Vector2d Ib_b_   = Eigen::Vector2d(0, 0);

// | ------------ controller limits and saturations ----------- |

bool   _tilt_angle_failsafe_enabled_;
double _tilt_angle_failsafe_;

double _throttle_saturation_;

// | ----------------- Link attitude State ----------------- |

double alpha                  = 0.0;
double alpha_dot              = 0.0;

Eigen::Vector3d     qb_link (0.0,0.0,-1.0);
Eigen::Vector3d qb_dot_link (0.0,0.0, 0.0);

Eigen::Vector3d     q_link (0.0,0.0,-1.0);
Eigen::Vector3d q_dot_link (0.0,0.0, 0.0);

// | ----------------- Desired link attitude State ----------------- |

double alpha_des      = 0.0;
double alpha_dot_des  = 0.0;

Eigen::Vector3d     q_d (0.0,0.0,-1.0);
Eigen::Vector3d q_d_dot (0.0,0.0, 0.0);

// | ----------------- Error in link attitude State ----------------- |

double e_alpha                = 0.0;
double e_alpha_dot            = 0.0;
Eigen::Vector3d     e_q       (0.0,0.0,0.0);
Eigen::Vector3d     e_q_dot   (0.0,0.0,0.0);

//  Simple bounded adaptive law for CAC gains
// ---------- persistent states (initialize once) ----------
static Eigen::Array3d k_p;
static Eigen::Array3d k_d;
static Eigen::Vector3d filt_e_q     = Eigen::Vector3d::Zero();
static Eigen::Vector3d filt_e_qdot  = Eigen::Vector3d::Zero();

// | ------------------------ profiler_ ------------------------ |

mrs_lib::Profiler profiler_;
bool              _profiler_enabled_ = false;

//}

namespace pratik_controller_one_link_constraint
{

namespace controller_one_link_constraint
{

typedef struct
{
  double kpx;           // Position gain in x
  double kpy;           // Position gain in y
  double kpz;           // Position gain in z
  double kdx;           // Velocity gain in x
  double kdy;           // Velocity gain in y
  double kdz;           // Velocity gain in z
  double km;            // mass estimator gain
  double km_lim;        // mass estimator limit
  double kq_roll_pitch; // pitch/roll attitude gain
  double kq_yaw;        // yaw attitude gain
  double kq_1;          // link attitude gains
  double kq_2;          // link attitude gains
  double kq_3;          // link attitude gains
  double kq_dot_1;      // link attitude gains
  double kq_dot_2;      // link attitude gains
  double kq_dot_3;      // link attitude gains
  double gamma_p_LAC;   // link attitude adaptation rate
  double gamma_d_LAC;   // link attitude adaptation rate
  double kappa_p_LAC;   // link attitude leakage rate
  double kappa_d_LAC;   // link attitude leakage rate
} Gains_t;

/* //{ class ControllerOneLinkConstraint */

class ControllerOneLinkConstraint : public mrs_uav_managers::Controller {

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
  float clipping_angle(float max_value, float current_angle);
  Eigen::Vector3d Matrix_vector_mul(Eigen::Matrix3d R, Eigen::Vector3d v);
  float clipping_net_thrust_force(float max_value, float current_thrust);
  Eigen::Vector3d clipping_e_x_q(Eigen::Vector3d e_x_q_vector);
  Eigen::Vector3d clipping_e_x_q_dot(Eigen::Vector3d e_x_q_dot_vector);
  Eigen::Vector3d Rotation_matrix_to_Euler_angle(Eigen::Matrix3d R);
  Eigen::Matrix3d so3transform(const Eigen::Vector3d& body_z, const ::Eigen::Vector3d& heading, const bool& preserve_heading);
  double getHeadingSafely(const mrs_msgs::UavState& uav_state, const mrs_msgs::TrackerCommand& tracker_command);
  std::optional<mrs_msgs::HwApiAttitudeRateCmd> attitudeController(const mrs_msgs::UavState& uav_state, const mrs_msgs::HwApiAttitudeCmd& reference,
                                                                 const Eigen::Vector3d& ff_rate, const Eigen::Vector3d& rate_saturation,
                                                                 const Eigen::Vector3d& gains, const bool& parasitic_heading_rate_compensation);
  Eigen::Vector3d orientationError(const Eigen::Matrix3d& R, const Eigen::Matrix3d& Rd);
  Eigen::Matrix3d hatmap(const Eigen::Vector3d &v);
  Eigen::VectorXd integrate_q_expmap_RK4(Eigen::Vector3d &q,
                            Eigen::Vector3d &q_dot,
                            const Eigen::Vector3d &q_ddot_des,
                            double dt);
  inline Eigen::Vector3d omega_dot(const Eigen::Vector3d &q,
                                 const Eigen::Vector3d &q_ddot_des,
                                 const Eigen::Vector3d &omega);
  inline Eigen::Vector3d expmap_rotate(const Eigen::Vector3d &theta,
                                     const Eigen::Vector3d &v);
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
  typedef pratik_controller_one_link_constraint::controller_one_link_constraintConfig DrsConfig_t;
  typedef dynamic_reconfigure::Server<DrsConfig_t>            Drs_t;
  boost::shared_ptr<Drs_t>                                    drs_;
  void                                                        callbackDrs(pratik_controller_one_link_constraint::controller_one_link_constraintConfig& config, uint32_t level);
  DrsConfig_t                                                 drs_params_;
  std::mutex                                                  mutex_drs_params_;

  // | ----------------------- constraints ---------------------- |

  mrs_msgs::DynamicsConstraints constraints_;
  std::mutex                    mutex_constraints_;

  // | --------- throttle generation and mass estimation -------- |

  double _uav_mass_;
  double uav_mass_difference_;

  Gains_t gains_;
  std::mutex mutex_gains_;       // locks the gains the are used and filtered

  ros::Timer timer_gains_;
  void       timerGains(const ros::TimerEvent& event);

  double _gain_filtering_rate_;

  double calculateGainChange(const double dt, const double current_value, const double desired_value, const bool bypass_rate, std::string name, bool& updated);

  double _gains_filter_change_rate_;
  double _gains_filter_min_change_rate_;

  // | ----------------------- gain muting ---------------------- |

  std::atomic<bool> mute_gains_            = false;
  std::atomic<bool> mute_gains_by_tracker_ = false;
  double            _gain_mute_coefficient_;

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
  mrs_lib::PublisherHandler<nav_msgs::Odometry> ph_CAC_gains_values;
  mrs_lib::PublisherHandler<visualization_msgs::Marker> ph_des_q_link;


  // | ---------------------- msg callbacks --------------------- |
  mrs_lib::SubscribeHandler<nav_msgs::Odometry>                  sh_link_states;
  void              callback_link_states(const nav_msgs::Odometry::ConstPtr msg);

};

//}



// --------------------------------------------------------------
// |                   controller's interface                   |
// --------------------------------------------------------------

/* //{ initialize() */

bool ControllerOneLinkConstraint::initialize(const ros::NodeHandle& nh, std::shared_ptr<mrs_uav_managers::control_manager::CommonHandlers_t> common_handlers,
                                   std::shared_ptr<mrs_uav_managers::control_manager::PrivateHandlers_t> private_handlers) {

  nh_ = nh;

  common_handlers_  = common_handlers;
  private_handlers_ = private_handlers;

  _uav_mass_ = common_handlers->getMass();

  last_update_time_ = ros::Time(0);

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

  success *= private_handlers->loadConfigFile(ros::package::getPath("pratik_controller_one_link_constraint") + "/config/controller_one_link_constraint.yaml");

  if (!success) {
    return false;
  }

  mrs_lib::ParamLoader param_loader(nh_, "ControllerOneLinkConstraint");

  const std::string yaml_namespace = "pratik_controller_one_link_constraint/default_gains/";

  private_handlers->param_loader->loadParam(yaml_namespace + "kpx", gains_.kpx);
  private_handlers->param_loader->loadParam(yaml_namespace + "kpy", gains_.kpy);
  private_handlers->param_loader->loadParam(yaml_namespace + "kpz", gains_.kpz);
  private_handlers->param_loader->loadParam(yaml_namespace + "kdx", gains_.kdx);
  private_handlers->param_loader->loadParam(yaml_namespace + "kdy", gains_.kdy);
  private_handlers->param_loader->loadParam(yaml_namespace + "kdz", gains_.kdz);

  private_handlers->param_loader->loadParam(yaml_namespace + "kq_1", gains_.kq_1);
  private_handlers->param_loader->loadParam(yaml_namespace + "kq_2", gains_.kq_2);
  private_handlers->param_loader->loadParam(yaml_namespace + "kq_3", gains_.kq_3);
  private_handlers->param_loader->loadParam(yaml_namespace + "kq_dot_1", gains_.kq_dot_1);
  private_handlers->param_loader->loadParam(yaml_namespace + "kq_dot_2", gains_.kq_dot_2);
  private_handlers->param_loader->loadParam(yaml_namespace + "kq_dot_3", gains_.kq_dot_3);

  // LAC - Link attitude controller auto gain tuner, parameter handles
  const std::string yaml_namespace_LAC_gain_tuner = "pratik_controller_one_link_constraint/LAC_gain_tuner/";
  private_handlers->param_loader->loadParam(yaml_namespace_LAC_gain_tuner + "gamma_p_LAC", gains_.gamma_p_LAC);
  private_handlers->param_loader->loadParam(yaml_namespace_LAC_gain_tuner + "gamma_d_LAC", gains_.gamma_d_LAC);
  private_handlers->param_loader->loadParam(yaml_namespace_LAC_gain_tuner + "kappa_p_LAC", gains_.kappa_p_LAC);
  private_handlers->param_loader->loadParam(yaml_namespace_LAC_gain_tuner + "kappa_d_LAC", gains_.kappa_d_LAC);
  private_handlers->param_loader->loadParam(yaml_namespace_LAC_gain_tuner + "LAC_Auto_Gain_Tuner_enabled",
                                            drs_params_.LAC_Auto_Gain_Tuner_enabled);

  ///////////////////

  private_handlers->param_loader->loadParam(yaml_namespace + "mass_estimator/km", gains_.km);
  private_handlers->param_loader->loadParam(yaml_namespace + "mass_estimator/km_lim", gains_.km_lim);

  // attitude gains
  private_handlers->param_loader->loadParam(yaml_namespace + "attitude/kq_roll_pitch", gains_.kq_roll_pitch);
  private_handlers->param_loader->loadParam(yaml_namespace + "attitude/kq_yaw", gains_.kq_yaw);

  // gain filtering
  const std::string yaml_namespace_gain = "pratik_controller_one_link_constraint/gain_filtering/";

  private_handlers->param_loader->loadParam(yaml_namespace_gain + "perc_change_rate", _gains_filter_change_rate_);
  private_handlers->param_loader->loadParam(yaml_namespace_gain + "min_change_rate", _gains_filter_min_change_rate_);
  private_handlers->param_loader->loadParam(yaml_namespace_gain + "rate", _gain_filtering_rate_);
  private_handlers->param_loader->loadParam(yaml_namespace_gain + "gain_mute_coefficient", _gain_mute_coefficient_);

  // angular rate feed forward
  const std::string yaml_namespace_ang_rate_ff = "pratik_controller_one_link_constraint/angular_rate_feedforward/";
  private_handlers->param_loader->loadParam(yaml_namespace_ang_rate_ff + "parasitic_pitch_roll",
                                            drs_params_.pitch_roll_heading_rate_compensation);
  private_handlers->param_loader->loadParam(yaml_namespace_ang_rate_ff + "jerk", drs_params_.jerk_feedforward);
  
  // constraints
  private_handlers->param_loader->loadParam("pratik_controller_one_link_constraint/constraints/tilt_angle_failsafe/enabled", _tilt_angle_failsafe_enabled_);
  private_handlers->param_loader->loadParam("pratik_controller_one_link_constraint/constraints/tilt_angle_failsafe/limit", _tilt_angle_failsafe_);

  _tilt_angle_failsafe_ = M_PI * (_tilt_angle_failsafe_ / 180.0);

  if (_tilt_angle_failsafe_enabled_ && std::abs(_tilt_angle_failsafe_) < 1e-3) {
    ROS_ERROR("[ControllerOneLinkConstraint]: constraints/tilt_angle_failsafe/enabled = 'TRUE' but the limit is too low");
    return false;
  }

  // | -------- initialize a publisher -------- |

  // | ------------------ finish loading params ----------------- |

  if (!param_loader.loadedSuccessfully()) {
    ROS_ERROR("[ControllerOneLinkConstraint]: could not load all parameters!");
    return false;
  }

  // | --------------- dynamic reconfigure server --------------- |
  drs_params_.kpx         = gains_.kpx;
  drs_params_.kpy         = gains_.kpy;
  drs_params_.kpz         = gains_.kpz;
  drs_params_.kdx         = gains_.kdx;
  drs_params_.kdy         = gains_.kdy;
  drs_params_.kdz         = gains_.kdz;

  drs_params_.kq_1          = gains_.kq_1;
  drs_params_.kq_2          = gains_.kq_2;
  drs_params_.kq_3          = gains_.kq_3;
  drs_params_.kq_dot_1      = gains_.kq_dot_1;
  drs_params_.kq_dot_2      = gains_.kq_dot_2;
  drs_params_.kq_dot_3      = gains_.kq_dot_3;

  drs_params_.gamma_p_LAC       = gains_.gamma_p_LAC;
  drs_params_.gamma_d_LAC       = gains_.gamma_d_LAC;
  drs_params_.kappa_p_LAC       = gains_.kappa_p_LAC;
  drs_params_.kappa_d_LAC       = gains_.kappa_d_LAC;

  drs_params_.km                = gains_.km;
  drs_params_.km_lim            = gains_.km_lim;
  drs_params_.kq_roll_pitch     = gains_.kq_roll_pitch;
  drs_params_.kq_yaw            = gains_.kq_yaw;
  drs_params_.jerk_feedforward  = true;

  drs_.reset(new Drs_t(mutex_drs_, nh_));
  drs_->updateConfig(drs_params_);
  Drs_t::CallbackType f = boost::bind(&ControllerOneLinkConstraint::callbackDrs, this, _1, _2);
  drs_->setCallback(f);

  // | ------------------------- timers ------------------------- |

  timer_gains_ = nh_.createTimer(ros::Rate(_gain_filtering_rate_), &ControllerOneLinkConstraint::timerGains, this, false, false);

  // | ------------------ initialize subscribers ----------------- |

  mrs_lib::SubscribeHandlerOptions shopts;
  shopts.nh                 = nh;
  shopts.node_name          = "ControllerOneLinkConstraint";
  shopts.no_message_timeout = ros::Duration(1.0);
  shopts.threadsafe         = true;
  shopts.autostart          = true;
  shopts.queue_size         = 10;
  shopts.transport_hints    = ros::TransportHints().tcpNoDelay();

  sh_link_states           = mrs_lib::SubscribeHandler<nav_msgs::Odometry>(shopts, "/sim_one_link_constraint/uav1/link_state",
                                                                                            &ControllerOneLinkConstraint::callback_link_states, this);

  ph_CAC_gains_values         = mrs_lib::PublisherHandler<nav_msgs::Odometry>(nh_, "/link_attitude_gains", 250, false);
  ph_des_q_link               = mrs_lib::PublisherHandler<visualization_msgs::Marker>(nh_, "/des_q_link", 50, false);

  // initialize the integrals
  uav_mass_difference_ = 0;

  // Simple bounded adaptive law for CAC gains
  
  k_p   << gains_.kq_1,     gains_.kq_2,     gains_.kq_3;
  k_d   << gains_.kq_dot_1, gains_.kq_dot_2, gains_.kq_dot_3;

  // | ----------------------- finish init ---------------------- |

  ROS_INFO("[ControllerOneLinkConstraint]: initialized");

  is_initialized_ = true;

  return true;
}

//}

/* //{ activate() */

bool ControllerOneLinkConstraint::activate(const ControlOutput& last_control_output) {

  activation_control_output_ = last_control_output;

  double activation_mass = _uav_mass_;

  if (activation_control_output_.diagnostics.mass_estimator) {
    uav_mass_difference_ = activation_control_output_.diagnostics.mass_difference;
    activation_mass += uav_mass_difference_;
    ROS_INFO("[ControllerOneLinkConstraint]: setting mass difference from the last control output: %.2f kg", uav_mass_difference_);
  }

  last_control_output_.diagnostics.controller_enforcing_constraints = false;

  first_iteration_ = true;
  mute_gains_      = true;

  timer_gains_.start();

  is_active_ = true;

  ROS_INFO("[ControllerOneLinkConstraint]: activated");

  return true;
}

//}

/* //{ deactivate() */

void ControllerOneLinkConstraint::deactivate(void) {

  is_active_            = false;
  first_iteration_      = false;
  uav_mass_difference_  = 0;
  timer_gains_.stop();

  ROS_INFO("[ControllerOneLinkConstraint]: deactivated");
}

//}

/* updateInactive() //{ */

void ControllerOneLinkConstraint::updateInactive(const mrs_msgs::UavState& uav_state, [[maybe_unused]] const std::optional<mrs_msgs::TrackerCommand>& tracker_command) {

  mrs_lib::set_mutexed(mutex_uav_state_, uav_state, uav_state_);

  last_update_time_ = uav_state.header.stamp;

  first_iteration_ = false;
}

//}

/* //{ updateActive() */

ControllerOneLinkConstraint::ControlOutput ControllerOneLinkConstraint::updateActive(const mrs_msgs::UavState& uav_state, const mrs_msgs::TrackerCommand& tracker_command) {

  auto drs_params = mrs_lib::get_mutexed(mutex_drs_params_, drs_params_);
  auto gains       = mrs_lib::get_mutexed(mutex_gains_, gains_);

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
    dt               = 0.01;
    // dt               = 0.004;
    first_iteration_ = false;
  } else {
    dt = (uav_state.header.stamp - last_update_time_).toSec();
  }

  last_update_time_ = uav_state.header.stamp;

  if (fabs(dt) < 0.001) {

    ROS_DEBUG("[ControllerOneLinkConstraint]: the last odometry message came too close (%.2f s)!", dt);
    dt = 0.01;
    // dt               = 0.004;
  }

  // | -------- check for the available output modalities ------- |

  // you can decide what to return, but it needs to be available
  if (common_handlers_->control_output_modalities.attitude) {
    ROS_INFO_THROTTLE(1.0, "[ControllerOneLinkConstraint]: desired attitude output modality is available");
  }

  // | ---------- extract the detailed model parameters --------- |

  if (common_handlers_->detailed_model_params) {

    mrs_uav_managers::control_manager::DetailedModelParams_t detailed_model_params = common_handlers_->detailed_model_params.value();

    // ROS_INFO_STREAM_THROTTLE(1.0, "[ControllerOneLinkConstraint]: UAV inertia is: " << detailed_model_params.inertia);
  }

  // | -------------- prepare the control reference ------------- |

  // --------------------------------------------------------------
  // |          load the control reference and estimates          |
  // --------------------------------------------------------------

  // Rp - position reference in global frame
  // Rv - velocity reference in global frame
  // Ra - velocity reference in global frame

  Eigen::Vector3d Rp = Eigen::Vector3d::Zero(3);
  Eigen::Vector3d Rv = Eigen::Vector3d::Zero(3);
  Eigen::Vector3d Ra = Eigen::Vector3d::Zero(3);

  if (tracker_command.use_position_vertical || tracker_command.use_position_horizontal) {

    if (tracker_command.use_position_horizontal) {
      Rp(0) = tracker_command.position.x;
      Rp(1) = tracker_command.position.y;
    } else {
      Rv(0) = 0;
      Rv(1) = 0;
    }

    if (tracker_command.use_position_vertical) {
      Rp(2) = tracker_command.position.z;
    } else {
      Rv(2) = 0;
    }
  }

  if (tracker_command.use_velocity_horizontal) {
    Rv(0) = tracker_command.velocity.x;
    Rv(1) = tracker_command.velocity.y;
  } else {
    Rv(0) = 0;
    Rv(1) = 0;
  }

  if (tracker_command.use_velocity_vertical) {
    Rv(2) = tracker_command.velocity.z;
  } else {
    Rv(2) = 0;
  }

  if (tracker_command.use_acceleration) {
    Ra << tracker_command.acceleration.x, tracker_command.acceleration.y, tracker_command.acceleration.z;
  } else {
    Ra << 0, 0, 0;
  }

  // | ----------------- Quadcopter State ----------------- |
  // Op - position in global frame
  // Ov - velocity in global frame

  Eigen::Vector3d   Op(uav_state.pose.position.x,uav_state.pose.position.y,uav_state.pose.position.z);
  Eigen::Vector3d   Ov(uav_state.velocity.linear.x,uav_state.velocity.linear.y,uav_state.velocity.linear.z);
  
  Eigen::Matrix3d R     = mrs_lib::AttitudeConverter(uav_state.pose.orientation);
  Eigen::Vector3d Omega (uav_state.velocity.angular.x,uav_state.velocity.angular.y,uav_state.velocity.angular.z);

  double uav_heading = getHeadingSafely(uav_state, tracker_command);

  // | ----------------- Error definition ----------------- |
  Eigen::Vector3d     Ep(0, 0, 0);
    // if (tracker_command.use_position_horizontal || tracker_command.use_position_vertical) {
    Ep = Rp - Op;
  // }

  Eigen::Vector3d     Ev(0, 0, 0);

  // if (tracker_command.use_velocity_horizontal || tracker_command.use_velocity_vertical ||
  //     tracker_command.use_position_vertical) {  // use_position_vertical = true, not a mistake, this provides dampening
    Ev = Rv - Ov;
  // }

  // | --------------------- load the gains --------------------- |

  mute_gains_by_tracker_ = tracker_command.disable_position_gains;

  Eigen::Array3d  Kp(0, 0, 0);
  Eigen::Array3d  Kv(0, 0, 0);
  Eigen::Array3d  Kq(0, 0, 0);
  Eigen::Array3d  kq_link(0, 0, 0);
  Eigen::Array3d  kq_dot_link(0, 0, 0);

  Kp(0) = gains.kpx;
  Kp(1) = gains.kpy;
  Kp(2) = gains.kpz;
  Kv(0) = gains.kdx;
  Kv(1) = gains.kdy;
  Kv(2) = gains.kdz;

  kq_link(0) = gains.kq_1;
  kq_link(1) = gains.kq_2;
  kq_link(2) = gains.kq_3;

  kq_dot_link(0) = gains.kq_dot_1;
  kq_dot_link(1) = gains.kq_dot_2;
  kq_dot_link(2) = gains.kq_dot_3;

  // Kp = Kp * (_uav_mass_ + uav_mass_difference_);
  // Kv = Kv * (_uav_mass_ + uav_mass_difference_);
  // | --------------------- Quadcopter position controller --------------------- |

  double total_mass                 = _uav_mass_ + uav_mass_difference_;
  Eigen::Vector3d feed_forward      = total_mass * (Eigen::Vector3d(0, 0, common_handlers_->g));
  Eigen::Vector3d position_feedback = Kp * Ep.array();
  Eigen::Vector3d velocity_feedback = Kv * Ev.array();

  // ROS_INFO_STREAM_THROTTLE(0.2,"Pos Feedback:" << position_feedback.transpose());
  // ROS_INFO_STREAM_THROTTLE(0.2,"Vel Feedback:" << velocity_feedback.transpose());

  Eigen::Vector3d u_quad_input      = position_feedback + velocity_feedback + feed_forward;
  
  /* mass estimatior //{ */

  // --------------------------------------------------------------
  // |                integrate the mass difference               |
  // --------------------------------------------------------------

  {
    std::scoped_lock lock(mutex_gains_);

    if (tracker_command.use_position_vertical && !rampup_active_) {
      uav_mass_difference_ += gains.km * Ep(2) * dt;
    }

    // saturate the mass estimator
    bool uav_mass_saturated = false;
    if (!std::isfinite(uav_mass_difference_)) {
      uav_mass_difference_ = 0;
      ROS_WARN_THROTTLE(1.0, "[Se3Controller]: NaN detected in variable 'uav_mass_difference_', setting it to 0 and returning!!!");
    } else if (uav_mass_difference_ > gains.km_lim) {
      uav_mass_difference_ = gains.km_lim;
      uav_mass_saturated   = true;
    } else if (uav_mass_difference_ < -gains.km_lim) {
      uav_mass_difference_ = -gains.km_lim;
      uav_mass_saturated   = true;
    }

    if (uav_mass_saturated) {
      ROS_WARN_THROTTLE(1.0, "[Se3Controller]: The UAV mass difference is being saturated to %.2f!", uav_mass_difference_);
    }
  }

  // ROS_INFO("x err: %2.2f, y err: %2.2f, z err: %2.2f", Ep(0), Ep(1), Ep(2));
  ROS_INFO_THROTTLE(0.5, "Mass Estimator: %2.2f", total_mass);

  // | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
  // | ---------------- Get the states of the link --------------- |

  qb_link(0)        = -sin(alpha);
  qb_link(1)        = 0.0;
  qb_link(2)        = -cos(alpha);
  
  qb_dot_link(0)    = -cos(alpha) * alpha_dot;
  qb_dot_link(1)    = 0.0;
  qb_dot_link(2)    =  sin(alpha) * alpha_dot;

  q_link            = R * qb_link;
  q_dot_link        = R * hatmap(Omega) * qb_link + R * qb_dot_link;

  // q_link            = qb_link;
  // q_dot_link        = qb_dot_link;

  // if (Ep.norm() > 1.0)
  // {
/////////////////////////////////////////////////////////////////////////////////////// 
              // // ---- compute u_perp ----
              // Eigen::Matrix3d I             = Eigen::Matrix3d::Identity();
              // Eigen::Vector3d q_vec         = q_link; // assume q is unit (normalize if needed)
              // Eigen::Matrix3d P_perp        = I - q_vec * q_vec.transpose();

              // double gamma                  = 0.5 / 3.0;
              // double mq                     = 3.0;
              // double l                      = 0.4;
              // Eigen::Vector3d u_perp        = (1.0 / (1.0 + gamma)) * P_perp * ( total_mass * Ra + total_mass * 9.81 * e3 );

              // // ---- desired q double dot (eq. (1)) ----
              // Eigen::Vector3d q_dot_norm2   = Eigen::Vector3d::Zero(); // only need norm
              // double qdot_norm2             = q_dot_link.squaredNorm();

              // Eigen::Vector3d q_ddot_des    = - (1.0 / (mq * l)) * u_perp - qdot_norm2 * q_vec;
              // // ROS_INFO_STREAM(" before qd " << q_d);

              // Eigen::VectorXd des_q_link_state;
              // des_q_link_state  = integrate_q_expmap_RK4(q_d, q_d_dot, q_ddot_des, dt);
              // q_d               = des_q_link_state.segment<3>(0);
              // q_d_dot           = des_q_link_state.segment<3>(3);

              // // ROS_INFO_STREAM(" after qd " << q_d);

              // Eigen::Vector3d des_pos_of_load;
              // des_pos_of_load           = Op + l * q_d;

              // visualization_msgs::Marker marker;

              // marker.header.frame_id = "uav1/world_origin" ;
              // marker.header.stamp    = ros::Time::now();
              // marker.pose.position.x = des_pos_of_load(0);
              // marker.pose.position.y = des_pos_of_load(1);
              // marker.pose.position.z = des_pos_of_load(2);
              // marker.pose.orientation = mrs_lib::AttitudeConverter(0, 0, 0);
              // marker.color.a = 1;
              // marker.color.r = 1;
              // marker.color.g = 1;
              // marker.color.b = 0;
              // marker.scale.x = 0.2;
              // marker.scale.y = 0.2;
              // marker.scale.z = 0.2;
              // marker.type         = visualization_msgs::Marker::SPHERE;

              // marker.ns       = "des q link";
              // ph_des_q_link.publish(marker);

    // compute tilde q_dot and project to tangent at q_d
    // Eigen::Vector3d qdot_tilde    = q_dot_link + q_ddot_des * dt;
    // Eigen::Vector3d qd_dot        = qdot_tilde - q_d * (q_d.dot(qdot_tilde)); // (I - qd*qd^T) * qdot_tilde

    // ---- desired angular rate omega_d such that qd_dot = omega_d x q_d ----
    // Eigen::Vector3d omega_d = q_d.cross(qd_dot);

    // ---- optional: desired angular acceleration omega_dot_d ----
    // Eigen::Vector3d r = q_ddot_des - omega_d.cross(omega_d.cross(q_d)); // r should be perpendicular to q_d
    // Eigen::Vector3d omega_dot_d = q_d.cross(r);

  /////////////////////////////////////////////////////////////////////////////////////// 

      // q_d           = u_quad_input.normalized();
    // Eigen::Matrix3d Rd_from_tracker = mrs_lib::AttitudeConverter(tracker_command.orientation);
    // q_d = -Rd_from_tracker.col(3);
    // q_d = -Rd_from_tracker.col(3);

  // }else{
  q_d         << 0.0,0.0,-1.0;
  q_d_dot     << 0.0,0.0,0.0;
  // }
  
  e_q               = q_link.cross(q_link.cross(q_d));
  e_q_dot           = q_dot_link - (q_d.cross(q_d_dot)).cross(q_link);

  // | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
  // | ---------------- prepare for link Attitude Controller --------------- |

  // false  - no autotune
  // true   - Simple bounded adaptive law

  Eigen::Vector3d u_link_input   (0.0,0.0,0.0);

  if (drs_params_.LAC_Auto_Gain_Tuner_enabled      == true)
  {
    // ---------- parameters (tune these) ----------

    Eigen::Array3d k_p0; 
    Eigen::Array3d k_d0;

    k_p0                            = kq_link;
    k_d0                            = kq_dot_link;

    Eigen::Array3d k_p_min; k_p_min << 0.0, 0.0, 0.0;
    Eigen::Array3d k_p_max; k_p_max << 50.0, 50.0, 50.0;
    Eigen::Array3d k_d_min; k_d_min << 0.0, 0.0, 0.0;
    Eigen::Array3d k_d_max; k_d_max << 50.0, 50.0, 50.0;

    Eigen::Array3d gamma_p; gamma_p << gains_.gamma_p_LAC, gains_.gamma_p_LAC, 0.0;   // adaptation rates
    Eigen::Array3d gamma_d; gamma_d << gains_.gamma_d_LAC, gains_.gamma_d_LAC, 0.0;
    Eigen::Array3d kappa_p; kappa_p << gains_.kappa_p_LAC, gains_.kappa_p_LAC, 0.0;   // leakage rates
    Eigen::Array3d kappa_d; kappa_d << gains_.kappa_d_LAC, gains_.kappa_d_LAC, 0.0;

    // filtering & rate limits
    double alpha_filt       = 0.85;       // smoother: 0.7-0.99
    double max_gain_rate    = 50.0;       // max change per second (elementwise cap)

    // ---------- inside control loop ----------
    // compute e_q (Vector3d) and e_q_dot (Vector3d) as you already do
    // e_q, e_q_dot are Eigen::Vector3d

    // low-pass filter the errors (sampled)
    filt_e_q            = alpha_filt * filt_e_q    + (1.0 - alpha_filt) * e_q;
    filt_e_qdot         = alpha_filt * filt_e_qdot + (1.0 - alpha_filt) * e_q_dot;

    // compute adaptation increments (elementwise)
    Eigen::Array3d inc_kp         = ( gamma_p * filt_e_q.array().square() )    - ( kappa_p * (k_p - k_p0) );
    Eigen::Array3d inc_kd         = ( gamma_d * filt_e_qdot.array().square() ) - ( kappa_d * (k_d - k_d0) );

    // discrete Euler update
    Eigen::Array3d k_p_new        = k_p + dt * inc_kp;
    Eigen::Array3d k_d_new        = k_d + dt * inc_kd;

    // rate limit (per-step) to avoid jumps
    Eigen::Array3d max_step = max_gain_rate * dt * Eigen::Array3d::Ones();
    for (int i=0;i<3;++i) {
        double delta_p = k_p_new[i] - k_p[i];
        if (delta_p > max_step[i]) k_p_new[i] = k_p[i] + max_step[i];
        if (delta_p < -max_step[i]) k_p_new[i] = k_p[i] - max_step[i];

        double delta_d = k_d_new[i] - k_d[i];
        if (delta_d > max_step[i]) k_d_new[i] = k_d[i] + max_step[i];
        if (delta_d < -max_step[i]) k_d_new[i] = k_d[i] - max_step[i];
    }

    // projection into bounds
    for (int i=0;i<3;++i) {
        if (k_p_new[i] < k_p_min[i]) k_p_new[i] = k_p_min[i];
        if (k_p_new[i] > k_p_max[i]) k_p_new[i] = k_p_max[i];
        if (k_d_new[i] < k_d_min[i]) k_d_new[i] = k_d_min[i];
        if (k_d_new[i] > k_d_max[i]) k_d_new[i] = k_d_max[i];
    }


    // commit updates
    k_p         = k_p_new;
    k_d         = k_d_new;

    // form control
    u_link_input = (k_p * e_q.array()).matrix() + (k_d * e_q_dot.array()).matrix();

    // Publish gains values to ROS topic
    nav_msgs::Odometry odom;

    odom.header.stamp    = ros::Time::now();

    odom.pose.pose.position.x = k_p(0);
    odom.pose.pose.position.y = k_p(1);
    odom.pose.pose.position.z = k_p(2);

    odom.twist.twist.linear.x = k_d(0);
    odom.twist.twist.linear.y = k_d(1);
    odom.twist.twist.linear.z = k_d(2);

    ph_CAC_gains_values.publish(odom);

  }
  else
  {
    // Eigen::Vector3d q_tangent_vec           = q_link.normalized();

    // Eigen::Array3d kq_link_custom           = kq_link     * q_link.array();
    // Eigen::Array3d kq_dot_link_custom       = kq_dot_link * q_link.array();

    u_link_input      = kq_link * e_q.array()  + kq_dot_link * e_q_dot.array();
  }

  // | -------------------------------------------------------------------------------------------------------------------------------------------------------- |

  // | ---------------- prepare the final control output --------------- |

  Eigen::Vector3d u_control_input;

  //  0 - standard
  //  1 - simple adaptive
  //  2 - Simple QP
  //  3 - Normlized error shaping

  int which_controller = 0;

  if (which_controller == 1)
  { 

    double eps          = 1e-6;
    double gamma        = 1.0;   // tune: >1 favors link, <1 favors quad
    double mu           = 1.0;   // relative importance of position

    double norm_eq      = e_q.norm();    // link attitude error
    double norm_ex      = Ep.norm();  // position error (choose e_pos or e_v)

    double lambda       = gamma * norm_eq / (norm_eq + mu*norm_ex + eps);
    lambda              = std::min(1.0, std::max(0.0, lambda));

    u_control_input = (1.0 - lambda) * u_quad_input + lambda * u_link_input;

    // ROS_INFO_STREAM_THROTTLE(0.5,"Lamda: " << lambda);

  }
  else if (which_controller == 2)
  {
      Eigen::Vector3d q_t;

      q_t(0)      = -cos(alpha);
      q_t(1)      = 0.0;
      q_t(2)      =  sin(alpha);

    Eigen::Vector3d q       = q_t.normalized();
    double alpha            = 10.0; // alpha large -> prioritize link, small -> prioritize quad
    double s                = q.dot(u_quad_input) + q.dot(u_link_input); // choice

    double factor           = alpha / (1.0 + alpha);
    double delta            = s - q.dot(u_quad_input);
    Eigen::Vector3d u_opt   = u_quad_input + factor * delta * q;
    u_control_input         = u_opt;

  }
  else if (which_controller == 3)
  {
    double gamma_p      = 1.0;
    double gamma_alpha  = 1.0;

    double norm_eq      = e_q.norm();    // link attitude error
    double norm_ex      = Ep.norm();  // position error (choose e_pos or e_v)

    double w_p          = 1.0 + gamma_p * (norm_ex*norm_ex) / (1.0 + norm_ex*norm_ex);
    double w_alpha      = 1.0 + gamma_alpha * (norm_eq*norm_eq) / (1.0 + norm_eq*norm_eq);

    u_control_input     = w_p * u_quad_input + w_alpha * u_link_input;

  }
  else if (which_controller == 0)
  {
    u_control_input     = u_quad_input + u_link_input;
  }

  if (u_control_input(2) < 0) {
    ROS_WARN_THROTTLE(1.0, "[ControllerOneLinkConstraint]: the calculated downwards desired force is negative (%.2f) -> mitigating flip", u_control_input(2));
    u_control_input << 0, 0, 1;
  }

  // | ------------------- desired orientation ------------------ |

  Eigen::Matrix3d Rd;

  if (tracker_command.use_orientation) {

    // fill in the desired orientation based on the desired orientation from the control command
    Rd = mrs_lib::AttitudeConverter(tracker_command.orientation);

    if (tracker_command.use_heading) {
      try {
        Rd = mrs_lib::AttitudeConverter(Rd).setHeading(tracker_command.heading);
      }
      catch (...) {
        ROS_ERROR_THROTTLE(1.0, "[ControllerOneLinkConstraint]: could not set the desired heading");
      }
    }

  } else {

    Eigen::Vector3d bxd;  // desired heading vector

    if (tracker_command.use_heading) {
      bxd << cos(tracker_command.heading), sin(tracker_command.heading), 0;
    } else {
      ROS_WARN_THROTTLE(10.0, "[ControllerOneLinkConstraint]: desired heading was not specified, using current heading instead!");
      bxd << cos(uav_heading), sin(uav_heading), 0;
    }

    Rd = so3transform(u_control_input, bxd, 1);
  }

  // | -------------------- desired throttle -------------------- |

  double desired_thrust_force     = u_control_input.dot(R.col(2));
  double throttle                 = 0.0;

  if (tracker_command.use_throttle) {

    // the throttle is overriden from the tracker command
    throttle = tracker_command.throttle;
  }
  else
  {
    if (desired_thrust_force >= 0) {
      throttle = mrs_lib::quadratic_throttle_model::forceToThrottle(common_handlers_->throttle_model, desired_thrust_force);
    } else {
      ROS_WARN_THROTTLE(1.0, "[ControllerOneLinkConstraint]: just so you know, the desired throttle force is negative (%.2f)", desired_thrust_force);
    }
  }

  Eigen::Vector3d unbiased_des_acc(0, 0, 0);

  {
    Eigen::Vector3d unbiased_des_acc_world = (position_feedback + velocity_feedback) / total_mass + Ra;

    geometry_msgs::Vector3Stamped world_accel;

    world_accel.header.stamp    = ros::Time::now();
    world_accel.header.frame_id = uav_state.header.frame_id;
    world_accel.vector.x        = unbiased_des_acc_world(0);
    world_accel.vector.y        = unbiased_des_acc_world(1);
    world_accel.vector.z        = unbiased_des_acc_world(2);

    auto res = common_handlers_->transformer->transformSingle(world_accel, "fcu");

    if (res) {
      unbiased_des_acc << res.value().vector.x, res.value().vector.y, res.value().vector.z;
    }
  }

  // | --------------- fill the resulting command --------------- |

  // fill the desired orientation for the tilt error check
  last_control_output_.desired_orientation = mrs_lib::AttitudeConverter(Rd);

  // fill the unbiased desired accelerations
  last_control_output_.desired_unbiased_acceleration = unbiased_des_acc;

  // | ----------------- fill in the diagnostics ---------------- |

  last_control_output_.diagnostics.ramping_up = rampup_active_;

  last_control_output_.diagnostics.mass_estimator  = true;
  last_control_output_.diagnostics.mass_difference = uav_mass_difference_;
  last_control_output_.diagnostics.total_mass      = total_mass;

  last_control_output_.diagnostics.disturbance_estimator = true;

  last_control_output_.diagnostics.disturbance_bx_b = -Ib_b_(0);
  last_control_output_.diagnostics.disturbance_by_b = -Ib_b_(1);

  last_control_output_.diagnostics.disturbance_bx_w = -Ib_w(0);
  last_control_output_.diagnostics.disturbance_by_w = -Ib_w(1);

  last_control_output_.diagnostics.disturbance_wx_w = -Iw_w_(0);
  last_control_output_.diagnostics.disturbance_wy_w = -Iw_w_(1);

  last_control_output_.diagnostics.controller_enforcing_constraints = false;

  last_control_output_.diagnostics.controller       = "ControllerOneLinkConstraint";

  // | ------------ construct the attitude reference ------------ |

  mrs_msgs::HwApiAttitudeCmd attitude_cmd;

  attitude_cmd.stamp       = ros::Time::now();
  attitude_cmd.orientation = mrs_lib::AttitudeConverter(Rd);
  attitude_cmd.throttle    = throttle;

  // | ----------------- set the control output ----------------- |
    last_control_output_.control_output = attitude_cmd;

  // | ----------------- Return last control input ---------- |
  return last_control_output_;

}

//}

/* orientationError() //{ */

Eigen::Vector3d ControllerOneLinkConstraint::orientationError(const Eigen::Matrix3d& R, const Eigen::Matrix3d& Rd) {

  // orientation error
  Eigen::Matrix3d R_error = 0.5 * (Rd.transpose() * R - R.transpose() * Rd);

  // vectorize the orientation error
  // clang-format off
    Eigen::Vector3d R_error_vec;
    R_error_vec << (R_error(1, 2) - R_error(2, 1)) / 2.0,
                   (R_error(2, 0) - R_error(0, 2)) / 2.0,
                   (R_error(0, 1) - R_error(1, 0)) / 2.0;
  // clang-format on

  return R_error_vec;
}

/* //{ getStatus() */

const mrs_msgs::ControllerStatus ControllerOneLinkConstraint::getStatus() {

  mrs_msgs::ControllerStatus controller_status;

  controller_status.active = is_active_;

  return controller_status;
}

//}

/* switchOdometrySource() //{ */

void ControllerOneLinkConstraint::switchOdometrySource([[maybe_unused]] const mrs_msgs::UavState& new_uav_state) {
}

//}

/* resetDisturbanceEstimators() //{ */

void ControllerOneLinkConstraint::resetDisturbanceEstimators(void) {
}

//}

/* setConstraints() //{ */

const mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr ControllerOneLinkConstraint::setConstraints([
    [maybe_unused]] const mrs_msgs::DynamicsConstraintsSrvRequest::ConstPtr& constraints) {

  if (!is_initialized_) {
    return mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr(new mrs_msgs::DynamicsConstraintsSrvResponse());
  }

  mrs_lib::set_mutexed(mutex_constraints_, constraints->constraints, constraints_);

  ROS_INFO("[ControllerOneLinkConstraint]: updating constraints");

  mrs_msgs::DynamicsConstraintsSrvResponse res;
  res.success = true;
  res.message = "constraints updated";

  return mrs_msgs::DynamicsConstraintsSrvResponse::ConstPtr(new mrs_msgs::DynamicsConstraintsSrvResponse(res));
}

//}

// --------------------------------------------------------------
// |                          callbacks                         |
// --------------------------------------------------------------

void ControllerOneLinkConstraint::callbackDrs(pratik_controller_one_link_constraint::controller_one_link_constraintConfig& config, [[maybe_unused]] uint32_t level) {

  mrs_lib::set_mutexed(mutex_drs_params_, config, drs_params_);

  ROS_INFO("[ControllerOneLinkConstraint]: dynamic reconfigure params updated");
}

// --------------------------------------------------------------
// |                           timers                           |
// --------------------------------------------------------------

/* timerGains() //{ */

void ControllerOneLinkConstraint::timerGains(const ros::TimerEvent& event) {

  mrs_lib::Routine    profiler_routine = profiler_.createRoutine("timerGains", _gain_filtering_rate_, 1.0, event);
  mrs_lib::ScopeTimer timer = mrs_lib::ScopeTimer("ControllerOneLinkConstraint::timerGains", common_handlers_->scope_timer.logger, common_handlers_->scope_timer.enabled);

  auto drs_params = mrs_lib::get_mutexed(mutex_drs_params_, drs_params_);
  auto gains      = mrs_lib::get_mutexed(mutex_gains_, gains_);

  // When muting the gains, we want to bypass the filter,
  // so it happens immediately.
  bool   bypass_filter = (mute_gains_ || mute_gains_by_tracker_);
  double gain_coeff    = (mute_gains_ || mute_gains_by_tracker_) ? _gain_mute_coefficient_ : 1.0;

  mute_gains_ = false;

  double dt = (event.current_real - event.last_real).toSec();

  if (!std::isfinite(dt) || (dt <= 0) || (dt > 5 * (1.0 / _gain_filtering_rate_))) {
    return;
  }

  bool updated = false;

  gains.kpx      = calculateGainChange(dt, gains.kpx, drs_params.kpx * gain_coeff, bypass_filter, "kpx", updated);
  gains.kpy      = calculateGainChange(dt, gains.kpy, drs_params.kpy * gain_coeff, bypass_filter, "kpy", updated);
  gains.kpz      = calculateGainChange(dt, gains.kpz, drs_params.kpz * gain_coeff, bypass_filter, "kpz", updated);
  gains.kdx      = calculateGainChange(dt, gains.kdx, drs_params.kdx * gain_coeff, bypass_filter, "kdx", updated);
  gains.kdy      = calculateGainChange(dt, gains.kdy, drs_params.kdy * gain_coeff, bypass_filter, "kdy", updated);
  gains.kdz      = calculateGainChange(dt, gains.kdz, drs_params.kdz * gain_coeff, bypass_filter, "kdz", updated);

  gains.kq_1     = calculateGainChange(dt, gains.kq_1, drs_params.kq_1 * gain_coeff, bypass_filter, "kq_1", updated);
  gains.kq_2     = calculateGainChange(dt, gains.kq_2, drs_params.kq_2 * gain_coeff, bypass_filter, "kq_2", updated);
  gains.kq_3     = calculateGainChange(dt, gains.kq_3, drs_params.kq_3 * gain_coeff, bypass_filter, "kq_3", updated);
  gains.kq_dot_1 = calculateGainChange(dt, gains.kq_dot_1, drs_params.kq_dot_1 * gain_coeff, bypass_filter, "kq_dot_1", updated);
  gains.kq_dot_2 = calculateGainChange(dt, gains.kq_dot_2, drs_params.kq_dot_2 * gain_coeff, bypass_filter, "kq_dot_2", updated);
  gains.kq_dot_3 = calculateGainChange(dt, gains.kq_dot_3, drs_params.kq_dot_3 * gain_coeff, bypass_filter, "kq_dot_3", updated);

  gains.gamma_p_LAC = calculateGainChange(dt, gains.gamma_p_LAC, drs_params.gamma_p_LAC * gain_coeff, bypass_filter, "gamma_p_LAC", updated);
  gains.gamma_d_LAC = calculateGainChange(dt, gains.gamma_d_LAC, drs_params.gamma_d_LAC * gain_coeff, bypass_filter, "gamma_d_LAC", updated);
  gains.kappa_p_LAC = calculateGainChange(dt, gains.kappa_p_LAC, drs_params.kappa_p_LAC * gain_coeff, bypass_filter, "kappa_p_LAC", updated);
  gains.kappa_d_LAC = calculateGainChange(dt, gains.kappa_d_LAC, drs_params.kappa_d_LAC * gain_coeff, bypass_filter, "kappa_d_LAC", updated);  

  gains.km       = calculateGainChange(dt, gains.km, drs_params.km * gain_coeff, bypass_filter, "km", updated);
  gains.kq_roll_pitch = calculateGainChange(dt, gains.kq_roll_pitch, drs_params.kq_roll_pitch * gain_coeff, bypass_filter, "kq_roll_pitch", updated);
  gains.kq_yaw        = calculateGainChange(dt, gains.kq_yaw, drs_params.kq_yaw * gain_coeff, bypass_filter, "kq_yaw", updated);

  // do not apply muting on these gains
  gains.km_lim   = calculateGainChange(dt, gains.km_lim, drs_params.km_lim, false, "km_lim", updated);

  mrs_lib::set_mutexed(mutex_gains_, gains, gains_);

  // set the gains back to dynamic reconfigure
  // and only do it when some filtering occurs
  if (updated) {

    drs_params.kpx      = gains.kpx;
    drs_params.kpy      = gains.kpy;
    drs_params.kpz      = gains.kpz;
    drs_params.kdx      = gains.kdx;
    drs_params.kdy      = gains.kdy;
    drs_params.kdz      = gains.kdz;
    drs_params.km       = gains.km;
    drs_params.km_lim   = gains.km_lim;

    drs_params.gamma_p_LAC  = gains.gamma_p_LAC;
    drs_params.gamma_d_LAC  = gains.gamma_d_LAC;
    drs_params.kappa_p_LAC  = gains.kappa_p_LAC;
    drs_params.kappa_d_LAC  = gains.kappa_d_LAC;

    drs_params.kq_1     = gains.kq_1;
    drs_params.kq_2     = gains.kq_2;
    drs_params.kq_3     = gains.kq_3;
    drs_params.kq_dot_1 = gains.kq_dot_1;
    drs_params.kq_dot_2 = gains.kq_dot_2;
    drs_params.kq_dot_3 = gains.kq_dot_3;

    drs_params.kq_roll_pitch = gains.kq_roll_pitch;
    drs_params.kq_yaw        = gains.kq_yaw;

    drs_->updateConfig(drs_params);

    ROS_INFO_THROTTLE(10.0, "[ControllerOneLinkConstraint]: gains have been updated");
  }
}

//}

/* attitudeController() //{ */

std::optional<mrs_msgs::HwApiAttitudeRateCmd> ControllerOneLinkConstraint::attitudeController(const mrs_msgs::UavState& uav_state, const mrs_msgs::HwApiAttitudeCmd& reference,
                                                                 const Eigen::Vector3d& ff_rate, const Eigen::Vector3d& rate_saturation,
                                                                 const Eigen::Vector3d& gains, const bool& parasitic_heading_rate_compensation) {

  Eigen::Matrix3d R  = mrs_lib::AttitudeConverter(uav_state.pose.orientation);
  Eigen::Matrix3d Rd = mrs_lib::AttitudeConverter(reference.orientation);

  // calculate the orientation error
  Eigen::Vector3d E = orientationError(R, Rd);

  Eigen::Vector3d rate_feedback = gains.array() * E.array() + ff_rate.array();

  // | ----------- parasitic heading rate compensation ---------- |

  if (parasitic_heading_rate_compensation) {

    ROS_DEBUG_THROTTLE(1.0, "[AttitudeController]: parasitic heading rate compensation enabled");

    // compensate for the parasitic heading rate created by the desired pitch and roll rate
    Eigen::Vector3d rp_heading_rate_compensation(0, 0, 0);

    Eigen::Vector3d q_feedback_yawless = rate_feedback;
    q_feedback_yawless(2)              = 0;  // nullyfy the effect of the original yaw feedback

    double parasitic_heading_rate = 0;

    try {
      parasitic_heading_rate = mrs_lib::AttitudeConverter(uav_state.pose.orientation).getHeadingRate(q_feedback_yawless);
    }
    catch (...) {
      ROS_ERROR("[AttitudeController]: exception caught while calculating the parasitic heading rate!");
    }

    try {
      rp_heading_rate_compensation(2) = mrs_lib::AttitudeConverter(uav_state.pose.orientation).getYawRateIntrinsic(-parasitic_heading_rate);
    }
    catch (...) {
      ROS_ERROR("[AttitudeController]: exception caught while calculating the parasitic heading rate compensation!");
    }

    rate_feedback += rp_heading_rate_compensation;
  }

  // | --------------- saturate the attitude rate --------------- |

  if (rate_feedback(0) > rate_saturation(0)) {
    rate_feedback(0) = rate_saturation(0);
  } else if (rate_feedback(0) < -rate_saturation(0)) {
    rate_feedback(0) = -rate_saturation(0);
  }

  if (rate_feedback(1) > rate_saturation(1)) {
    rate_feedback(1) = rate_saturation(1);
  } else if (rate_feedback(1) < -rate_saturation(1)) {
    rate_feedback(1) = -rate_saturation(1);
  }

  if (rate_feedback(2) > rate_saturation(2)) {
    rate_feedback(2) = rate_saturation(2);
  } else if (rate_feedback(2) < -rate_saturation(2)) {
    rate_feedback(2) = -rate_saturation(2);
  }

  // | ------------ fill in the attitude rate command ----------- |

  mrs_msgs::HwApiAttitudeRateCmd cmd;

  cmd.stamp = ros::Time::now();

  cmd.body_rate.x = rate_feedback(0);
  cmd.body_rate.y = rate_feedback(1);
  cmd.body_rate.z = rate_feedback(2);

  cmd.throttle = reference.throttle;

  return cmd;
}

void ControllerOneLinkConstraint::callback_link_states(const nav_msgs::Odometry::ConstPtr msg) {

  alpha                         = msg->pose.pose.position.x;
  // alpha_dot                     = msg->twist.twist.linear.x;

  double alpha_dot_filter       = 0.85;       // smoother: 0.7-0.99
  alpha_dot                     = alpha_dot_filter * alpha_dot    + (1.0 - alpha_dot_filter) * msg->twist.twist.linear.x;;

}

float ControllerOneLinkConstraint::clipping_angle(float max_value, float current_angle){
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

float ControllerOneLinkConstraint::clipping_net_thrust_force(float max_value, float current_thrust){
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

Eigen::Vector3d ControllerOneLinkConstraint::clipping_e_x_q(Eigen::Vector3d e_x_q_vector){
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

Eigen::Vector3d ControllerOneLinkConstraint::clipping_e_x_q_dot(Eigen::Vector3d e_x_q_dot_vector){
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

float ControllerOneLinkConstraint::distance_bt_two_pts(Eigen::Vector3d A, Eigen::Vector3d B){

  float norm__  = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
  norm__ = sqrtf(norm__ );
  return norm__;

}

Eigen::Vector3d ControllerOneLinkConstraint::Matrix_vector_mul(Eigen::Matrix3d R, Eigen::Vector3d v){
  Eigen::Vector3d mul_vector (R(0,0)*v[0] + R(0,1)*v[1] + R(0,2)*v[2], R(1,0)*v[0] + R(1,1)*v[1] + R(1,2)*v[2],  R(2,0)*v[0] + R(2,1)*v[1] + R(2,2)*v[2]);
  return mul_vector;
}

// helper: expmap rotation (Rodrigues)
inline Eigen::Vector3d ControllerOneLinkConstraint::expmap_rotate(const Eigen::Vector3d &theta,
                                     const Eigen::Vector3d &v) {
    double th = theta.norm();
    if (th < 1e-9) {
        // small angle expansion
        return v + hatmap(theta) * v + 0.5 * hatmap(theta) * hatmap(theta) * v;
    }
    Eigen::Vector3d k = theta / th;
    Eigen::Matrix3d K = hatmap(k);
    return v * cos(th) + K * v * sin(th) + k * (k.dot(v)) * (1.0 - cos(th));
}

// compute omega_dot from q, q_dot, omega, and q_ddot_des
inline Eigen::Vector3d ControllerOneLinkConstraint::omega_dot(const Eigen::Vector3d &q,
                                 const Eigen::Vector3d &q_ddot_des,
                                 const Eigen::Vector3d &omega) {
    Eigen::Vector3d term = omega.cross(omega.cross(q));
    return q.cross(q_ddot_des - term);
}

// RK4 integrator for q, omega
inline Eigen::VectorXd ControllerOneLinkConstraint::integrate_q_expmap_RK4(Eigen::Vector3d &q,
                            Eigen::Vector3d &q_dot,
                            const Eigen::Vector3d &q_ddot_des,
                            double dt) 
{
  
    // ensure unit vector
    q.normalize();

    // current angular velocity
    Eigen::Vector3d omega = q.cross(q_dot);

    // RK4 integration of omega_dot
    Eigen::Vector3d k1 = omega_dot(q, q_ddot_des, omega);
    
    Eigen::Vector3d k2 = omega_dot(q, q_ddot_des, omega + 0.5 * dt * k1);
    
    Eigen::Vector3d k3 = omega_dot(q, q_ddot_des, omega + 0.5 * dt * k2);
    
    Eigen::Vector3d k4 = omega_dot(q, q_ddot_des, omega + dt * k3);

    Eigen::Vector3d omega_new = omega + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);

    // average omega for rotation update (midpoint)
    Eigen::Vector3d omega_avg = 0.5 * (omega + omega_new);

    // update q using exponential map
    Eigen::Vector3d q_new = expmap_rotate(omega_avg * dt, q);
    q_new.normalize();

    // update q_dot consistent with omega_new
    Eigen::Vector3d q_dot_new = omega_new.cross(q_new);

    // assign back
    // q       = q_new;
    // q_dot   = q_dot_new;
    Eigen::VectorXd output(6);
    
    output.segment<3>(0) = q_new;
    
    output.segment<3>(3) = q_dot_new;
    
    return output;
}


Eigen::Vector3d ControllerOneLinkConstraint::Rotation_matrix_to_Euler_angle(Eigen::Matrix3d R){

  float sy      = sqrtf( R(0,0) * R(0,0) +  R(1,0) * R(1,0) );
  float ph_des  = atan2f(R(2,1) , R(2,2));
  float th_des  = atan2f(-R(2,0), sy);
  float ps_des  = atan2f(R(1,0), R(0,0));

  Eigen::Vector3d des_roll_pitch_yaw(ph_des, th_des, ps_des);
  return des_roll_pitch_yaw;
}

// hatmap map (skew-symmetric) for a 3-vector
Eigen::Matrix3d ControllerOneLinkConstraint::hatmap(const Eigen::Vector3d &v) {
  Eigen::Matrix3d H;
  H <<     0.0, -v.z(),  v.y(),
        v.z(),    0.0, -v.x(),
       -v.y(),  v.x(),    0.0;
  return H;
}

Eigen::Matrix3d ControllerOneLinkConstraint::so3transform(const Eigen::Vector3d& body_z, const ::Eigen::Vector3d& heading, const bool& preserve_heading) {

  Eigen::Vector3d body_z_normed = body_z.normalized();

  Eigen::Matrix3d Rd;

  if (preserve_heading) {

    ROS_DEBUG_THROTTLE(1.0, "[SO3Transform]: using Baca's method");

    // | ------------------------- body z ------------------------- |
    Rd.col(2) = body_z_normed;

    // | ------------------------- body x ------------------------- |

    // construct the oblique projection
    Eigen::Matrix3d projector_body_z_compl = (Eigen::Matrix3d::Identity(3, 3) - body_z_normed * body_z_normed.transpose());

    // create a basis of the body-z complement subspace
    Eigen::MatrixXd A = Eigen::MatrixXd(3, 2);
    A.col(0)          = projector_body_z_compl.col(0);
    A.col(1)          = projector_body_z_compl.col(1);

    // create the basis of the projection null-space complement
    Eigen::MatrixXd B = Eigen::MatrixXd(3, 2);
    B.col(0)          = Eigen::Vector3d(1, 0, 0);
    B.col(1)          = Eigen::Vector3d(0, 1, 0);

    // oblique projector to <range_basis>
    Eigen::MatrixXd Bt_A               = B.transpose() * A;
    Eigen::MatrixXd Bt_A_pseudoinverse = ((Bt_A.transpose() * Bt_A).inverse()) * Bt_A.transpose();
    Eigen::MatrixXd oblique_projector  = A * Bt_A_pseudoinverse * B.transpose();

    Rd.col(0) = oblique_projector * heading;
    Rd.col(0).normalize();

    // | ------------------------- body y ------------------------- |

    Rd.col(1) = Rd.col(2).cross(Rd.col(0));
    Rd.col(1).normalize();

  } else {

    ROS_DEBUG_THROTTLE(1.0, "[SO3Transform]: using Lee's method");

    Rd.col(2) = body_z_normed;
    Rd.col(1) = Rd.col(2).cross(heading);
    Rd.col(1).normalize();
    Rd.col(0) = Rd.col(1).cross(Rd.col(2));
    Rd.col(0).normalize();
  }

  return Rd;
}

/* calculateGainChange() //{ */

double ControllerOneLinkConstraint::calculateGainChange(const double dt, const double current_value, const double desired_value, const bool bypass_rate, std::string name,
                                          bool& updated) {

  double change = desired_value - current_value;

  double gains_filter_max_change = _gains_filter_change_rate_ * dt;
  double gains_filter_min_change = _gains_filter_min_change_rate_ * dt;

  if (!bypass_rate) {

    // if current value is near 0...
    double change_in_perc;
    double saturated_change;

    if (std::abs(current_value) < 1e-6) {
      change *= gains_filter_max_change;
    } else {

      saturated_change = change;

      change_in_perc = ((current_value + saturated_change) / current_value) - 1.0;

      if (change_in_perc > gains_filter_max_change) {
        saturated_change = current_value * gains_filter_max_change;
      } else if (change_in_perc < -gains_filter_max_change) {
        saturated_change = current_value * -gains_filter_max_change;
      }

      if (std::abs(saturated_change) < std::abs(change) * gains_filter_min_change) {
        change *= gains_filter_min_change;
      } else {
        change = saturated_change;
      }
    }
  }

  if (std::abs(change) > 1e-3) {
    ROS_DEBUG("[ControllerOneLinkConstraint]: changing gain '%s' from %.2f to %.2f", name.c_str(), current_value, desired_value);
    updated = true;
  }

  return current_value + change;
}

//}

double ControllerOneLinkConstraint::getHeadingSafely(const mrs_msgs::UavState& uav_state, const mrs_msgs::TrackerCommand& tracker_command) {

  try {
    return mrs_lib::AttitudeConverter(uav_state.pose.orientation).getHeading();
  }
  catch (...) {
  }

  try {
    return mrs_lib::AttitudeConverter(uav_state.pose.orientation).getYaw();
  }
  catch (...) {
  }

  if (tracker_command.use_heading) {
    return tracker_command.heading;
  }

  return 0;
}

// void ControllerOneLinkConstraint::Publisher_QuadState(const mrs_msgs::UavState& uav_state){

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

// void ControllerOneLinkConstraint::timerPublishQuadState([[maybe_unused]] const ros::TimerEvent& te, const mrs_msgs::UavState& uav_state) {

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

}  // namespace controller_one_link_constraint

}  // namespace pratik_controller_one_link_constraint

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(pratik_controller_one_link_constraint::controller_one_link_constraint::ControllerOneLinkConstraint, mrs_uav_managers::Controller)
