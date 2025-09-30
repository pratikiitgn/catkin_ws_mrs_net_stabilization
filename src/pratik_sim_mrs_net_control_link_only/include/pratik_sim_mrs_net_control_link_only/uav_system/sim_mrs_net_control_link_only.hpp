/**
 * @brief Dynamic simulation of a Quadcopter with MRS NET Mechanism
 *
 * Acknowledgement:
 * * https://github.com/HKUST-Aerial-Robotics/Fast-Planner
 */

#ifndef MULTIROTOR_MODEL_H
#define MULTIROTOR_MODEL_H

// Number of states
// xq               - 3
// xq_dot           - 3
// R                - 9
// omega            - 3
// alpha            - 1
// alpha_dot        - 1

// Total states     = 18 + 2 = 20

// #define N_INTERNAL_STATES 20

#include <boost/array.hpp>
#include "ode/boost/numeric/odeint.hpp"
#include "controllers/references.hpp"

#include <Eigen/Dense>

using namespace Eigen;

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

Eigen::Vector3d rho_vector(0.0,0.0,-0.01);    // Vertical offset of the net mechanism underneath the drone

constexpr int NO_OF_CAP           = 3;        // (m)
constexpr int N_INTERNAL_STATES   = 18 + 2;
constexpr int alpha_index_start   = 18;
constexpr int N_states_in_EOM     = 3 + 3 + 1; // total DOF (3 position + 3 attitude + 1 link)

Eigen::VectorXd X_dot_dot_custom(N_INTERNAL_STATES,1);
Eigen::VectorXd X_dot_dot(N_INTERNAL_STATES,1);
Eigen::VectorXd rhs(N_states_in_EOM,1);
Eigen::VectorXd X_ddot(N_INTERNAL_STATES,1);

// Other system parameters
// system parameters
constexpr float mr                    = 0.05;            // Mass of net rod
constexpr float lr                    = 1.0;             // half length of rod
constexpr float m_control_link        = 0.1;             // mass of control link
constexpr double l_control_link       = 0.4;             // Length of control link
double distance_between_each_CAP      = 2.0 * lr / (NO_OF_CAP - 1);
MatrixXd control_link_mass_array      = MatrixXd::Constant(NO_OF_CAP, 1, m_control_link);
double c1                             = (control_link_mass_array.sum()) * l_control_link;

std::vector<Eigen::Vector3d> rho_i_vectors(NO_OF_CAP);

Eigen::Vector3d e2(0.0,1.0,0.0);
Eigen::Vector3d e3(0.0,0.0,1.0);

namespace odeint = boost::numeric::odeint;

namespace pratik_sim_mrs_net_control_link_only
{

class MultirotorModel {

public:
  class ModelParams {
  public:
    ModelParams() {

      // | -------- default parameters of the x500 quadrotor -------- |
      // it is recommended to load them through the setParams() method

      n_motors             = 4;
      g                    = 9.81;
      mass                 = 2.0;
      kf                   = 0.00000027087;
      km                   = 0.07;
      prop_radius          = 0.15;
      arm_length           = 0.25;
      body_height          = 0.1;
      motor_time_constant  = 0.03;
      max_rpm              = 7800;
      min_rpm              = 1170;
      air_resistance_coeff = 0.30;

      // First Rigid Link Only
      mp                    = 0.1;           // mass of load
      l                     = 0.4;             // length of link

      J       = Eigen::Matrix3d::Zero();
      J(0, 0) = mass * (3.0 * arm_length * arm_length + body_height * body_height) / 12.0;
      J(1, 1) = mass * (3.0 * arm_length * arm_length + body_height * body_height) / 12.0;
      J(2, 2) = (mass * arm_length * arm_length) / 2.0;
      
      allocation_matrix = Eigen::MatrixXd::Zero(4, 4);

      // clang-format off
      allocation_matrix <<
      -0.707, 0.707, 0.707,  -0.707,
      -0.707, 0.707, -0.707, 0.707,
      -1,     -1,    1,      1,
      1,      1,     1,      1;
      // clang-format on

      allocation_matrix.row(0) *= arm_length * kf;
      allocation_matrix.row(1) *= arm_length * kf;
      allocation_matrix.row(2) *= km * (3.0 * prop_radius) * kf;
      allocation_matrix.row(3) *= kf;

      ground_enabled        = false;
      takeoff_patch_enabled = true;

    }

    int    n_motors;
    double g;
    double mass;
    double kf;
    double km;
    double prop_radius;
    double arm_length;
    double body_height;
    double motor_time_constant;
    double max_rpm;
    double min_rpm;
    double air_resistance_coeff;

    // First Rigid Link Only
    double mp;
    double l;

    Eigen::Matrix3d J;
    Eigen::MatrixXd allocation_matrix;

    bool   ground_enabled;
    double ground_z;

    bool takeoff_patch_enabled;
  };

  struct State
  {
    Eigen::Vector3d x;
    Eigen::Vector3d v;
    Eigen::Vector3d v_prev;
    Eigen::Matrix3d R;
    Eigen::Vector3d omega;
    Eigen::VectorXd motor_rpm;
    double alpha;
    double alpha_dot;

  };

  MultirotorModel();

  MultirotorModel(const ModelParams& params, const Eigen::Vector3d& spawn_pos, const double spawn_heading);

  const MultirotorModel::State& getState(void) const;

  void setState(const MultirotorModel::State& state);

  void applyForce(const Eigen::Vector3d& state);

  void setStatePos(const Eigen::Vector3d& pos, const double heading);

  const Eigen::Vector3d& getExternalForce(void) const;
  void                   setExternalForce(const Eigen::Vector3d& force);

  const Eigen::Vector3d& getExternalMoment(void) const;
  void                   setExternalMoment(const Eigen::Vector3d& moment);

  void setInput(const reference::Actuators& input);

  void step(const double& dt);

  typedef boost::array<double, N_INTERNAL_STATES> InternalState;

  void operator()(const MultirotorModel::InternalState& x, MultirotorModel::InternalState& dxdt, const double t);

  Eigen::Vector3d getImuAcceleration() const;

  inline Eigen::Matrix3d hatmap(const Eigen::Vector3d &v);
  inline double wrapToPi(double angle);

  double M_00(const MatrixXd &control_link_mass_array, double mq, double mr);
  Vector3d A1_vector(const MatrixXd &control_link_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, double l_control_link, double alpha);
  Matrix3d A2_matrix(double mr, double g, const Vector3d &rho, const Matrix3d &R, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m);
  Matrix3d A3_matrix(const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m);
  Vector3d A4_vector(const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, double alpha, const MatrixXd &control_link_mass_array, int m);
  Matrix3d J_bar(const Matrix3d &Jq, double mr, const Vector3d &rho, double lr, const Vector3d &e2, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m);

  Vector3d qb_vec(double alpha);
  Vector3d qt_vec(double alpha);

  VectorXd EOM_mrs_net_control_link_only(
  const VectorXd &XK,             // big state vector
  double mq,
  const Matrix3d &Jq,
  int m,                          // number of CAP
  const Vector3d &rho,
  const std::vector<Vector3d> &rho_i_vectors, // size m
  double lr,
  double g,
  double m_I,
  double intruder_position,
  int drone_hitting_i_th_cable_attachment_point,
  int drone_hitting_j_th_point_mass,
  const Vector3d &e3,
  double mr,
  const Vector3d &e2,
  const Eigen::Vector3d &M,
  double f, double l_control_link, const MatrixXd &control_link_mass_array);

  ModelParams getParams(void);
  void        setParams(const ModelParams& params);

  void initializeState(void);

private:

  void updateInternalState(void);

  MultirotorModel::State state_;

  Eigen::Vector3d imu_acceleration_;

  Eigen::VectorXd input_;
  Eigen::Vector3d external_force_;
  Eigen::Vector3d external_moment_;

  Eigen::Vector3d _initial_pos_;

  ModelParams params_;

  InternalState internal_state_;
};

/* constructor MultirotorModel //{ */

MultirotorModel::MultirotorModel(void) {

  initializeState();
  updateInternalState();

}

MultirotorModel::MultirotorModel(const MultirotorModel::ModelParams& params, const Eigen::Vector3d& spawn_pos, const double spawn_heading) {

  params_ = params;

  initializeState();

  _initial_pos_           = spawn_pos;
  state_.x                = spawn_pos;
  state_.R                = Eigen::AngleAxis(-spawn_heading, Eigen::Vector3d(0, 0, 1));
  
  state_.alpha            = 0.0;
  state_.alpha_dot        = 0.0;

  updateInternalState();
}

//}

/* initializeState() //{ */

void MultirotorModel::initializeState(void) {

  state_.x                = Eigen::Vector3d::Zero();
  state_.v                = Eigen::Vector3d::Zero();
  state_.v_prev           = Eigen::Vector3d::Zero();
  state_.R                = Eigen::Matrix3d::Identity();
  state_.omega            = Eigen::Vector3d::Zero();

  state_.alpha                                =  0.0;
  state_.alpha_dot                            =  0.0;

  imu_acceleration_       = Eigen::Vector3d::Zero();

  state_.motor_rpm        = Eigen::VectorXd::Zero(params_.n_motors);
  input_                  = Eigen::VectorXd::Zero(params_.n_motors);

  external_force_.setZero();
  external_moment_.setZero();
}

//}

/* updatedInternalState() //{ */

void MultirotorModel::updateInternalState(void) {

  for (int i = 0; i < 3; ++i)
  {
    internal_state_.at(0 + i)                 = state_.x(i);
    internal_state_.at(3 + i)                 = state_.v(i);
    internal_state_.at(6 + i)                 = state_.R(i, 0);
    internal_state_.at(9 + i)                 = state_.R(i, 1);
    internal_state_.at(12 + i)                = state_.R(i, 2);
    internal_state_.at(15 + i)                = state_.omega(i);
  }

  internal_state_.at(alpha_index_start)       = wrapToPi(state_.alpha);
  internal_state_.at(alpha_index_start + 1)   = state_.alpha_dot;

}

//}

/* step() //{ */

void MultirotorModel::step(const double& dt) {

  // ROS_INFO("First hi from step() inside .hpp file");

  auto save = internal_state_;

  // for (int i = 0; i < N_INTERNAL_STATES; ++i) {
  //   ROS_INFO("%d th -> Internal state from step() inside .hpp file %f",i, internal_state_.at(i));
  // }

  // boost::numeric::odeint::runge_kutta_fehlberg78<InternalState> rk;
  boost::numeric::odeint::runge_kutta4<InternalState> rk;
  
  // ROS_INFO("First hi from step() inside .hpp file");

  odeint::integrate_n_steps(rk, boost::ref(*this), internal_state_, 0.0, dt, 1);

  // ROS_INFO("N_INTERNAL_STATES -> %d from step() inside .hpp file", N_INTERNAL_STATES);

  for (int i = 0; i < N_INTERNAL_STATES; ++i) {
    if (std::isnan(internal_state_.at(i))) {
      internal_state_ = save;
      break;
    }
  }

  for (int i = 0; i < 3; ++i) {
    state_.x(i)     = internal_state_.at(0 + i);
    state_.v(i)     = internal_state_.at(3 + i);
    state_.R(i, 0)  = internal_state_.at(6 + i);
    state_.R(i, 1)  = internal_state_.at(9 + i);
    state_.R(i, 2)  = internal_state_.at(12 + i);
    state_.omega(i) = internal_state_.at(15 + i);
  }

  state_.alpha      = wrapToPi(internal_state_.at(alpha_index_start));
  state_.alpha_dot  = internal_state_.at(alpha_index_start + 1);

  // ROS_INFO("x(0)      from step() inside .hpp file %f", state_.x(0));
  // ROS_INFO("v(0)      from step() inside .hpp file %f", state_.v(0));
  // ROS_INFO("R (0, 0)  from step() inside .hpp file %f", state_.R(0,0));
  // ROS_INFO("omega (0) from step() inside .hpp file %f", state_.omega(0));
  // ROS_INFO("alpha     from step() inside .hpp file %f", state_.alpha);
  // ROS_INFO("alpha dot from step() inside .hpp file %f", state_.alpha_dot);

  double filter_const = exp((-dt) / (params_.motor_time_constant));

  state_.motor_rpm = filter_const * state_.motor_rpm + (1.0 - filter_const) * input_;

  Matrix3d Rt = state_.R.transpose();

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(Rt * state_.R);
  Eigen::Matrix3d P = llt.matrixL();
  Eigen::Matrix3d R = state_.R * P.inverse();
  state_.R          = R;

  // simulate the ground
  if (params_.ground_enabled) {
    if (state_.x(2) < params_.ground_z && state_.v(2) < 0) {
      state_.x(2)  = params_.ground_z;
      state_.v     = Eigen::Vector3d::Zero();
      state_.omega = Eigen::Vector3d::Zero();
    }
  }

  if (params_.takeoff_patch_enabled) {

    const double hover_rpm = sqrt((params_.mass * params_.g) / (params_.n_motors * params_.kf));
    if (input_.mean() <= 0.90 * hover_rpm) {

      if (state_.x(2) < _initial_pos_(2) && state_.v(2) < 0) {
        state_.x(2)  = _initial_pos_(2);
        state_.v     = Eigen::Vector3d::Zero();
        state_.omega = Eigen::Vector3d::Zero();
      }
    } else {
      params_.takeoff_patch_enabled = false;
    }
  }

  // fabricate what the onboard accelerometer would feel
  imu_acceleration_ = Rt * (((state_.v - state_.v_prev) / dt) + Eigen::Vector3d(0, 0, params_.g));
  state_.v_prev     = state_.v;

  // simulate the takeoff patch
  updateInternalState();
}

//}

/* applyModel() //{ */

void MultirotorModel::applyForce(const Eigen::Vector3d& force) {

  external_force_ = force;
}

//}

/* operator() //{ */

void MultirotorModel::operator()(const MultirotorModel::InternalState& x, MultirotorModel::InternalState& dxdt, [[maybe_unused]] const double t) {

  // ROS_INFO("First hi from operator()");

  // for (int i = 0; i < N_INTERNAL_STATES; ++i) {
  //   ROS_INFO("%d th -> x state from operator() inside .hpp file %f",i, x.at(i));
  // }

  State cur_state;

  for (int i = 0; i < 3; ++i) {
    cur_state.x(i)     = x.at(0 + i);
    cur_state.v(i)     = x.at(3 + i);
    cur_state.R(i, 0)  = x.at(6 + i);
    cur_state.R(i, 1)  = x.at(9 + i);
    cur_state.R(i, 2)  = x.at(12 + i);
    cur_state.omega(i) = x.at(15 + i);
  }

  cur_state.alpha       = wrapToPi(x.at(alpha_index_start));
  cur_state.alpha_dot   = x.at(alpha_index_start + 1);
  // ROS_INFO("Third hi from operator()");

  Matrix3d Rt           = cur_state.R.transpose();

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(Rt * cur_state.R);
  Eigen::Matrix3d             P = llt.matrixL();
  cur_state.R                   = cur_state.R * P.inverse();

  // ROS_INFO("Fourth hi from operator()");

  Eigen::VectorXd XK(N_INTERNAL_STATES,1);

  for (int i = 0; i < 3; i++) {
    XK(0 + i)   = cur_state.x(i);
    XK(3 + i)   = cur_state.v(i);
    XK(6 + i)   = cur_state.R(i, 0);
    XK(9 + i)   = cur_state.R(i, 1);
    XK(12 + i)  = cur_state.R(i, 2);
    XK(15 + i)  = cur_state.omega(i);
  }
  // ROS_INFO("Sixth hi from operator()");
  XK(alpha_index_start)         = cur_state.alpha;
  XK(alpha_index_start + 1)     = cur_state.alpha_dot;
  // ROS_INFO("Seven hi from operator()");

  // for (int i = 0; i < N_INTERNAL_STATES; ++i) {
  //   ROS_INFO("%d th -> XK state from operator() inside .hpp file %f",i, XK(i));
  // }

  Eigen::VectorXd motor_rpm_sq  = state_.motor_rpm.array().square();
  Eigen::Vector4d torque_thrust = params_.allocation_matrix * motor_rpm_sq;
  double          thrust        = torque_thrust(3);

  // double resistance       = params_.air_resistance_coeff * M_PI * (params_.arm_length) * (params_.arm_length) * cur_state.v.norm() * cur_state.v.norm();
  // Eigen::Vector3d vnorm   = cur_state.v;
  // if (vnorm.norm() != 0) {
  //   vnorm.normalize();
  // }

  for (int i_rho = 0; i_rho < NO_OF_CAP; ++i_rho) {
    double distance_along_second_axis   = lr - i_rho * distance_between_each_CAP;
    rho_i_vectors[i_rho]                = Vector3d(0.0, distance_along_second_axis, 0.0);
  }

  ROS_INFO_STREAM("rho_i_vectors:");
  for (size_t i = 0; i < rho_i_vectors.size(); ++i) {
      ROS_INFO_STREAM(" [" << i << "] = " << rho_i_vectors[i].transpose());
  }

  double m_I                                    = 0.5;        // Mass of intruder drone
  int drone_hitting_i_th_cable_attachment_point = 1;
  int drone_hitting_j_th_point_mass             = 1;
  double intruder_position                      = 0.5;

  X_dot_dot_custom = EOM_mrs_net_control_link_only(
  XK,             // current state vector
  params_.mass,
  params_.J,
  NO_OF_CAP,                    // number of cable attachment points
  rho_vector,
  rho_i_vectors, // size m
  lr,
  params_.g,
  m_I,
  intruder_position,
  drone_hitting_i_th_cable_attachment_point,
  drone_hitting_j_th_point_mass,
  e3,
  mr,
  e2,
  torque_thrust.topRows(3), thrust, l_control_link, control_link_mass_array);

  // v_dot = v_dot + external_force_ / params_.mass - resistance * vnorm / params_.mass;

  for (int i = 0; i < N_INTERNAL_STATES; ++i) {
    dxdt.at(i)  = X_dot_dot_custom(i);
    // dxdt.at(i)  = 0;
  }

  for (int i = 0; i < N_INTERNAL_STATES; ++i) {
    if (std::isnan(dxdt.at(i))) {
      dxdt.at(i) = 0;
    }
  }

}

//}

// | ------------------- setters and getters ------------------ |

/* setParams() //{ */

void MultirotorModel::setParams(const MultirotorModel::ModelParams& params) {

  params_ = params;
}

//}

/* getParams() //{ */

MultirotorModel::ModelParams MultirotorModel::getParams(void) {

  return params_;
}

//}

/* setInput() //{ */

void MultirotorModel::setInput(const reference::Actuators& input) {

  for (int i = 0; i < params_.n_motors; ++i) {

    double val = input.motors(i);

    if (!std::isfinite(val)) {
      val = 0;
    }

    if (val < 0.0) {
      val = 0.0;
    } else if (val > 1.0) {
      val = 1.0;
    }

    input_(i) = params_.min_rpm + (params_.max_rpm - params_.min_rpm) * val;
  }
}

//}

/* getState() //{ */

const MultirotorModel::State& MultirotorModel::getState(void) const {
  return state_;
}

//}

/* setState() //{ */

void MultirotorModel::setState(const MultirotorModel::State& state) {

  state_.x         = state.x;
  state_.v         = state.v;
  state_.R         = state.R;
  state_.omega     = state.omega;
  state_.motor_rpm = state.motor_rpm;

  updateInternalState();
}

//}

/* setStatePos() //{ */

void MultirotorModel::setStatePos(const Eigen::Vector3d& pos, const double heading) {

  _initial_pos_           = pos;
  state_.x                = pos;
  state_.R                = Eigen::AngleAxis(-heading, Eigen::Vector3d(0, 0, 1));

  updateInternalState();
}

//}

/* getExternalForce() //{ */

const Eigen::Vector3d& MultirotorModel::getExternalForce(void) const {
  return external_force_;
}

//}

/* setExternalForce() //{ */

void MultirotorModel::setExternalForce(const Eigen::Vector3d& force) {
  external_force_ = force;
}

//}

/* getExternalMoment() //{ */

const Eigen::Vector3d& MultirotorModel::getExternalMoment(void) const {
  return external_moment_;
}

//}

/* setExternalMoment() //{ */

void MultirotorModel::setExternalMoment(const Eigen::Vector3d& moment) {
  external_moment_ = moment;
}

//}

/* getImuAcceleration() //{ */

Eigen::Vector3d MultirotorModel::getImuAcceleration() const {
  return imu_acceleration_;
}

//}
// hat map (skew-symmetric) for a 3-vector
Eigen::Matrix3d MultirotorModel::hatmap(const Eigen::Vector3d &v) {
  Eigen::Matrix3d H;
  H <<     0.0, -v.z(),  v.y(),
        v.z(),    0.0, -v.x(),
       -v.y(),  v.x(),    0.0;
  return H;
}

double MultirotorModel::wrapToPi(double angle) {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0)
        angle += 2.0 * M_PI;
    return angle - M_PI;
}

// qb_vec and qt_vec as used in MATLAB version
Vector3d MultirotorModel::qb_vec(double alpha) {
  return Eigen::Vector3d(-std::sin(alpha), 0.0, -std::cos(alpha));
}

Vector3d MultirotorModel::qt_vec(double alpha) {
  return Eigen::Vector3d(-std::cos(alpha), 0.0, std::sin(alpha));
}

double MultirotorModel::M_00(const MatrixXd &control_link_mass_array, double mq, double mr)
{
  return mq + mr + control_link_mass_array.sum();
}

Vector3d MultirotorModel::A1_vector(const MatrixXd &control_link_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, double l_control_link, double alpha)
{ 
  Vector3d out = mr * rho;
  for (int i = 1; i <= m; ++i) {
      out += control_link_mass_array(i-1, 0) * (rho + rho_i_vectors[i-1] + l_control_link * qb_vec(alpha));
  }
  return out;
}

Matrix3d MultirotorModel::A2_matrix(double mr, double g, const Vector3d &rho, const Matrix3d &R, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m)
{
  Matrix3d Rt = R.transpose();

  Matrix3d out = mr * g * hatmap(rho) * Rt;
  for (int i = 1; i <= m; ++i) {
      Matrix3d temp_matrix = hatmap(rho + rho_i_vectors[i-1] + l_control_link * qb_vec(alpha));
      out += g * control_link_mass_array(i-1, 0) * temp_matrix * Rt;
  }
  return out;
}

Matrix3d MultirotorModel::A3_matrix(const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m)
{
  Matrix3d out;
  for (int i = 1; i <= m; ++i) {
      Matrix3d temp_matrix = hatmap(rho + rho_i_vectors[i-1] + l_control_link * qb_vec(alpha));
      out += control_link_mass_array(i-1, 0) * temp_matrix * l_control_link;
  }
  return out;
}

Vector3d MultirotorModel::A4_vector(const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, double alpha, const MatrixXd &control_link_mass_array, int m)
{
  Vector3d out;
  for (int i = 1; i <= m; ++i) {
      Vector3d temp_matrix = rho + rho_i_vectors[i-1] + l_control_link * qb_vec(alpha);
      out += control_link_mass_array(i-1, 0) * temp_matrix;
  }
  return out;
}

Matrix3d MultirotorModel::J_bar(const Matrix3d &Jq, double mr, const Vector3d &rho, double lr, const Vector3d &e2, const std::vector<Vector3d> &rho_i_vectors, double l_control_link, double alpha, const MatrixXd &control_link_mass_array, int m)
{
  Matrix3d Jbar    = Jq;
  Jbar            -= 0.5 * mr * hatmap(rho) * hatmap(rho);
  Jbar            -= (1.0/6.0) * mr * lr * lr * hatmap(e2) * hatmap(e2);

  for (int i = 1; i <= m; ++i) {
      Matrix3d h = hatmap(rho + rho_i_vectors[i-1] + l_control_link * qb_vec(alpha));
      Jbar -= 0.5 *control_link_mass_array(i-1,0) * h * h;
  }
  return Jbar;
}

// ---------- EOM main function ----------
VectorXd MultirotorModel::EOM_mrs_net_control_link_only(
  const VectorXd &XK,             // big state vector
  double mq,
  const Matrix3d &Jq,
  int m,                          // number of branches
  const Vector3d &rho,
  const std::vector<Vector3d> &rho_i_vectors, // size m
  double lr,
  double g,
  double m_I,
  double intruder_position,
  int drone_hitting_i_th_cable_attachment_point,
  int drone_hitting_j_th_point_mass,
  const Vector3d &e3,
  double mr,
  const Vector3d &e2,
  const Eigen::Vector3d &M,
  double f, double l_control_link, const MatrixXd &control_link_mass_array)
  {
  // Basic checks
  assert(rho_i_vectors.size()       == static_cast<size_t>(m));

  // ROS_INFO("Tenth hi from EOM()");

  // safety
  if (XK.size() < 18) {
    std::cerr << "XK length too small: " << XK.size() << "\n";
    return VectorXd::Zero(1);
  }

  // for (int i = 0; i < N_INTERNAL_STATES; ++i) {
  //   ROS_INFO("%d th -> XK state from EOM() inside .hpp file %f",i, XK(i));
  // }

  Vector3d xq;
  xq      << XK(0), XK(1), XK(2);

  Vector3d xq_dot;
  xq_dot  << XK(3), XK(4), XK(5);

  Matrix3d R;
  R << XK(6),  XK(9),  XK(12),
       XK(7), XK(10),  XK(13),
       XK(8), XK(11),  XK(14);

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(R.transpose() * R);
  Eigen::Matrix3d P = llt.matrixL();
  R = R * P.inverse();

  // ROS_INFO("Eleven hi from EOM()");

  Vector3d omega;
  omega   << XK(15), XK(16), XK(17);

  Matrix3d R_dot = R * hatmap(omega);

  // ROS_INFO("15 hi from EOM()");

  // Control link state variables
  double alpha_in_EOM             = 0.0;
  double alpha_dot_in_EOM         = 0.0;

  alpha_in_EOM                    = XK(alpha_index_start);
  alpha_dot_in_EOM                = XK(alpha_index_start + 1);

  Vector3d u                      = f * R * e3;

  // ---------- General Defination ----------
  Matrix3d Rt                     = R.transpose();
  Eigen::Matrix3d omega_hat       = hatmap(omega);

  // ---------- Build Inertia Matrix ----------
  MatrixXd Final_Inertia = MatrixXd::Zero(N_states_in_EOM, N_states_in_EOM);

    // (1) Quadcotper Translational Inertia matrix
        MatrixXd Inertia_matrix_quad_position_raw                     = MatrixXd::Zero(3, N_states_in_EOM);

        // 1.1 --- (3x3) matrix
        Inertia_matrix_quad_position_raw.block(0, 0, 3, 3)            = Matrix3d::Identity() * M_00(control_link_mass_array, mq, mr);
        int col_offset                                                = 3;

        // 1.2 --- (3x3) matrix
        Vector3d A1_vector_ = A1_vector(control_link_mass_array, mr, rho, rho_i_vectors, m, l_control_link, alpha_in_EOM);
        Inertia_matrix_quad_position_raw.block(0, col_offset, 3, 3)   = -R * hatmap(A1_vector_);
        col_offset += 3;

        // 1.3 --- (3x1) matrix
        Vector3d last_vec_inertia_mat                                 = c1 * R * qt_vec(alpha_in_EOM);
        Inertia_matrix_quad_position_raw.block(0, col_offset, 3, 1)   = last_vec_inertia_mat;

        Final_Inertia.block(0, 0, 3, N_states_in_EOM)                 = Inertia_matrix_quad_position_raw;

    // (2) Quadcotper Attitude Inertia matrix
        MatrixXd quadcopter_attitude_inertia_matrix                   = MatrixXd::Zero(3, N_states_in_EOM);
        
        // 2.1 --- (3x3) matrix - quad pos
        quadcopter_attitude_inertia_matrix.block(0, 0, 3, 3)          = A2_matrix(mr, g, rho, R, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m);
        int col_index_quad_inertia_mat                                = 3;

        // 2.2 --- (3x3) matrix - quad att
        Matrix3d Jbar = J_bar(Jq, mr, rho, lr, e2, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m);
        quadcopter_attitude_inertia_matrix.block(0, col_index_quad_inertia_mat, 3, 3) = Jbar;
        col_index_quad_inertia_mat += 3; // move to next 3-column block

        // 2.3 --- (3x3) matrix - control link
        quadcopter_attitude_inertia_matrix.block(0, col_index_quad_inertia_mat, 3, 1) = A3_matrix(rho, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m) * qt_vec(alpha_in_EOM);

        Final_Inertia.block(3, 0, 3, N_states_in_EOM)                 = quadcopter_attitude_inertia_matrix;

    // (3) Control Link Inertia Matrix

        MatrixXd control_link_attitude_inertia_matrix = MatrixXd::Zero(1, N_states_in_EOM);

        // 3.1 --- (1x3) matrix - quad pos
        control_link_attitude_inertia_matrix.block(0, 0, 1, 3) = c1 * qt_vec(alpha_in_EOM).transpose() * Rt;
        int col_index_first_link_mat    = 3;

        // 3.2 --- (1x3) matrix - quad att
        control_link_attitude_inertia_matrix.block(0, col_index_first_link_mat, 1, 3) = -l_control_link * qt_vec(alpha_in_EOM).transpose() * hatmap(A4_vector(rho, rho_i_vectors, alpha_in_EOM, control_link_mass_array, m));
        col_index_first_link_mat       += 3;

        // 3.3 --- (1x1) matrix - control link
        control_link_attitude_inertia_matrix(0, col_index_first_link_mat) = c1 * l_control_link;
  
        Final_Inertia.block(6, 0, 1, N_states_in_EOM) = control_link_attitude_inertia_matrix;

  // ---------- Coriolis/centripeta and gravity term ----------
  VectorXd coriolis_final             = VectorXd::Zero(N_states_in_EOM);

    // (1) Quadcotper Translational Coriolis/centripeta and gravity term
      Vector3d coriolis_quad_pos      = R * omega_hat * omega_hat * A1_vector(control_link_mass_array, mr, rho, rho_i_vectors, m, l_control_link, alpha_in_EOM);
      coriolis_quad_pos              += 2 * c1 * alpha_dot_in_EOM * R * omega_hat * qt_vec(alpha_in_EOM);
      coriolis_quad_pos              +=    -c1 * alpha_dot_in_EOM * alpha_dot_in_EOM * R * qb_vec(alpha_in_EOM);
      coriolis_quad_pos              +=   M_00(control_link_mass_array, mq, mr) * g * e3;
      coriolis_final.segment(0, 3)    = coriolis_quad_pos;

    // (2) Quadcotper Attitude      Coriolis/centripeta and gravity term
        Matrix3d Jbar_for_quad_coriolis     = J_bar(Jq, mr, rho, lr, e2, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m);
        Vector3d coriolis_quad_att          = omega_hat * Jbar_for_quad_coriolis * omega;

        Matrix3d Jbar_dot;
        for (int j_J_bar_dor  = 1; j_J_bar_dor <= m; ++j_J_bar_dor) {
          Jbar_dot            += - (control_link_mass_array(j_J_bar_dor-1, 0)) * l_control_link * alpha_dot_in_EOM * hatmap(rho + rho_i_vectors[j_J_bar_dor-1] + l_control_link * qb_vec(alpha_in_EOM) ) * hatmap(qt_vec(alpha_in_EOM));
        }

        coriolis_quad_att   += Jbar_dot * omega;
        coriolis_quad_att   += - alpha_dot_in_EOM * alpha_dot_in_EOM * A3_matrix(rho, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m) * qb_vec(alpha_in_EOM);
        coriolis_quad_att   += alpha_dot_in_EOM * omega_hat * A3_matrix(rho, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m) * qt_vec(alpha_in_EOM);
        coriolis_quad_att   += A2_matrix(mr, g, rho, R, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m) * e3;

        coriolis_final.segment(3, 3) = coriolis_quad_att;

    // (3) Quadcotper Control Link  Coriolis/centripeta and gravity term
        double scalar1              = 0.5 * omega.dot(A3_matrix(rho, rho_i_vectors, l_control_link, alpha_in_EOM, control_link_mass_array, m) * hatmap(qt_vec(alpha_in_EOM)) * omega);
        double scalar2              = c1 * g * qt_vec(alpha_in_EOM).dot(Rt * e3);
        double coriolis_first_link  = scalar1 + scalar2;

        coriolis_final.segment(6, 1).setConstant(coriolis_first_link);

  // ROS_INFO_STREAM_THROTTLE(0.01,"Final_Inertia \n" << Final_Inertia);

  // ---------- Control input column ----------
      VectorXd control_input_column                   = VectorXd::Zero(N_states_in_EOM);
      control_input_column.segment(0, 3)              = u;
      control_input_column.segment(3, 3)              = M;

  // ---------- Solve for X_dot_dot ----------
  rhs                                   = control_input_column - coriolis_final;
  Eigen::PartialPivLU<MatrixXd> lu(Final_Inertia);
  X_dot_dot                             = lu.solve(rhs);

  // ROS_INFO_STREAM_THROTTLE(1,"X_dot_dot \n" << X_dot_dot.transpose());

  // VectorXd rhs       = VectorXd::Zero(N_states_in_EOM);

  // xq dot
  X_ddot.segment(0, 3)                  = xq_dot;
  // xq dot dot 
  X_ddot.segment(3, 3)                  = X_dot_dot.segment(0, 3);
  // R_dot
  X_ddot.segment(6, 3)                  = R_dot.col(0);
  X_ddot.segment(9, 3)                  = R_dot.col(1);
  X_ddot.segment(12, 3)                 = R_dot.col(2);
  // Omega dot
  X_ddot.segment(15, 3)                 = X_dot_dot.segment(3, 3);
  // alpha dot
  X_ddot(18)                            = alpha_dot_in_EOM;
  // alpha dot dot
  X_ddot(19)                            = X_dot_dot(6);

  // ---------- Debug text ----------
  // ROS_INFO_STREAM_THROTTLE(0.01,"State Vector \n"         << XK.transpose());
  // ROS_INFO_STREAM_THROTTLE(0.01,"Final_Inertia \n"        << Final_Inertia);
  // ROS_INFO_STREAM_THROTTLE(0.01,"coriolis_final \n"       << coriolis_final.transpose());
  // ROS_INFO_STREAM_THROTTLE(0.01,"control_input_column \n" << control_input_column.transpose());
  // ROS_INFO_STREAM_THROTTLE(1,"X_ddot -- final stuff 5 \n" << X_ddot.transpose());

  return X_ddot;

}  // namespace pratik_sim_mrs_net_control_link_only

}
#endif
