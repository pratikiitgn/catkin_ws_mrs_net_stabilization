/**
 * @brief Dynamic simulation of a multirotor helicopter.
 *
 * Acknowledgement:
 * * https://github.com/HKUST-Aerial-Robotics/Fast-Planner
 */

#ifndef MULTIROTOR_MODEL_H
#define MULTIROTOR_MODEL_H

// Number of states
// quad_x           - 3
// quad_v           - 3
// quad_R           - 9
// quad_omega       - 3
// link_alpha       - 1
// link_alpha_dot   - 1
// Total states     = 20

#define N_INTERNAL_STATES 20

#include <boost/array.hpp>
#include "ode/boost/numeric/odeint.hpp"
#include "controllers/references.hpp"

// small regularization to stabilize near-singular Inertia_matrix matrices
static constexpr double REG_EPS = 1e-9;
Eigen::Vector3d rho_vector(0.0,0.0,-0.1);

namespace odeint = boost::numeric::odeint;

namespace pratik_sim_one_link_constraint
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

    // First Rigid Link Only
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

  Eigen::VectorXd compute_Xdd(double mq, double mp, double l, double g, double f,
                     const Eigen::Vector3d &rho, const Eigen::Matrix3d &Jq, const Eigen::Vector3d &M,
                     const Eigen::VectorXd &XK, const Eigen::Vector3d &e3);

  static inline Eigen::Vector3d qb_vec(double alpha);
  static inline Eigen::Matrix3d hatmap(const Eigen::Vector3d &v);
  static inline Eigen::Vector3d qt_vec(double alpha);
  inline double wrapToPi(double angle);

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

// --------------------------------------------------------------
// |                       implementation                       |
// --------------------------------------------------------------

/* constructor MultirotorModel //{ */

MultirotorModel::MultirotorModel(void) {

  initializeState();

  updateInternalState();
}

MultirotorModel::MultirotorModel(const MultirotorModel::ModelParams& params, const Eigen::Vector3d& spawn_pos, const double spawn_heading) {

  params_ = params;

  initializeState();

  _initial_pos_ = spawn_pos;

  state_.x = spawn_pos;
  state_.R = Eigen::AngleAxis(-spawn_heading, Eigen::Vector3d(0, 0, 1));

  // First Rigid Link Only
  state_.alpha      = 0.0;

  updateInternalState();
}

//}

/* initializeState() //{ */

void MultirotorModel::initializeState(void) {

  state_.x      = Eigen::Vector3d::Zero();
  state_.v      = Eigen::Vector3d::Zero();
  state_.v_prev = Eigen::Vector3d::Zero();
  state_.R      = Eigen::Matrix3d::Identity();
  state_.omega  = Eigen::Vector3d::Zero();

  // First Rigid Link Only
  state_.alpha      = 0.0;
  state_.alpha_dot  = 0.0;

  imu_acceleration_   = Eigen::Vector3d::Zero();

  state_.motor_rpm    = Eigen::VectorXd::Zero(params_.n_motors);
  input_              = Eigen::VectorXd::Zero(params_.n_motors);

  external_force_.setZero();
  external_moment_.setZero();
}

//}

/* updatedInternalState() //{ */

void MultirotorModel::updateInternalState(void) {

  for (int i = 0; i < 3; i++) {
    internal_state_.at(0 + i)  = state_.x(i);
    internal_state_.at(3 + i)  = state_.v(i);
    internal_state_.at(6 + i)  = state_.R(i, 0);
    internal_state_.at(9 + i)  = state_.R(i, 1);
    internal_state_.at(12 + i) = state_.R(i, 2);
    internal_state_.at(15 + i) = state_.omega(i);
  }
  // First Rigid Link Only
    internal_state_.at(18) = wrapToPi(state_.alpha);
    internal_state_.at(19) = state_.alpha_dot;
}

//}

/* step() //{ */

void MultirotorModel::step(const double& dt) {

  auto save = internal_state_;

  boost::numeric::odeint::runge_kutta4<InternalState> rk;

  odeint::integrate_n_steps(rk, boost::ref(*this), internal_state_, 0.0, dt, 1);

  for (int i = 0; i < N_INTERNAL_STATES; ++i) {
    if (std::isnan(internal_state_.at(i))) {
      internal_state_ = save;
      break;
    }
  }

  for (int i = 0; i < 3; i++) {
    state_.x(i)     = internal_state_.at(0 + i);
    state_.v(i)     = internal_state_.at(3 + i);
    state_.R(i, 0)  = internal_state_.at(6 + i);
    state_.R(i, 1)  = internal_state_.at(9 + i);
    state_.R(i, 2)  = internal_state_.at(12 + i);
    state_.omega(i) = internal_state_.at(15 + i);
  }

    // First Rigid Link Only
    state_.alpha       = wrapToPi(internal_state_.at(18));
    state_.alpha_dot   = internal_state_.at(19);

  double filter_const = exp((-dt) / (params_.motor_time_constant));

  state_.motor_rpm = filter_const * state_.motor_rpm + (1.0 - filter_const) * input_;

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(state_.R.transpose() * state_.R);

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
  imu_acceleration_ = state_.R.transpose() * (((state_.v - state_.v_prev) / dt) + Eigen::Vector3d(0, 0, params_.g));
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

  State cur_state;

  for (int i = 0; i < 3; i++) {
    cur_state.x(i)     = x.at(0 + i);
    cur_state.v(i)     = x.at(3 + i);
    cur_state.R(i, 0)  = x.at(6 + i);
    cur_state.R(i, 1)  = x.at(9 + i);
    cur_state.R(i, 2)  = x.at(12 + i);
    cur_state.omega(i) = x.at(15 + i);
  }

  // First Rigid Link Only
    cur_state.alpha       = wrapToPi(x.at(18));
    cur_state.alpha_dot   = x.at(19);

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(cur_state.R.transpose() * cur_state.R);
  Eigen::Matrix3d             P = llt.matrixL();
  Eigen::Matrix3d             R = cur_state.R * P.inverse();

  Eigen::Vector3d x_dot;
  Eigen::Vector3d v_dot;
  Eigen::Vector3d omega_dot;
  Eigen::Matrix3d R_dot;
  // First Rigid Link Only
  double alpha_dot_dot;

  Eigen::Matrix3d omega_tensor(Eigen::Matrix3d::Zero());

  omega_tensor(2, 1) =  cur_state.omega(0);
  omega_tensor(1, 2) = -cur_state.omega(0);
  omega_tensor(0, 2) =  cur_state.omega(1);
  omega_tensor(2, 0) = -cur_state.omega(1);
  omega_tensor(1, 0) =  cur_state.omega(2);
  omega_tensor(0, 1) = -cur_state.omega(2);

  Eigen::VectorXd motor_rpm_sq = state_.motor_rpm.array().square();

  Eigen::Vector4d torque_thrust = params_.allocation_matrix * motor_rpm_sq;
  double          thrust        = torque_thrust(3);

  double resistance = params_.air_resistance_coeff * M_PI * (params_.arm_length) * (params_.arm_length) * cur_state.v.norm() * cur_state.v.norm();

  Eigen::Vector3d vnorm = cur_state.v;
  if (vnorm.norm() != 0) {
    vnorm.normalize();
  }

  x_dot = cur_state.v;
  // Eigen::Vector3d u_vector = thrust * R.col(2);

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////

  Eigen::Vector3d e3(0.0,0.0,1.0);

  Eigen::VectorXd XK(20,1);
  
  for (int i = 0; i < 3; i++) {
    XK(0 + i)   = x.at(0 + i);
    XK(3 + i)   = x.at(3 + i);
    XK(6 + i)   = x.at(6 + i);
    XK(9 + i)   = x.at(9 + i);
    XK(12 + i)  = x.at(12 + i);
    XK(15 + i)  = x.at(15 + i);
  }

  // First Rigid Link Only
    XK(18)      = x.at(18);
    XK(19)      = x.at(19);

  // Calling the custom dynamics of the system
  Eigen::Matrix<double, 7, 1> X_dot_dot_custom;
  
  // double thrust_pratik = (params_.mass + params_.mp)*params_.g;
  // Eigen::Vector3d M_pratik (0.0,0.0,0.0);

  // ROS_INFO_THROTTLE("Thrust force: [%2.2f]",thrust_pratik);
  // ROS_INFO_THROTTLE("Control Momemet: [%2.2f,%2.2f,%2.2f]",M_pratik(0),M_pratik(1),M_pratik(2));
  // ROS_INFO_THROTTLE("Initial State: [%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f]",XK(0),XK(1),XK(2),XK(3),XK(4),XK(5),XK(6),XK(7),XK(8),XK(9),XK(10),XK(11),XK(12),XK(13),XK(14),XK(15),XK(16),XK(17),XK(18),XK(19));

  X_dot_dot_custom = compute_Xdd(params_.mass, params_.mp, params_.l, params_.g, thrust,
                     rho_vector, params_.J, torque_thrust.topRows(3), XK, e3);

  // ROS_INFO_THROTTLE(0.5,"Xq dot dot   : [%2.2f,%2.2f,%2.2f]",X_dot_dot_custom(0),X_dot_dot_custom(1),X_dot_dot_custom(2));
  // ROS_INFO_THROTTLE(0.5,"Omega dot dot: [%2.2f,%2.2f,%2.2f]",X_dot_dot_custom(3),X_dot_dot_custom(4),X_dot_dot_custom(5));
  // ROS_INFO_THROTTLE(0.5,"Alpha dot dot: [%2.2f]",X_dot_dot_custom(6));

  for (int jj = 0; jj < 3; jj++) {
    // v_dot = -Eigen::Vector3d(0, 0, params_.g) + thrust * R.col(2) / params_.mass + external_force_ / params_.mass - resistance * vnorm / params_.mass;
    v_dot(jj)       = X_dot_dot_custom(0 + jj);
    // omega_dot = params_.J.inverse() * (torque_thrust.topRows(3) - cur_state.omega.cross(params_.J * cur_state.omega) + external_moment_);
    omega_dot(jj)   = X_dot_dot_custom(3 + jj);
  }

  v_dot = v_dot + external_force_ / params_.mass - resistance * vnorm / params_.mass;

  R_dot = R * omega_tensor;

  // ROS_INFO_THROTTLE("u_vector state -> [%2.3f, %2.3f, %2.3f]",u_vector(0), u_vector(1), u_vector(2));
  alpha_dot_dot    = X_dot_dot_custom(6);

  for (int i = 0; i < 3; i++) {
    dxdt.at(0 + i)  = x_dot(i);
    dxdt.at(3 + i)  = v_dot(i);
    dxdt.at(6 + i)  = R_dot(i, 0);
    dxdt.at(9 + i)  = R_dot(i, 1);
    dxdt.at(12 + i) = R_dot(i, 2);
    dxdt.at(15 + i) = omega_dot(i);
  }
    // First Rigid Link Only
    dxdt.at(18) = cur_state.alpha_dot;
    dxdt.at(19) = alpha_dot_dot;

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

  for (int i = 0; i < params_.n_motors; i++) {

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

  // First Rigid Link Only
  state_.alpha        = wrapToPi(state.alpha);
  state_.alpha_dot    = state.alpha_dot;

  updateInternalState();
}

//}

/* setStatePos() //{ */

void MultirotorModel::setStatePos(const Eigen::Vector3d& pos, const double heading) {

  _initial_pos_ = pos;
  state_.x      = pos;
  state_.R      = Eigen::AngleAxis(-heading, Eigen::Vector3d(0, 0, 1));
  state_.alpha  = 0.0;

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


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// hat map (skew-symmetric) for a 3-vector
Eigen::Matrix3d MultirotorModel::hatmap(const Eigen::Vector3d &v) {
  Eigen::Matrix3d H;
  H <<     0.0, -v.z(),  v.y(),
        v.z(),    0.0, -v.x(),
       -v.y(),  v.x(),    0.0;
  return H;
}

inline double MultirotorModel::wrapToPi(double angle) {
    angle = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0)
        angle += 2.0 * M_PI;
    return angle - M_PI;
}

// qb_vec and qt_vec as used in MATLAB version
Eigen::Vector3d MultirotorModel::qb_vec(double alpha) {
  return Eigen::Vector3d(-std::sin(alpha), 0.0, -std::cos(alpha));
}

Eigen::Vector3d MultirotorModel::qt_vec(double alpha) {
  return Eigen::Vector3d(-std::cos(alpha), 0.0, std::sin(alpha));
}

// compute X_ddot (7x1): [ xq_ddot(3); Omega_dot(3); alpha_ddot(1) ]
// Inputs:
//   mq, mp, l, g, f  - scalars
//   rho (3x1), Jq (3x3), M (3x1), XK (20x1), e3 (3x1)
// XK layout matches the mapping you gave:
//
// XK indices (1-based in matlab) -> 0-based in C++:
//  0..2:   xq
//  3..5:   xq_dot
//  6..14:  R entries (mapping below)
// 15..17:  Omega
// 18:      alpha
// 19:      alpha_dot
Eigen::VectorXd MultirotorModel::compute_Xdd(double mq, double mp, double l, double g, double f,
                     const Eigen::Vector3d &rho, const Eigen::Matrix3d &Jq, const Eigen::Vector3d &M,
                     const Eigen::VectorXd &XK, const Eigen::Vector3d &e3)
{
  // input sanity
  if (XK.size() != 20) {
    throw std::runtime_error("XK must be length 20");
  }

  // extract states
  Eigen::Vector3d xq (XK(0),XK(1),XK(2));
  Eigen::Vector3d xq_dot (XK(3),XK(4),XK(5));

  // reconstruct R from given components as in MATLAB code (watch order)
  // MATLAB used this mapping:
  // R = [ XK(7)   XK(10)  XK(13);
  //       XK(8)   XK(11)  XK(14);
  //       XK(9)   XK(12)  XK(15) ];
  // In 0-based indexing this corresponds to XK(6), XK(9), XK(12), etc.
  Eigen::Matrix3d R;
  R << XK(6),  XK(9),  XK(12),
       XK(7), XK(10),  XK(13),
       XK(8), XK(11),  XK(14);

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(R.transpose() * R);
  Eigen::Matrix3d             P = llt.matrixL();
  R = R * P.inverse();

  Eigen::Vector3d Omega (XK(15),XK(16),XK(17));

  double alpha        = wrapToPi(XK(18));
  double alpha_dot    = XK(19);

  // helper matrices/vectors
  Eigen::Matrix3d hat_rho_lqb  = hatmap(rho + l * qb_vec(alpha)); // hat(rho + l*qb)
  Eigen::Matrix3d hat_qt       = hatmap(qt_vec(alpha));
  Eigen::Matrix3d hat_qb       = hatmap(qb_vec(alpha));
  Eigen::Matrix3d hatOmega     = hatmap(Omega);

  // Build Inertia_matrix sub-blocks (as in MATLAB)
  Eigen::Matrix3d I11 = (mq * Eigen::Matrix3d::Identity()) + (mp * Eigen::Matrix3d::Identity());
  Eigen::Matrix3d I12 = - mp * R * hat_rho_lqb;
  Eigen::Vector3d I13 = mp * l * (R * qt_vec(alpha));             // 3x1

  Eigen::Matrix3d I21 = mp * hat_rho_lqb * R.transpose();
  Eigen::Matrix3d I22 = (Jq - mp * (hat_rho_lqb * hat_rho_lqb));  // Jq - mp * hat(...)^2
  Eigen::Vector3d I23 = - mp * l * (hat_qt * (rho + l * qb_vec(alpha)));

  // bottom blocks
  Eigen::RowVector3d I31 = (mp * l * (qt_vec(alpha).transpose() * R.transpose())); // 1x3
  Eigen::RowVector3d I32 = (mp * l * (rho + l * qb_vec(alpha)).transpose() * hat_qt); // 1x3
  double I33 = mp * l * l;

  // assemble full 7x7 Inertia_matrix matrix
  Eigen::MatrixXd Inertia_matrix = Eigen::MatrixXd::Zero(7,7);
  
  Inertia_matrix.block<3,3>(0,0) = I11;
  Inertia_matrix.block<3,3>(0,3) = I12;
  Inertia_matrix.block<3,1>(0,6) = I13;

  Inertia_matrix.block<3,3>(3,0) = I21;
  Inertia_matrix.block<3,3>(3,3) = I22;
  Inertia_matrix.block<3,1>(3,6) = I23;

  Inertia_matrix.block<1,3>(6,0) = I31;
  Inertia_matrix.block<1,3>(6,3) = I32;
  Inertia_matrix(6,6) = I33;

// ROS_INFO_THROTTLE(0.5,"Inertia_matrix: [%f, %f, %f; %f, %f, %f; %f, %f, %f]",
//          Inertia_matrix(0,0), Inertia_matrix(0,1), Inertia_matrix(0,2),
//          Inertia_matrix(1,0), Inertia_matrix(1,1), Inertia_matrix(1,2),
//          Inertia_matrix(2,0), Inertia_matrix(2,1), Inertia_matrix(2,2));

  // Build coriolis + gravity term (7x1)
  Eigen::Vector3d C1 = mp * R * (hatOmega * hatOmega) * (rho + l * qb_vec(alpha))
               + 2.0 * mp * l * R * hatOmega * qt_vec(alpha) * alpha_dot
               - mp * l * R * qb_vec(alpha) * (alpha_dot * alpha_dot)
               + (mq + mp) * g * e3;

  Eigen::Vector3d C2 = hatOmega * ( Jq - mp * ( hat_rho_lqb * hat_rho_lqb ) ) * Omega
               - mp * l * alpha_dot * hat_rho_lqb * hat_qt * Omega
               + mp * l * (alpha_dot * alpha_dot) * hat_qb * rho
               - mp * l * alpha_dot * hatOmega * hat_qt * ( rho + l * qb_vec(alpha) )
               + mp * g * ( hatmap(rho) + l * hat_qb ) * R.transpose() * e3;

  double C3 = mp * l * ( (Omega.transpose() * hat_rho_lqb * hat_qt * Omega)(0,0) )
            + mp * l * ( (e3.transpose() * (R * qt_vec(alpha)))(0,0) );

  Eigen::VectorXd coriolis_gravity(7);

  coriolis_gravity.segment<3>(0) = C1;
  coriolis_gravity.segment<3>(3) = C2;

  coriolis_gravity(6) = C3;

  // control input column [ f * R * e3; M; 0 ]
  Eigen::VectorXd control(7);
  Eigen::Vector3d fRe3  = f * R * e3;
  control.segment<3>(0) = fRe3;
  control.segment<3>(3) = M;
  control(6) = 0.0;

  // RHS
  Eigen::VectorXd rhs = control - coriolis_gravity;

  // Compute LU decomposition
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(Inertia_matrix);

  // Solve for X_dot_dot
  Eigen::VectorXd X_dot_dot = lu.solve(rhs);

  // Optional: check for numerical issues
  if ((X_dot_dot.array() != X_dot_dot.array()).any()) { // checks for NaN
      std::cerr << "Warning: X_dot_dot contains NaN values. Matrix may be singular!" << std::endl;
  }

  return X_dot_dot;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




}  // namespace pratik_sim_one_link_constraint

#endif
