/**
 * @brief Dynamic simulation of a multirotor helicopter.
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
// Omega            - 3
// q_ij's           - 3 * n * m
// q_dot_ij's       - 3 * n * m

// Total states     = 18 + 6 * n *m

// N_INTERNAL_STATES = 18 + 6 * NO_OF_LINKS * NO_OF_CAP ;
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

Eigen::Vector3d rho_vector(0.0,0.0,-0.01); // Vertical offset of the net mechanism underneath the drone

constexpr int NO_OF_LINKS         = 2;        // (n)
constexpr int NO_OF_CAP           = 2;        // (m)
constexpr int link_dof            = 3 * NO_OF_LINKS * NO_OF_CAP;
constexpr int q_index_start       = 18;       // MATLAB q_index_start (1-based 19) -> C++ index 18
constexpr int q_dot_index_start   = q_index_start + link_dof;
constexpr int N_INTERNAL_STATES   = 18 + 2 * link_dof;

std::vector<std::vector<Vector3d>> q_each_link_local(NO_OF_CAP, std::vector<Vector3d>(NO_OF_LINKS));
std::vector<std::vector<Vector3d>> q_dot_each_link_local(NO_OF_CAP, std::vector<Vector3d>(NO_OF_LINKS));
std::vector<std::vector<Vector3d>> q_dot_dot_each_link_local(NO_OF_CAP, std::vector<Vector3d>(NO_OF_LINKS));

// Other system parameters
// system parameters
constexpr float m_pt        = 0.01;            // Mass of each point mass
constexpr float l_pt        = 1.0;            // Length of each links
constexpr float mr          = 0.05;            // Mass of net rod
constexpr float lr          = 0.2;            // half length of rod

// Assume point_mass_array and link_lengths_array are MatrixXd of size (m, n)
MatrixXd point_mass_array     = MatrixXd::Constant(NO_OF_CAP, NO_OF_LINKS, m_pt);
MatrixXd link_lengths_array   = MatrixXd::Constant(NO_OF_CAP, NO_OF_LINKS, l_pt);

std::vector<Eigen::Vector3d> rho_i_vectors(NO_OF_CAP);

Eigen::Vector3d e2(0.0,1.0,0.0);
Eigen::Vector3d e3(0.0,0.0,1.0);

// float Ts                  = 0.01;         // sampling time [s]

namespace odeint = boost::numeric::odeint;

namespace pratik_sim_mrs_net_free_ends
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
    // mrs_net_free_ends
    Eigen::Matrix<double, link_dof, 1> q_for_all_link;
    Eigen::Matrix<double, link_dof, 1> q_dot_for_all_link;

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

  // static inline Eigen::Vector3d qb_vec(double alpha);
  static inline Eigen::Matrix3d hatmap(const Eigen::Vector3d &v);
  // static inline Eigen::Vector3d qt_vec(double alpha);
  inline double wrapToPi(double angle);

  static double M_00(const MatrixXd &point_mass_array, double mq, double mr);
  static double M_0ij(int ii, int jj, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array);
  static double M_ijk(int ii, int jj, int kk, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array);
  static Vector3d K1_vector(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, int n);
  static Matrix3d K2_vector(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, int n);
  static Matrix3d J_bar(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors,
                      const Vector3d &e2, const Matrix3d &Jq, int n, int m, double lr);
  static Matrix3d K_ij(int ii, int jj, const Matrix3d &R, const Vector3d &q_each_link_ii_jj, int n,
                     const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array,
                     const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors);
  static Matrix3d K_R_ij_matrix(int ii, int jj, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array,
                              const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors);
  VectorXd EOM_mrs_net_free_ends(
  const VectorXd &XK,             // big state vector
  double mq,
  const Matrix3d &Jq,
  int n,                          // number of link segments per branch
  int m,                          // number of branches
  double m_pt,
  double l_pt,
  const MatrixXd &link_lengths_array,           // m x n
  const Vector3d &rho,
  const std::vector<Vector3d> &rho_i_vectors,   // size m
  double lr,
  double g,
  double m_I,
  double intruder_position,
  MatrixXd point_mass_array,                    // m x n
  int drone_hitting_i_th_cable_attachment_point,
  int drone_hitting_j_th_point_mass,
  const Vector3d &e3,
  double mr,
  const Vector3d &e2,
  const Eigen::Vector3d &M,
  double f
 );


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

  _initial_pos_           = spawn_pos;
  state_.x                = spawn_pos;
  state_.R                = Eigen::AngleAxis(-spawn_heading, Eigen::Vector3d(0, 0, 1));

  for (int i6 = 0; i6 < NO_OF_CAP; ++i6) {
    for (int j6 = 0; j6 < NO_OF_LINKS; ++j6) {
      int start     = (i6) * 3 * NO_OF_LINKS + (3 * j6);

      // assign (make sure indices are inside vector)
      state_.q_for_all_link(start)           =  0.0;
      state_.q_for_all_link(start + 1)       =  0.0;
      state_.q_for_all_link(start + 2)       = -1.0;
    }
  }
  
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

  for (int i6 = 0; i6 < NO_OF_CAP; ++i6) {
    for (int j6 = 0; j6 < NO_OF_LINKS; ++j6) {
      int start = (i6) * 3 * NO_OF_LINKS + (3 * j6);

      // assign (make sure indices are inside vector)
      state_.q_for_all_link(start)           =  0.0;
      state_.q_for_all_link(start + 1)       =  0.0;
      state_.q_for_all_link(start + 2)       = -1.0;

      state_.q_dot_for_all_link(start)      =  0.0;
      state_.q_dot_for_all_link(start + 1)  =  0.0;
      state_.q_dot_for_all_link(start + 2)  =  0.0;

    }
  }

  imu_acceleration_       = Eigen::Vector3d::Zero();

  state_.motor_rpm        = Eigen::VectorXd::Zero(params_.n_motors);
  input_                  = Eigen::VectorXd::Zero(params_.n_motors);

  external_force_.setZero();
  external_moment_.setZero();
}

//}

/* updatedInternalState() //{ */

void MultirotorModel::updateInternalState(void) {

  for (int i = 0; i < 3; ++i) {
    internal_state_.at(0 + i)     = state_.x(i);
    internal_state_.at(3 + i)     = state_.v(i);
    internal_state_.at(6 + i)     = state_.R(i, 0);
    internal_state_.at(9 + i)     = state_.R(i, 1);
    internal_state_.at(12 + i)    = state_.R(i, 2);
    internal_state_.at(15 + i)    = state_.omega(i);
  }

  for (int i = 0; i < link_dof; ++i) {
    internal_state_.at(q_index_start + i)     = state_.q_for_all_link(i);
    internal_state_.at(q_dot_index_start + i) = state_.q_dot_for_all_link(i);
  }

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

  for (int i = 0; i < 3; ++i) {
    state_.x(i)     = internal_state_.at(0 + i);
    state_.v(i)     = internal_state_.at(3 + i);
    state_.R(i, 0)  = internal_state_.at(6 + i);
    state_.R(i, 1)  = internal_state_.at(9 + i);
    state_.R(i, 2)  = internal_state_.at(12 + i);
    state_.omega(i) = internal_state_.at(15 + i);
  }

  for (int i = 0; i < link_dof; ++i) {
    state_.q_for_all_link(i)        = internal_state_.at(q_index_start + i);
    state_.q_dot_for_all_link(i)    = internal_state_.at(q_dot_index_start + i);
  }

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

  for (int i = 0; i < 3; ++i) {
    cur_state.x(i)     = x.at(0 + i);
    cur_state.v(i)     = x.at(3 + i);
    cur_state.R(i, 0)  = x.at(6 + i);
    cur_state.R(i, 1)  = x.at(9 + i);
    cur_state.R(i, 2)  = x.at(12 + i);
    cur_state.omega(i) = x.at(15 + i);
  }

  for (int i = 0; i < link_dof; ++i) {
    cur_state.q_for_all_link(i)           = x.at(q_index_start      + i);
    cur_state.q_dot_for_all_link(i)       = x.at(q_dot_index_start  + i);
  }

  // // Organize q_for_all_link and q_dot_for_all_link in each link seperately
  // std::vector<std::vector<Vector3d>> q_each_link_local(NO_OF_CAP, std::vector<Vector3d>(NO_OF_LINKS));
  // std::vector<std::vector<Vector3d>> q_dot_each_link_local(NO_OF_CAP, std::vector<Vector3d>(NO_OF_LINKS));

  for (int i6 = 0; i6 < NO_OF_CAP; ++i6) {
    for (int j6 = 0; j6 < NO_OF_LINKS; ++j6) {

      int start = (i6) * 3 * NO_OF_LINKS + (3      * j6);

      q_each_link_local[i6][j6]       = cur_state.q_for_all_link.segment(start, 3);

      double norm_q = q_each_link_local[i6][j6].norm();

      if (norm_q > 0) q_each_link_local[i6][j6] /= norm_q;

      q_dot_each_link_local[i6][j6]   = cur_state.q_dot_for_all_link.segment(start, 3);

    }
  }

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(cur_state.R.transpose() * cur_state.R);
  Eigen::Matrix3d             P = llt.matrixL();
  cur_state.R                   = cur_state.R * P.inverse();

  // First time derivatives of states

  Eigen::Vector3d x_dot_local;
  Eigen::Vector3d x_dot_dot_local;
  Eigen::Matrix3d R_dot_local;
  Eigen::Vector3d omega_dot_local;

  Eigen::Matrix3d omega_tensor(Eigen::Matrix3d::Zero());
  omega_tensor(2, 1)            =  cur_state.omega(0);
  omega_tensor(1, 2)            = -cur_state.omega(0);
  omega_tensor(0, 2)            =  cur_state.omega(1);
  omega_tensor(2, 0)            = -cur_state.omega(1);
  omega_tensor(1, 0)            =  cur_state.omega(2);
  omega_tensor(0, 1)            = -cur_state.omega(2);

  Eigen::VectorXd motor_rpm_sq  = state_.motor_rpm.array().square();
  Eigen::Vector4d torque_thrust = params_.allocation_matrix * motor_rpm_sq;
  double          thrust        = torque_thrust(3);

  // double resistance       = params_.air_resistance_coeff * M_PI * (params_.arm_length) * (params_.arm_length) * cur_state.v.norm() * cur_state.v.norm();
  Eigen::Vector3d vnorm   = cur_state.v;
  if (vnorm.norm() != 0) {
    vnorm.normalize();
  }

  // ROS_INFO_STREAM_THROTTLE(1,"Point Mass Array: \n" << point_mass_array);
  // ROS_INFO_STREAM_THROTTLE(1,"Link Length Array:\n" << link_lengths_array);

  double distance_between_each_cable_attachment_point;
  if (NO_OF_CAP > 1){
    distance_between_each_cable_attachment_point = 2.0 * lr / (NO_OF_CAP - 1);
  }else{
    distance_between_each_cable_attachment_point = 0;
  }

  for (int i_rho = 0; i_rho < NO_OF_CAP; ++i_rho) {
    double distance_along_second_axis = lr - i_rho * distance_between_each_cable_attachment_point;
    rho_i_vectors[i_rho] = Vector3d(0.0, distance_along_second_axis, 0.0);
  }

  // ROS_INFO_STREAM("rho_i_vectors:");
  // for (size_t i = 0; i < rho_i_vectors.size(); ++i) {
  //     ROS_INFO_STREAM(" [" << i << "] = " << rho_i_vectors[i].transpose());
  // }

  double m_I                                    = 0.5;        // Mass of intruder drone
  int drone_hitting_i_th_cable_attachment_point = 1;
  int drone_hitting_j_th_point_mass             = 1;
  double intruder_position                      = 0.5;

  // Check feasibility
  int check_feasibility_for_n = drone_hitting_i_th_cable_attachment_point - NO_OF_CAP;
  int check_feasibility_for_m = drone_hitting_j_th_point_mass - NO_OF_LINKS;

  if (check_feasibility_for_n > 0 || check_feasibility_for_m > 0) {
      throw std::runtime_error("Check the values of n and m parameters");
  }

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////

  // ROS_INFO_STREAM("Point Mass Array ->\n" << point_mass_array);

  Eigen::VectorXd XK(N_INTERNAL_STATES,1);

  for (int i = 0; i < N_INTERNAL_STATES; ++i) {
    XK(i) = x.at(i);
  }

  // ROS_INFO_STREAM_THROTTLE(1,"Quad Trans. states Operator(): " << XK.segment(0,6).transpose());

  // ROS_INFO_STREAM_THROTTLE(1,"Quad Rotat. states Operator(): " << XK.segment(6,12).transpose());

  // ROS_INFO_STREAM_THROTTLE(1,"Link states Operator():        " << XK.segment(q_index_start,link_dof).transpose());

  // ROS_INFO_STREAM_THROTTLE(1,"Links dot states Operator():   " << XK.segment(q_dot_index_start,link_dof).transpose());

  // Calling the custom dynamics of the system
  Eigen::VectorXd X_dot_dot_custom(N_INTERNAL_STATES,1);

  X_dot_dot_custom = EOM_mrs_net_free_ends(
  XK,             // current state vector
  params_.mass,
  params_.J,
  NO_OF_LINKS,                  // number of link per cable
  NO_OF_CAP,                    // number of cable attachment points
  m_pt,
  l_pt,
  link_lengths_array,   // m x n
  rho_vector,
  rho_i_vectors, // size m
  lr,
  params_.g,
  m_I,
  intruder_position,
  point_mass_array,      // m x n
  drone_hitting_i_th_cable_attachment_point,
  drone_hitting_j_th_point_mass,
  e3,
  mr,
  e2,
  torque_thrust.topRows(3), thrust);

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

  state_.q_for_all_link        = state.q_for_all_link;
  state_.q_dot_for_all_link    = state.q_dot_for_all_link;

  updateInternalState();
}

//}

/* setStatePos() //{ */

void MultirotorModel::setStatePos(const Eigen::Vector3d& pos, const double heading) {

  _initial_pos_           = pos;
  state_.x                = pos;
  state_.R                = Eigen::AngleAxis(-heading, Eigen::Vector3d(0, 0, 1));

  for (int i6 = 0; i6 < NO_OF_CAP; ++i6) {
    for (int j6 = 0; j6 < NO_OF_LINKS; ++j6) {
      int start     = (i6) * 3 * NO_OF_LINKS + (3 * j6);

      ROS_INFO_THROTTLE(0.1, "Start-> %d", start);

      // assign (make sure indices are inside vector)
      state_.q_for_all_link(start)           =  0.0;
      state_.q_for_all_link(start + 1)       =  0.0;
      state_.q_for_all_link(start + 2)       = -1.0;

    }
  }
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

// // qb_vec and qt_vec as used in MATLAB version
// Eigen::Vector3d MultirotorModel::qb_vec(double alpha) {
//   return Eigen::Vector3d(-std::sin(alpha), 0.0, -std::cos(alpha));
// }

// Eigen::Vector3d MultirotorModel::qt_vec(double alpha) {
//   return Eigen::Vector3d(-std::cos(alpha), 0.0, std::sin(alpha));
// }

// M_00 = mq + mr + sum(point_mass_array, "all")
double MultirotorModel::M_00(const MatrixXd &point_mass_array, double mq, double mr) {
  return mq + mr + point_mass_array.sum();
}

// M_0ij  (note: 1-based MATLAB converted to 0-based C++)

double MultirotorModel::M_0ij(int ii, int jj, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array) {
  // MATLAB: for a = jj:n, sum point_mass_array(ii,a)*link_lengths_array(ii,jj)
  // ii and jj are assumed 1..m in MATLAB; in C++ we expect 0..m-1 / 0..n-1
  double val = 0.0;
  for (int a = jj; a <= n; ++a) {
    val += point_mass_array(ii-1, a-1) * link_lengths_array(ii-1, jj-1);
  }
  return val;
}

// M_ijk
double MultirotorModel::M_ijk(int ii, int jj, int kk, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array) {
  double val = 0.0;
  if (jj != 0 && kk != 0) {
    int start = std::max(jj, kk);
    for (int a = start; a <= n; ++a) {
      val += point_mass_array(ii-1, a-1) * link_lengths_array(ii-1, jj-1) * link_lengths_array(ii-1, kk-1);
    }
  }
  return val;
}

// K1_vector
Vector3d MultirotorModel::K1_vector(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, int n) {
  Vector3d out = mr * rho;
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      out += point_mass_array(i-1, j-1) * (rho + rho_i_vectors[i-1]);
    }
  }
  return out;
}

// K2_vector
Matrix3d MultirotorModel::K2_vector(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors, int m, int n) {
  Matrix3d out = mr * hatmap(rho);
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      out += point_mass_array(i-1, j-1) * hatmap(rho + rho_i_vectors[i-1]);
    }
  }
  return out;
}

// J_bar
Matrix3d MultirotorModel::J_bar(const MatrixXd &point_mass_array, double mr, const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors,
                      const Vector3d &e2, const Matrix3d &Jq, int n, int m, double lr) {
  Matrix3d Jbar = Jq;
  Jbar -= 0.5 * mr * hatmap(rho) * hatmap(rho);
  Jbar -= (1.0/6.0) * mr * lr * lr * hatmap(e2) * hatmap(e2);
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      Matrix3d h = hatmap(rho + rho_i_vectors[i-1]);
      Jbar -= point_mass_array(i-1,j-1) * h * h;
    }
  }
  return Jbar;
}

// K_ij function returning 3x3 (see MATLAB)
Matrix3d MultirotorModel::K_ij(int ii, int jj, const Matrix3d &R, const Vector3d &q_each_link_ii_jj, int n,
                     const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array,
                     const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors) {

  Matrix3d out = Matrix3d::Zero();
  for (int a = 1; a <= n; ++a) {
    out += point_mass_array(ii-1, a-1) * link_lengths_array(ii-1, jj-1) * hatmap(q_each_link_ii_jj) * hatmap(q_each_link_ii_jj) * R * hatmap(rho + rho_i_vectors[ii-1]);
  }
  return out;
}

// K_R_ij_matrix (MATLAB version)
Matrix3d MultirotorModel::K_R_ij_matrix(int ii, int jj, int n, const MatrixXd &point_mass_array, const MatrixXd &link_lengths_array,
                              const Vector3d &rho, const std::vector<Vector3d> &rho_i_vectors) {
  Matrix3d out = Matrix3d::Zero();

  for (int a = jj; a <= n; ++a) {
    // ROS_INFO_STREAM("out \n" << out);
    out += point_mass_array(ii-1, a-1) * link_lengths_array(ii-1, jj-1) * hatmap(rho + rho_i_vectors[ii-1]);
  }

  return out;
}


// ---------- EOM main function ----------
// A single function that mirrors the MATLAB EOM skeleton. Many details are implemented,
// but the final block-matrix assembly is left in a clear TODO so you can verify shapes.
VectorXd MultirotorModel::EOM_mrs_net_free_ends(
  const VectorXd &XK,             // big state vector
  double mq,
  const Matrix3d &Jq,
  int n,                          // number of link segments per branch
  int m,                          // number of branches
  double m_pt,
  double l_pt,
  const MatrixXd &link_lengths_array,   // m x n
  const Vector3d &rho,
  const std::vector<Vector3d> &rho_i_vectors, // size m
  double lr,
  double g,
  double m_I,
  double intruder_position,
  MatrixXd point_mass_array,      // m x n
  int drone_hitting_i_th_cable_attachment_point,
  int drone_hitting_j_th_point_mass,
  const Vector3d &e3,
  double mr,
  const Vector3d &e2,
  const Eigen::Vector3d &M,
  double f)
  {
  // Basic checks
  assert(rho_i_vectors.size()       == static_cast<size_t>(m));
  assert(point_mass_array.rows()    == m && point_mass_array.cols()   == n);
  assert(link_lengths_array.rows()  == m && link_lengths_array.cols() == n);

  // ---------- Recover states from XK ----------
  // MATLAB used 1-based indexing. Convert to 0-based.
  // Expect XK length >= 18 + 6*n*m (as per your final assembly)

  // safety
  if (XK.size() < 18) {
    std::cerr << "XK length too small: " << XK.size() << "\n";
    return VectorXd::Zero(1);
  }

  Vector3d xq;
  xq << XK(0), XK(1), XK(2);

  Vector3d xq_dot;
  xq_dot << XK(3), XK(4), XK(5);

  Matrix3d R;
  // MATLAB assembled R column-wise from indices 7..15 (1-based)
  // MATLAB:
  // R = [ XK(7) XK(10) XK(13);
  //       XK(8) XK(11) XK(14);
  //       XK(9) XK(12) XK(15)];
  R << XK(6),  XK(9),  XK(12),
       XK(7), XK(10),  XK(13),
       XK(8), XK(11),  XK(14);

  // Re-orthonormalize R (polar decomposition)
  Eigen::LLT<Eigen::Matrix3d> llt(R.transpose() * R);
  Eigen::Matrix3d             P = llt.matrixL();
  R = R * P.inverse();

  Vector3d Omega;
  Omega << XK(15), XK(16), XK(17);

  Matrix3d R_dot = R * hatmap(Omega);

  // Now extract q_for_all_link and q_dot_for_all_link
  if (static_cast<int>(XK.size()) < q_index_start + 2*link_dof) {
    std::cerr << "XK size doesn't contain full link state arrays. expected at least " << (q_index_start + 2*link_dof) << " got " << XK.size() << "\n";
    // return zero vector of reasonable size (same as MATLAB X_ddot)
    VectorXd Xddot_fail = VectorXd::Zero(18 + 6 * n * m);
    return Xddot_fail;
  }

  // first define local variables for q_for_all_link_local and q_dot_for_all_link_local variables
  Eigen::Matrix<double, link_dof, 1> q_for_all_link_local;
  Eigen::Matrix<double, link_dof, 1> q_dot_for_all_link_local;

  for (int i = 0; i < link_dof; ++i) {
    q_for_all_link_local(i)           = XK(q_index_start      + i);
    q_dot_for_all_link_local(i)       = XK(q_dot_index_start  + i);
  }

  // ROS_INFO_STREAM_THROTTLE(1,"q_for_all_link_local EOM(): \n" << q_for_all_link_local.transpose());
  // ROS_INFO_STREAM_THROTTLE(1,"q_dot_for_all_link_local EOM(): \n" << q_dot_for_all_link_local.transpose());

  for (int i6 = 0; i6 < m; ++i6) {
    for (int j6 = 0; j6 < n; ++j6) {

      int start = (i6) * 3 * n + (3 * j6);

      q_each_link_local[i6][j6]       = q_for_all_link_local.segment(start, 3);

      double norm_q = q_each_link_local[i6][j6].norm();

      if (norm_q > 0) q_each_link_local[i6][j6] /= norm_q;

      q_dot_each_link_local[i6][j6]   = q_dot_for_all_link_local.segment(start, 3);

      // Eigen::Vector3d v1;
      // v1 = q_each_link_local[i6][j6];
      // Eigen::Vector3d v2;
      // v2 = q_dot_each_link_local[i6][j6];

      // ROS_INFO("q_each_link_debug from EOM [%d][%d] = [%f, %f, %f]",
      //                   i6, j6, v1.x(), v1.y(), v1.z());
      // ROS_INFO("q_dot_each_link_debug from EOM [%d][%d] = [%f, %f, %f]",
      //                   i6, j6, v2.x(), v2.y(), v2.z());
    // ROS_INFO_STREAM_THROTTLE(1,"q_each_link_local EOM(): " << q_each_link_local[i6][j6].transpose());
    // ROS_INFO_STREAM_THROTTLE(1,"q_dot_each_link_local EOM(): " << q_dot_each_link_local[i6][j6].transpose());
    }
  }


  // --------- Intruder mass update (MATLAB logic) ----------
  if (xq(0) >= intruder_position) {
    point_mass_array(drone_hitting_i_th_cable_attachment_point, drone_hitting_j_th_point_mass) = m_I;
  }
  
  Vector3d u_pos                = f * R * e3;
  Vector3d u                    = u_pos;

  // ---------- Build inertia/coriolis blocks ----------
  // sizes:
  int N_states_in_EOM = 6 + 3 * n * m; // total DOF (3 position + 3 attitude + 3*n*m link attitudes)
  // Final_Intertia_matrix_of_System should be N_states_in_EOM x N_states_in_EOM
  MatrixXd Final_Inertia = MatrixXd::Zero(N_states_in_EOM, N_states_in_EOM);

  // Build Inertia_matrix_quad_position_raw (3 x N_states_in_EOM)
  // MATLAB code builds it by concatenating 3x3 blocks: first M_00 * I3, then for each branch/segment M_0ij*I3, and finally -R*hatmap(K1_vector)
  MatrixXd Inertia_matrix_quad_position_raw = MatrixXd::Zero(3, N_states_in_EOM);
  // Column block 0 (first 3 columns)

  // ROS_INFO_STREAM("Inertia_matrix_quad_position_raw \n" << Inertia_matrix_quad_position_raw);

  // // Next blocks: columns correspond to links (3 columns per link?), we need to carefully place them.
  // // TODO: finish assembly exactly as you want. Below is a sketch:

  // column offset after the 3 position columns
  int col_offset = 3;

  for (int i31 = 0; i31 <= m + 1 ; ++i31) {

          if ( i31 == 0 )
          {
                            Inertia_matrix_quad_position_raw.block(0, 0, 3, 3) = Matrix3d::Identity() * M_00(point_mass_array, mq, mr);
          }
          else if ( i31 >= 1 && i31 <= m )
          {
                            for (int j32 = 1; j32 <= n; ++j32) {
                              
                              double val = M_0ij(i31, j32, n, point_mass_array, link_lengths_array);
                              
                              Inertia_matrix_quad_position_raw.block(0, col_offset, 3, 3) = Matrix3d::Identity() * val;
                              col_offset += 3;
                            }
          }
          else if ( i31 == m + 1 )
          {
                              Vector3d K1 = K1_vector(point_mass_array, mr, rho, rho_i_vectors, m, n);
                              Inertia_matrix_quad_position_raw.block(0, col_offset, 3, 3) = -R * hatmap(K1);
          }

  }

  // ROS_INFO_STREAM("Inertia_matrix_quad_position_raw \n" << Inertia_matrix_quad_position_raw);

  Final_Inertia.block(0, 0, 3, N_states_in_EOM) = Inertia_matrix_quad_position_raw;

  // ROS_INFO_STREAM("Inertia_matrix_quad_position_raw \n" << Inertia_matrix_quad_position_raw);

  // ---------- Links attitude inertia assembly ----------
  MatrixXd Final_intertial_matrix_link_attitudes = MatrixXd::Zero(3 * n * m, N_states_in_EOM);
  int row_idx_for_CAP = 0;
  int row_idx = 0;
  int write_col = 3;

  for (int j40 = 1; j40 <= m; ++j40) {

    MatrixXd Inertia_matrix_link_attitude_raw_matrix = MatrixXd::Zero(3*n, N_states_in_EOM);
    row_idx = 0;

    for (int j42 = 1; j42 <= n; ++j42) {

      MatrixXd Inertia_matrix_link_attitude_column_start = MatrixXd::Zero(3, N_states_in_EOM);
      write_col = 3;

                  for (int i41 = 0; i41 <= m + 1; ++i41) {

                    if (i41 == 0)
                    {
                                  // Build column by column into a temporary 3 x N_states_in_EOM block
                                  // i41 == 0 block:
                                  Inertia_matrix_link_attitude_column_start.block(0, 0, 3, 3) = M_0ij(j40, j42, n, point_mass_array, link_lengths_array) * hatmap(q_each_link_local[j40-1][j42-1]) * hatmap(q_each_link_local[j40-1][j42-1]);
                    }
                    else if ( i41 >= 1 && i41 <= m )
                    {

                            if (i41 == j40) {
                              
                                    // append for j43=1..n
                                    for (int j43 = 1; j43 <= n; ++j43) {

                                          if (j43 == j42) 
                                          {
                                            Inertia_matrix_link_attitude_column_start.block(0, write_col, 3, 3) = -M_ijk(j40, j42, j43, n, point_mass_array, link_lengths_array) * Matrix3d::Identity();
                                          } else {
                                            Inertia_matrix_link_attitude_column_start.block(0, write_col, 3, 3) =  M_ijk(j40, j42, j43, n, point_mass_array, link_lengths_array) * hatmap(q_each_link_local[j40-1][j42-1]) * hatmap(q_each_link_local[j40-1][j42-1]);
                                          }
                                          write_col += 3;
                                    }

                            } else 
                            { 
                                    for (int j43 = 1; j43 <= n; ++j43) {

                                          Inertia_matrix_link_attitude_column_start.block(0, write_col, 3, 3) =  Matrix3d::Zero();
                                          write_col += 3;
                                    }
                            }
                    }
                    else if (i41 == m+1)
                    {
                        Inertia_matrix_link_attitude_column_start.block(0, write_col, 3, 3) = -K_ij(j40, j42, R, q_each_link_local[j40-1][j42-1], n, point_mass_array, link_lengths_array, rho, rho_i_vectors);
                    }
                  }

          Inertia_matrix_link_attitude_raw_matrix.block(row_idx,0,3,N_states_in_EOM) = Inertia_matrix_link_attitude_column_start;
          row_idx += 3;
    }

    Final_intertial_matrix_link_attitudes.block(row_idx_for_CAP,0,3*n,N_states_in_EOM) = Inertia_matrix_link_attitude_raw_matrix;
    row_idx_for_CAP += 3*n;
  }

  // place link attitude block into Final_Inertia
  Final_Inertia.block(3, 0, 3 * n * m, N_states_in_EOM) = Final_intertial_matrix_link_attitudes;

  // ROS_INFO_STREAM("Final_Inertia \n" << Final_Inertia);

  // ---------- Quadcopter attitude inertia matrix (3 x N_states_in_EOM) ----------
  MatrixXd quadcopter_attitude_inertia_matrix = MatrixXd::Zero(3, N_states_in_EOM);

  int col_index = 0;

  // ----------------- i51 == 0 -----------------
  quadcopter_attitude_inertia_matrix.block(0, col_index, 3, 3) = K2_vector(point_mass_array, mr, rho, rho_i_vectors, m, n) * R.transpose();
  col_index += 3;

  // ----------------- i51 >= 1 && i51 <= m -----------------
  for (int i51 = 1; i51 <= m; ++i51) {
      for (int j53 = 1; j53 <= n; ++j53) {
          // ROS_INFO_STREAM("before quadcopter_attitude_inertia_matrix \n" << quadcopter_attitude_inertia_matrix);

          quadcopter_attitude_inertia_matrix.block(0, col_index, 3, 3) = K_R_ij_matrix(i51, j53, n, point_mass_array, link_lengths_array, rho, rho_i_vectors) * R.transpose();
            // ROS_INFO_STREAM("after quadcopter_attitude_inertia_matrix \n" << quadcopter_attitude_inertia_matrix);

          col_index += 3; // move to next 3-column block
      }
  }

  // ROS_INFO_STREAM("quadcopter_attitude_inertia_matrix \n" << quadcopter_attitude_inertia_matrix);

  // ----------------- i51 == m+1 -----------------
  Matrix3d Jbar = J_bar(point_mass_array, mr, rho, rho_i_vectors, e2, Jq, n, m, lr);
  quadcopter_attitude_inertia_matrix.block(0, col_index, 3, 3) = Jbar;

  // ----------------- copy to final inertia matrix -----------------
  Final_Inertia.block(3 + 3 * n * m, 0, 3, quadcopter_attitude_inertia_matrix.cols()) = quadcopter_attitude_inertia_matrix;

  // ROS_INFO_STREAM_THROTTLE(1,"Final_Inertia \n" << Final_Inertia);

  // ---------- Coriolis/Centripetal terms ----------
  // Build the Coriolis/centripetal term vectors to match the vertical concatenation in MATLAB
  VectorXd coriolis_final = VectorXd::Zero(N_states_in_EOM);

  // quad position term
  Vector3d coriolis_quad_pos = R * hatmap(Omega) * hatmap(Omega) * K1_vector(point_mass_array, mr, rho, rho_i_vectors, m, n);
  coriolis_quad_pos += M_00(point_mass_array, mq, mr) * g * e3;
  coriolis_final.segment(0, 3) = coriolis_quad_pos;

  // links attitude coriolis: mechanical assembly (sketch)
  VectorXd coriolis_links = VectorXd::Zero(3 * n * m);
  row_idx = 0;
  for (int j60 = 1; j60 <= m; ++j60) {

    double sum_of_m_ij_temp = 0.0;

    for (int j62_a = 1; j62_a <= n; ++j62_a)
    {
      sum_of_m_ij_temp += point_mass_array(j60-1, j62_a-1);
    }

    for (int j62 = 1; j62 <= n; ++j62) {

      Vector3d term = - M_ijk(j60, j62, j62, n, point_mass_array, link_lengths_array) * q_each_link_local[j60-1][j62-1] * q_dot_each_link_local[j60-1][j62-1].norm() * q_dot_each_link_local[j60-1][j62-1].norm();

      term += sum_of_m_ij_temp * link_lengths_array(j60-1, j62-1) * hatmap(q_each_link_local[j60-1][j62-1]) * hatmap(q_each_link_local[j60-1][j62-1]) * R * hatmap(Omega) * hatmap(Omega) * (rho + rho_i_vectors[j60-1]);

      term += M_0ij(j60, j62, n, point_mass_array, link_lengths_array) * g * hatmap(q_each_link_local[j60-1][j62-1]) * hatmap(q_each_link_local[j60-1][j62-1]) * e3;

      coriolis_links.segment(row_idx, 3) = term;

      row_idx += 3;

    }

  }

  coriolis_final.segment(3, 3 * n * m) = coriolis_links;


  // quadcopter attitude coriolis
  Vector3d coriolis_quad_att = hatmap(Omega) * J_bar(point_mass_array, mr, rho, rho_i_vectors, e2, Jq, n, m, lr) * Omega;

  for (int j63 = 1; j63 <= m; ++j63) {

    for (int j64 = 1; j64 <= n; ++j64) {

      Matrix3d KR = K_R_ij_matrix(j63, j64, n, point_mass_array, link_lengths_array, rho, rho_i_vectors);
      coriolis_quad_att -= KR * hatmap(Omega) * R.transpose() * q_dot_each_link_local[j63-1][j64-1];
    }
  
  }

  coriolis_quad_att += K2_vector(point_mass_array, mr, rho, rho_i_vectors, m, n) * g * R.transpose() * e3;
  coriolis_final.segment(3 + 3 * n * m, 3) = coriolis_quad_att;

  // ROS_INFO_STREAM("coriolis_final \n" << coriolis_final);

  // ---------- Control input column ----------
  VectorXd control_input_column = VectorXd::Zero(N_states_in_EOM);
  // first 3 entries = u (force)
  control_input_column.segment(0, 3) = u;
  // next 3*n*m entries = zeros
  // last 3 entries = M (moment)
  control_input_column.segment(N_states_in_EOM - 3, 3) = M;

  // ROS_INFO_STREAM("control_input_column \n" << control_input_column);

  // ---------- Solve for X_dot_dot ----------
  // Try to solve Final_Inertia * X_dot_dot = control_input_column - coriolis_final
  VectorXd rhs = VectorXd::Zero(N_states_in_EOM);
  rhs    = control_input_column - coriolis_final;
  // VectorXd rhs       = VectorXd::Zero(N_states_in_EOM);
  // rhs    = control_input_column;

  // LU solve
  Eigen::PartialPivLU<MatrixXd> lu(Final_Inertia);

  // Compute rank using FullPivLU
  FullPivLU<MatrixXd> lu_decomp(Final_Inertia);
  int rank_of_Final_Inertia_matrix = lu_decomp.rank();
  // ROS_INFO_STREAM_THROTTLE(1,"Rank of Inertia Matrix: " << rank_of_Final_Inertia_matrix);

  VectorXd X_dot_dot;
  // If the matrix is singular, this will produce NaNs or infs; caller should check
  X_dot_dot                           = lu.solve(rhs);

  // ROS_INFO_STREAM_THROTTLE(1,"X_dot_dot \n" << X_dot_dot.transpose());

  VectorXd X_ddot                     = VectorXd::Zero(N_INTERNAL_STATES);
  // VectorXd rhs       = VectorXd::Zero(N_states_in_EOM);

  // xq derivative
  X_ddot.segment(0, 3)                = xq_dot;

  // translational acceleration
  X_ddot.segment(3, 3)                = X_dot_dot.segment(0, 3);

  // ROS_INFO_STREAM("X_ddot -- final stuff  1 \n" << X_ddot);

  // R_dot columns
  X_ddot.segment(6, 3)                = R_dot.col(0);
  X_ddot.segment(9, 3)                = R_dot.col(1);
  X_ddot.segment(12, 3)               = R_dot.col(2);
  // ROS_INFO_STREAM("X_ddot -- final stuff 2 \n" << X_ddot);

  // quad attitude accelerations (last 3 of X_dot_dot)
  X_ddot.segment(15, 3)               = X_dot_dot.segment(N_states_in_EOM - 3, 3);

  // ROS_INFO_STREAM("X_ddot -- final stuff 3 \n" << X_ddot);

  X_ddot.segment(18, link_dof)        = q_dot_for_all_link_local;

  // ROS_INFO_STREAM("X_ddot -- final stuff 4 \n" << X_ddot);

  // remaining accelerations (link attitude accelerations)
  // X_dot_dot(4:end-3) in MATLAB corresponds to indices 4..(N_states_in_EOM-4) (1-based)
  int src_start                           = 3;            // C++ index of 4th element
  int src_len                             = N_states_in_EOM - 6;         // remove first 3 and last 3 -> N_states_in_EOM - 6 elements
  if (src_len > 0) {
    // put this into X_ddot starting at index (18 + 3*n*m)
    int dst_start                         = 18 + link_dof;
    X_ddot.segment(dst_start, link_dof)   = X_dot_dot.segment(src_start, src_len);
  }

  // ROS_INFO_STREAM_THROTTLE(1,"X_ddot -- final stuff 5 \n" << X_ddot.transpose());

  return X_ddot;

}  // namespace pratik_sim_mrs_net_free_ends

}
#endif
