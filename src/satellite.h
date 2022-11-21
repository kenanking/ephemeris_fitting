/*
 * File: satellite.h
 * File Created: Friday, 18th November 2022 10:50:53 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Monday, 21st November 2022 6:22:36 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#ifndef SATELLITE_H
#define SATELLITE_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "constants.h"
#include "gnsstime.h"
#include "sp3.h"
#include "utils.h"

using namespace gnsstime;
using namespace sp3;

const int PARAM_NUM = 15;

// Orbit parameters
struct OrbitParam {
  double toe = 0.0;       // Time of ephemeris
  double sqrtA = 0.0;     // square root of semi-major axis
  double e = 0.0;         // eccentricity
  double omega = 0.0;     // argument of perigee (at time toe)
  double Delta_n = 0.0;   // mean motion correction
  double M0 = 0.0;        // mean anomaly (at time toe)
  double OMEGA = 0.0;     // longtitude of ascending node
  double OMEGA_DOT = 0.0; // rate of change of longitude of the ascending node
  double i0 = 0.0;        // inclination angle (at time toe)
  double IDOT = 0.0;      // rate of change of inclination angle
  double Cuc = 0.0; // amplitude of cosine correction to argument of latitude
  double Cus = 0.0; // amplitude of sine correction to argument of latitude
  double Crc = 0.0; // amplitude of cosine correction to orbital radius
  double Crs = 0.0; // amplitude of sine correction to orbital radius (m)
  double Cic = 0.0; // amplitude of cosine correction to inclination angle
  double Cis = 0.0; // amplitude of sine correction to inclination angle

  OrbitParam() {}
  OrbitParam(const std::vector<double> &param) {
    toe = param[0];
    sqrtA = param[1];
    e = param[2];
    omega = param[3];
    Delta_n = param[4];
    M0 = param[5];
    OMEGA = param[6];
    OMEGA_DOT = param[7];
    i0 = param[8];
    IDOT = param[9];
    Cuc = param[10];
    Cus = param[11];
    Crc = param[12];
    Crs = param[13];
    Cic = param[14];
    Cis = param[15];
  }

  std::vector<double> GetParam() {
    std::vector<double> param(PARAM_NUM + 1, 0.0);
    param[0] = toe;
    param[1] = sqrtA;
    param[2] = e;
    param[3] = omega;
    param[4] = Delta_n;
    param[5] = M0;
    param[6] = OMEGA;
    param[7] = OMEGA_DOT;
    param[8] = i0;
    param[9] = IDOT;
    param[10] = Cuc;
    param[11] = Cus;
    param[12] = Crc;
    param[13] = Crs;
    param[14] = Cic;
    param[15] = Cis;
    return param;
  }
};

struct DeltaOrbitParam {
  /**
   * 初值设置参考下面的论文:
   * 邢志成, 王解先. GPS广播星历参数拟合的雅可比矩阵数值导数计算方法[J].
   * 全球定位系统, 2012, 37(1): 28–31.
   */
  double d_sqrtA = 1e-3;
  double d_e = 1e-7;
  double d_omega = 1e-7;
  double d_Delta_n = 1e-10;
  double d_M0 = 1e-7;
  double d_OMEGA = 1e-7;
  double d_OMEGA_DOT = 1e-10;
  double d_i0 = 1e-6;
  double d_IDOT = 1e-6;
  double d_Cuc = 1e-7;
  double d_Cus = 1e-6;
  double d_Crc = 10;
  double d_Crs = 10;
  double d_Cic = 1e-6;
  double d_Cis = 1e-6;

  DeltaOrbitParam() {}
  DeltaOrbitParam(const std::vector<double> &param) {
    d_sqrtA = param[0];
    d_e = param[1];
    d_omega = param[2];
    d_Delta_n = param[3];
    d_M0 = param[4];
    d_OMEGA = param[5];
    d_OMEGA_DOT = param[6];
    d_i0 = param[7];
    d_IDOT = param[8];
    d_Cuc = param[9];
    d_Cus = param[10];
    d_Crc = param[11];
    d_Crs = param[12];
    d_Cic = param[13];
    d_Cis = param[14];
  }

  std::vector<double> GetParam() {
    std::vector<double> param(PARAM_NUM, 0.0);
    param[0] = d_sqrtA;
    param[1] = d_e;
    param[2] = d_omega;
    param[3] = d_Delta_n;
    param[4] = d_M0;
    param[5] = d_OMEGA;
    param[6] = d_OMEGA_DOT;
    param[7] = d_i0;
    param[8] = d_IDOT;
    param[9] = d_Cuc;
    param[10] = d_Cus;
    param[11] = d_Crc;
    param[12] = d_Crs;
    param[13] = d_Cic;
    param[14] = d_Cis;
    return param;
  }
};

struct FittingOptions {
  int maxIter;             // 最大迭代次数
  double tol;              // 收敛阈值
  int deltaF;              // deltaF的精度：10^deltaF米
  double delta;            // 插值得到toe+Δ时刻的位置
  bool fitKepler;          // 是否使用Kepler初轨作为初值
  bool fixDeltaOrbitParam; // 是否固定数值求导的微小量
  bool verbose;            // 是否打印迭代过程

  FittingOptions() {
    maxIter = 30;     // 默认最大迭代30次
    tol = 1e-3;       // 默认收敛阈值1e-3米
    deltaF = 0;       // 默认使用米级精度
    delta = 1e-6;     // 默认插值1微秒
    fitKepler = true; // 默认使用Kepler初轨作为初值
    fixDeltaOrbitParam = true; // 默认不固定数值求导的微小量（在迭代中更新）
    verbose = true; // 默认打印迭代过程
  }
};

Eigen::Matrix3d GetRotationMatrix(double angle, int axis);

Eigen::Matrix3d GetRotationMatrixDerivative(double angle, int axis);

std::vector<double> CalculateKeplerOrbit(Eigen::Vector3d &x,
                                         Eigen::Vector3d &v);

double Distance(std::vector<double> &x, std::vector<double> &y);

std::vector<double> CalculateSatellitePosition(double t,
                                               const OrbitParam &param);

std::vector<double>
CalculateSatellitePosition(double t, const std::vector<double> &param_array);

std::vector<std::vector<double>> CalculateNumericalDerivative(
    std::vector<double> (*func)(double, const OrbitParam &), double t,
    const std::vector<double> &param, const std::vector<double> &delta_param);

double FitEphemeris(std::vector<PosAndClock> &data, OrbitParam &param,
                    DeltaOrbitParam &delta_param,
                    const FittingOptions &options = FittingOptions());

#endif // SATELLITE_H