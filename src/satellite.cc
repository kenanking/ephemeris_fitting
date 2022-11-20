/*
 * File: satellite.cc
 * File Created: Friday, 18th November 2022 10:50:59 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 6:56:44 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#include <cmath>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "satellite.h"

Eigen::Matrix3d GetRotationMatrix(double angle, int axis) {
  Eigen::Matrix3d R;
  switch (axis) {
  case 1: // x轴
    R << 1, 0, 0, 0, cos(angle), sin(angle), 0, -sin(angle), cos(angle);
    break;
  case 2: // y轴
    R << cos(angle), 0, -sin(angle), 0, 1, 0, sin(angle), 0, cos(angle);
    break;
  case 3: // z轴
    R << cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0, 0, 0, 1;
    break;
  default:
    std::cout << "Invalid axis!" << std::endl;
    break;
  }
  return R;
}

Eigen::Matrix3d GetRotationMatrixDerivative(double angle, int axis) {
  Eigen::Matrix3d dR;
  switch (axis) {
  case 1:
    dR << 0, 0, 0, 0, -sin(angle), cos(angle), 0, -cos(angle), -sin(angle);
    break;
  case 2:
    dR << -sin(angle), 0, -cos(angle), 0, 0, 0, cos(angle), 0, -sin(angle);
    break;
  case 3:
    dR << -sin(angle), cos(angle), 0, -cos(angle), -sin(angle), 0, 0, 0, 0;
    break;
  default:
    std::cout << "Invalid axis!" << std::endl;
    break;
  }
  return dR;
}

double LagrangeInterpolate(std::vector<double> &x, std::vector<double> &y,
                           double x0) {
  double mult, sum = 0.0;
  int n = (int)x.size();
  for (int i = 0; i < n; i++) {
    mult = 1.0;
    for (int j = 0; j < n; j++) {
      if (j != i) {
        mult *= (x0 - x[j]) / (x[i] - x[j]);
      }
    }
    sum += mult * y[i];
  }

  return sum;
}

std::vector<double> CalculateKeplerOrbit(Eigen::Vector3d &x,
                                         Eigen::Vector3d &v) {
  Eigen::Vector3d zv = Eigen::Vector3d(0, 0, 1);
  Eigen::Vector3d h = x.cross(v);
  Eigen::Vector3d ev = v.cross(h);
  double r = sqrt(x.dot(x));
  ev = ev / WGS84::GM - x / r;
  double e = sqrt(ev.dot(ev));
  double a = 1 / (2 / r - v.dot(v) / WGS84::GM);

  Eigen::Vector3d vnu = zv.cross(h);
  Eigen::Vector3d eta = h.cross(vnu);

  double w = atan2(ev.dot(eta) / sqrt(h.dot(h)), ev.dot(vnu));
  double Omega = atan2(vnu[1], vnu[0]);

  if (w < 0)
    w += 2 * pi;

  if (Omega < 0)
    Omega += 2 * pi;

  double i0 = atan2(h[0] * h[0] + h[1] * h[1], h[2]);
  if (i0 < 0)
    i0 += 2 * pi;

  double E0 =
      atan2(x.dot(v) / sqrt(WGS84::GM * a), r * v.dot(v) / WGS84::GM - 1);
  double M0 = E0 - e * sin(E0);
  if (M0 < 0)
    M0 += 2 * pi;

  std::vector<double> kepler{a, e, w, Omega, i0, M0};

  return kepler;
}

std::vector<double> CalculateSatellitePosition(double t,
                                               const OrbitParam &param) {
  double tk = t - param.toe;
  if (abs(tk) > 7200)
    std::cout << "Warning: time difference is too large!" << std::endl;

  double a = param.sqrtA * param.sqrtA;
  double n0 = sqrt(WGS84::GM / (a * a * a));
  double n = n0 + param.Delta_n;
  double Mk = param.M0 + n * tk;

  double Ek = Mk;
  double Ek_old = 0.0;
  for (int i = 0; i < 10 && abs(Ek - Ek_old) > 1e-12; i++) {
    Ek_old = Ek;
    Ek = Mk + param.e * sin(Ek_old);
  }

  double vkc = (cos(Ek) - param.e) / (1 - param.e * cos(Ek));
  double vks =
      (sqrt(1 - param.e * param.e) * sin(Ek)) / (1 - param.e * cos(Ek));
  double vk = atan2(vks, vkc);

  double uk = vk + param.omega;

  double delta_uk = param.Cuc * cos(2 * uk) + param.Cus * sin(2 * uk);
  double delta_rk = param.Crc * cos(2 * uk) + param.Crs * sin(2 * uk);
  double delta_ik = param.Cic * cos(2 * uk) + param.Cis * sin(2 * uk);

  double u = uk + delta_uk;
  double r = a * (1 - param.e * cos(Ek)) + delta_rk;
  double i = param.i0 + param.IDOT * tk + delta_ik;

  double x = r * cos(u);
  double y = r * sin(u);

  double lamb = param.OMEGA + (param.OMEGA_DOT - WGS84::OMEGA_E) * tk -
                WGS84::OMEGA_E * param.toe;

  double X = r * (cos(u) * cos(lamb) - sin(u) * sin(lamb) * cos(i));
  double Y = r * (cos(u) * sin(lamb) + sin(u) * cos(lamb) * cos(i));
  double Z = r * (sin(u) * sin(i));

  return {X, Y, Z};
}

std::vector<double>
CalculateSatellitePosition(double t, const std::vector<double> &param_array) {
  OrbitParam param = OrbitParam(param_array);
  return CalculateSatellitePosition(t, param);
}

std::vector<std::vector<double>> CalculateNumericalDerivative(
    std::vector<double> (*func)(double, const OrbitParam &), double t,
    const std::vector<double> &param, const std::vector<double> &delta_param) {
  std::vector<double> coord_0 = func(t, param);
  double x0 = coord_0[0], y0 = coord_0[1], z0 = coord_0[2];

  std::vector<std::vector<double>> df_dparam;
  for (int i = 0; i < PARAM_NUM; i++) {
    std::vector<double> param_tmp = param;
    param_tmp[i + 1] += delta_param[i];
    std::vector<double> coord_1 = func(t, param_tmp);
    double x1 = coord_1[0], y1 = coord_1[1], z1 = coord_1[2];
    df_dparam.push_back({(x1 - x0) / delta_param[i], (y1 - y0) / delta_param[i],
                         (z1 - z0) / delta_param[i]});
  }

  return df_dparam;
}

double Distance(std::vector<double> &x, std::vector<double> &y) {
  double sum = 0.0;
  for (int i = 0; i < x.size(); i++) {
    sum += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return sqrt(sum);
}

void FitEphemeris(std::vector<PosAndClock> &data, OrbitParam &param,
                  DeltaOrbitParam &delta_param, const FittingOptions &options) {
  int n = (int)data.size(); // 拟合数据大小
  double toe = param.toe;   // 星历参考时刻

  // 将参数结构体变为数组，便于使用循环操作
  std::vector<double> param_array = param.GetParam();
  std::vector<double> delta_param_array = delta_param.GetParam();

#pragma region Kepler_parameter
  if (options.fitKepler) { // 拟合Kepler轨道参数
    std::vector<double> ts(n);
    std::vector<double> xs(n);
    std::vector<double> ys(n);
    std::vector<double> zs(n);
    for (int i = 0; i < n; i++) {
      ts[i] = data[i].toe;
      xs[i] = data[i].x * 1000.0; // 单位：m
      ys[i] = data[i].y * 1000.0; // 单位：m
      zs[i] = data[i].z * 1000.0; // 单位：m
    }

    // 插值得到toe时刻的位置
    double x0 = LagrangeInterpolate(ts, xs, toe);
    double y0 = LagrangeInterpolate(ts, ys, toe);
    double z0 = LagrangeInterpolate(ts, zs, toe);

    // 插值得到toe+Δ时刻的位置
    double x1 = LagrangeInterpolate(ts, xs, toe + options.delta);
    double y1 = LagrangeInterpolate(ts, ys, toe + options.delta);
    double z1 = LagrangeInterpolate(ts, zs, toe + options.delta);

    // 得到地心地固系下的位置和速度
    Eigen::Vector3d xp(x0, y0, z0);
    Eigen::Vector3d xv((x1 - x0) / options.delta, (y1 - y0) / options.delta,
                       (z1 - z0) / options.delta);

    // 计算地球自转角
    double angle = WGS84::OMEGA_E * toe;
    Eigen::Matrix3d R, dR;
    R = GetRotationMatrix(-angle, 3);
    dR = GetRotationMatrixDerivative(-angle, 3);

    // 将位置和速度转换到J2000惯性系
    Eigen::Vector3d pos_tmp = xp, vel_tmp = xv;
    xp = R * pos_tmp;
    xv = R * vel_tmp - WGS84::OMEGA_E * dR * pos_tmp;

    // 通过toe时刻的位置和速度，计算开普勒初轨
    std::vector<double> kepler_param = CalculateKeplerOrbit(xp, xv);

    param_array[1] = sqrt(kepler_param[0]); // sqrtA
    param_array[2] = kepler_param[1];       // e
    param_array[3] = kepler_param[2];       // omega
    param_array[6] = kepler_param[3];       // OMEGA
    param_array[8] = kepler_param[4];       // i0
    param_array[5] = kepler_param[5];       // M0

    if (options.verbose) {
      std::cout << std::endl;
      std::cout << "- Using Kepler orbit to initialize the orbit parameters."
                << std::endl;
      std::cout << "- Initial values of Kepler parameters:" << std::endl;
      std::cout << " - sqrtA: " << param_array[1] << std::endl;
      std::cout << " - e: " << param_array[2] << std::endl;
      std::cout << " - omega: " << param_array[3] << std::endl;
      std::cout << " - OMEGA: " << param_array[6] << std::endl;
      std::cout << " - i0: " << param_array[8] << std::endl;
      std::cout << " - M0: " << param_array[5] << std::endl;
      std::cout << std::endl;
    }
  } else {
    if (options.verbose) {
      std::cout << std::endl;
      std::cout << "- Using the user-defined orbit parameters." << std::endl;
      std::cout << "- Initial values of orbit parameters:" << std::endl;
      std::cout << " - sqrtA: " << param_array[1] << std::endl;
      std::cout << " - e: " << param_array[2] << std::endl;
      std::cout << " - omega: " << param_array[3] << std::endl;
      std::cout << " - OMEGA: " << param_array[6] << std::endl;
      std::cout << " - i0: " << param_array[8] << std::endl;
      std::cout << " - M0: " << param_array[5] << std::endl;
      std::cout << std::endl;
    }
  }
#pragma endregion

#pragma region Small_quantities
  if (!options.fixDeltaOrbitParam) {
    std::vector<double> coords_0, coords_1;
    coords_0 = CalculateSatellitePosition(toe + 100, param_array);

    for (int i = 0; i < PARAM_NUM; i++) {
      std::vector<double> param_array_tmp = param_array;
      param_array_tmp[i + 1] += delta_param_array[i];
      coords_1 = CalculateSatellitePosition(toe + 100, param_array_tmp);
      double dist = Distance(coords_0, coords_1);
      delta_param_array[i] *= pow(10, (options.deltaF - (int)log10(dist)));
    }

    if (options.verbose) {
      delta_param = DeltaOrbitParam(delta_param_array);
      std::cout << "- Ajusted delta orbit parameters:" << std::endl;
      std::cout << "- Initial Delta parameters:" << std::endl;
      std::cout << std::setprecision(0) << std::scientific;
      std::cout << " - d_sqrtA: " << delta_param.d_sqrtA << std::endl;
      std::cout << " - d_e: " << delta_param.d_e << std::endl;
      std::cout << " - d_omega: " << delta_param.d_omega << std::endl;
      std::cout << " - d_Delta_n: " << delta_param.d_Delta_n << std::endl;
      std::cout << " - d_M0: " << delta_param.d_M0 << std::endl;
      std::cout << " - d_OMEGA: " << delta_param.d_OMEGA << std::endl;
      std::cout << " - d_OMEGA_DOT: " << delta_param.d_OMEGA_DOT << std::endl;
      std::cout << " - d_i0: " << delta_param.d_i0 << std::endl;
      std::cout << " - d_IDOT: " << delta_param.d_IDOT << std::endl;
      std::cout << " - d_Cuc: " << delta_param.d_Cuc << std::endl;
      std::cout << " - d_Cus: " << delta_param.d_Cus << std::endl;
      std::cout << " - d_Crc: " << delta_param.d_Crc << std::endl;
      std::cout << " - d_Crs: " << delta_param.d_Crs << std::endl;
      std::cout << " - d_Cic: " << delta_param.d_Cic << std::endl;
      std::cout << " - d_Cis: " << delta_param.d_Cis << std::endl;
      std::cout << std::endl;
    }
  } else {
    delta_param = DeltaOrbitParam(delta_param_array);
    if (options.verbose) {
      std::cout << "- Using the user-defined delta orbit parameters."
                << std::endl;
      std::cout << "- Initial Delta parameters:" << std::endl;
      std::cout << std::setprecision(0) << std::scientific;
      std::cout << " - d_sqrtA: " << delta_param.d_sqrtA << std::endl;
      std::cout << " - d_e: " << delta_param.d_e << std::endl;
      std::cout << " - d_omega: " << delta_param.d_omega << std::endl;
      std::cout << " - d_Delta_n: " << delta_param.d_Delta_n << std::endl;
      std::cout << " - d_M0: " << delta_param.d_M0 << std::endl;
      std::cout << " - d_OMEGA: " << delta_param.d_OMEGA << std::endl;
      std::cout << " - d_OMEGA_DOT: " << delta_param.d_OMEGA_DOT << std::endl;
      std::cout << " - d_i0: " << delta_param.d_i0 << std::endl;
      std::cout << " - d_IDOT: " << delta_param.d_IDOT << std::endl;
      std::cout << " - d_Cuc: " << delta_param.d_Cuc << std::endl;
      std::cout << " - d_Cus: " << delta_param.d_Cus << std::endl;
      std::cout << " - d_Crc: " << delta_param.d_Crc << std::endl;
      std::cout << " - d_Crs: " << delta_param.d_Crs << std::endl;
      std::cout << " - d_Cic: " << delta_param.d_Cic << std::endl;
      std::cout << " - d_Cis: " << delta_param.d_Cis << std::endl;
      std::cout << std::endl;
    }
  }
#pragma endregion

#pragma region Optimize
  double sigma_0_old = INFINITY, sigma_0 = 0.0;
  std::cout << "- Optimize:" << std::endl;
  for (int iter = 0; iter < options.maxIter; iter++) {
    Eigen::MatrixXd B(3 * n, 15); // 系数矩阵 (3n, 15)
    Eigen::VectorXd L(3 * n);     // 常数向量 (3n, 1)
    Eigen::VectorXd X(15);        // 改正数向量 (15, 1)
    std::vector<std::vector<double>> df_dparam;

    // 计算系数矩阵和常数向量
    for (int i = 0; i < n; i++) {
      // 计算数值导数
      PosAndClock d = data[i];
      df_dparam = CalculateNumericalDerivative(
          CalculateSatellitePosition, d.toe, param_array, delta_param_array);

      // 计算系数矩阵
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 15; k++) {
          B(i * 3 + j, k) = df_dparam[k][j];
        }
      }

      // 计算常数向量
      std::vector<double> coords =
          CalculateSatellitePosition(d.toe, param_array);
      L(i * 3) = d.x * 1e3 - coords[0];
      L(i * 3 + 1) = d.y * 1e3 - coords[1];
      L(i * 3 + 2) = d.z * 1e3 - coords[2];
    }

    // 计算改正数向量
    X = (B.transpose() * B).inverse() * B.transpose() * L;

    // 更新参数
    for (int i = 0; i < PARAM_NUM; i++) {
      param_array[i + 1] += X(i);
    }

    // 计算残差
    Eigen::VectorXd V = B * X - L;
    double residual = V.dot(V);
    sigma_0 = sqrt(residual / (3 * n - 15));

    if (options.verbose) {
      std::cout << std::setprecision(6) << std::fixed;
      std::cout << " - Iter: " << std::setw(4) << iter
                << ", Residual: " << std::setw(16) << residual
                << ", sigma_0: " << std::setw(16) << sigma_0 << "m"
                << std::endl;
    }

    // 判断是否满足停止迭代阈值
    if (abs(sigma_0 - sigma_0_old) < options.tol) {
      std::cout << std::endl;
      std::cout << " - Converged!" << std::endl;
      break;
    } else {
      sigma_0_old = sigma_0;
    }

    if (!options.fixDeltaOrbitParam) {
      // 更新delta_param_array
      std::vector<double> coords_0, coords_1;
      coords_0 = CalculateSatellitePosition(toe + 100, param_array);

      for (int i = 0; i < PARAM_NUM; i++) {
        std::vector<double> param_array_tmp = param_array;
        param_array_tmp[i + 1] += delta_param_array[i];
        coords_1 = CalculateSatellitePosition(toe + 100, param_array_tmp);
        double dist = Distance(coords_0, coords_1);
        delta_param_array[i] *= pow(10, (options.deltaF - (int)log10(dist)));
      }
    }
  }
  std::cout << std::endl;

#pragma endregion

  param = OrbitParam(param_array);
  delta_param = DeltaOrbitParam(delta_param_array);

  if (options.verbose) {
    std::cout << "- Delta parameters after optimization:" << std::endl;
    std::cout << std::setprecision(0) << std::scientific;
    std::cout << " - d_sqrtA: " << delta_param.d_sqrtA << std::endl;
    std::cout << " - d_e: " << delta_param.d_e << std::endl;
    std::cout << " - d_omega: " << delta_param.d_omega << std::endl;
    std::cout << " - d_Delta_n: " << delta_param.d_Delta_n << std::endl;
    std::cout << " - d_M0: " << delta_param.d_M0 << std::endl;
    std::cout << " - d_OMEGA: " << delta_param.d_OMEGA << std::endl;
    std::cout << " - d_OMEGA_DOT: " << delta_param.d_OMEGA_DOT << std::endl;
    std::cout << " - d_i0: " << delta_param.d_i0 << std::endl;
    std::cout << " - d_IDOT: " << delta_param.d_IDOT << std::endl;
    std::cout << " - d_Cuc: " << delta_param.d_Cuc << std::endl;
    std::cout << " - d_Cus: " << delta_param.d_Cus << std::endl;
    std::cout << " - d_Crc: " << delta_param.d_Crc << std::endl;
    std::cout << " - d_Crs: " << delta_param.d_Crs << std::endl;
    std::cout << " - d_Cic: " << delta_param.d_Cic << std::endl;
    std::cout << " - d_Cis: " << delta_param.d_Cis << std::endl;
    std::cout << std::endl;
  }

  return;
}