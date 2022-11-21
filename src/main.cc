/*
 * File: main.cc
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Monday, 21st November 2022 6:57:54 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <json.hpp>

#include "gnsstime.h"
#include "satellite.h"
#include "sp3.h"

using namespace std;
using namespace Eigen;

using namespace sp3;
using namespace gnsstime;

int main(int, char **) {

  // 读取IGS精密星历文件
  string sp3_filepath = "E:\\Sync\\ephemeris_fitting\\data\\igs22296.sp3";
  SP3File sp3File;
  sp3File.readFile(sp3_filepath);

  // 设置拟合星历参考时间
  GNSSTime fit_time = GNSSTime(2022, 10, 1, 12, 0, 0, TIME_SYSTEM::GPST);
  int week;
  double sow;
  fit_time.ToGPSTime(&week, &sow);

  // 设置拟合的卫星PRN号
  string prn = "G01";

  // 找到距离toe时刻最近的K个数据用于拟合
  // 如果存在两个或多个相同的最小值，num_data可能大于给定值
  int num_data = 8;
  vector<PosAndClock> sp3Data =
      sp3File.GetKNearestNeighbors(prn, fit_time, num_data);
  num_data = (int)sp3Data.size();

  // 获取sp3Data中的起始时间和结束时间
  GNSSTime start_time = fit_time - (sow - sp3Data[0].toe);
  GNSSTime end_time = fit_time + (sp3Data[num_data - 1].toe - sow);

  cout << "========================== Initial Data =========================="
       << endl;
  cout << endl;
  cout << "- Fitting time: " << fit_time.ToString(TIME_SYSTEM::GPST) << endl;
  cout << "- GPS week: " << week << ", seconds of week: " << sow << endl;
  cout << "- Start time: " << start_time.ToString(TIME_SYSTEM::GPST) << endl;
  cout << "- End time: " << end_time.ToString(TIME_SYSTEM::GPST) << endl;
  cout << "- Interval: " << sp3Data[num_data - 1].toe - sp3Data[0].toe << "s"
       << endl;
  cout << "- PRN: " << prn << endl;
  cout << "- Number of data: " << sp3Data.size() << endl;
  cout << "- Data List: " << endl;
  for (auto &data : sp3Data) {
    cout << " - Toe: " << setw(14) << setprecision(6) << fixed << data.toe
         << ", X: " << setw(14) << setprecision(6) << fixed << data.x
         << ", Y: " << setw(14) << setprecision(6) << fixed << data.y
         << ", Z: " << setw(14) << setprecision(6) << fixed << data.z
         << ", Clock: " << setw(14) << setprecision(6) << fixed << data.clock
         << endl;
  }
  cout << endl;
  cout << "================================================================="
       << endl;
  cout << endl;

  // 初始化参数，toe必须是已知的
  OrbitParam param;
  param.toe = sow;

  // 如果不使用开普勒初轨，需要设置初始的轨道参数
  //   param.sqrtA = 5153.661;
  //   param.e = 0.01201;
  //   param.omega = 0.94435;
  //   param.OMEGA = 0.36391;
  //   param.i0 = 0.98890;
  //   param.M0 = 2.82278;

  // 拟合相关设置
  FittingOptions options;
  options.maxIter = 150; // 最大迭代次数
  options.tol = 1e-5; // 迭代收敛阈值，前后两次中误差的变化小于该值时
  options.deltaF = 0;       // ΔF采用米级精度
  options.delta = 1e-6;     // 计算toe时刻速度使用的微小量
  options.fitKepler = true; // 是否使用Kepler初轨作为初始值
  options.fixDeltaOrbitParam =
      false; // 是否固定数值求导的微小量（在迭代中更新）
  options.verbose = true; // 输出迭代过程

  // 设置数值导数计算时的微小量初值
  DeltaOrbitParam delta_param = DeltaOrbitParam();

  // 拟合星历
  cout << "========================== Fitting Orbit ========================="
       << endl;

  double sigma_0 = FitEphemeris(sp3Data, param, delta_param, options);

  cout << "================================================================="
       << endl;
  cout << endl;

  // 输出拟合结果
  cout << "========================== Fitting Result ========================"
       << endl;
  if (sigma_0 == -1) {
    cout << "- Fitting failed!" << endl;
  } else {
    cout << "- Fitting success!" << endl;
    cout << " - Sigma_0: " << sigma_0 << endl;
    cout << setprecision(11) << scientific << endl;
    cout << "- Orbit Parameters: " << endl;
    cout << " - Toe: " << param.toe << endl;
    cout << " - sqrtA: " << param.sqrtA << endl;
    cout << " - e: " << param.e << endl;
    cout << " - omega: " << param.omega << endl;
    cout << " - Delta_n: " << param.Delta_n << endl;
    cout << " - M0: " << param.M0 << endl;
    cout << " - OMEGA: " << param.OMEGA << endl;
    cout << " - OMEGA_DOT: " << param.OMEGA_DOT << endl;
    cout << " - i0: " << param.i0 << endl;
    cout << " - IDOT: " << param.IDOT << endl;
    cout << " - Cuc: " << param.Cuc << endl;
    cout << " - Cus: " << param.Cus << endl;
    cout << " - Crc: " << param.Crc << endl;
    cout << " - Crs: " << param.Crs << endl;
    cout << " - Cic: " << param.Cic << endl;
    cout << " - Cis: " << param.Cis << endl;
  }
  cout << endl;
  cout << "================================================================="
       << endl;
  cout << endl;

  // 验证拟合结果的精度
  cout << "========================== Verify Result ========================="
       << endl;
  cout << endl;

  // 验证内插精度
  // 对精密星历数据按每30秒取一个点进行插值，然后用拟合参数计算
  vector<double> toe_list;
  vector<double> xs_interp, ys_interp, zs_interp;
  vector<double> xs_interp_fit, ys_interp_fit, zs_interp_fit;
  GNSSTime t_interp = start_time;
  double interval = 30;
  while (t_interp < end_time) {
    // 计算内插结果
    PosAndClock d = sp3File.GetInterpolatedData(prn, t_interp, 10);
    toe_list.push_back(d.toe);
    xs_interp.push_back(d.x * 1000.0); // 单位转换为米
    ys_interp.push_back(d.y * 1000.0); // 单位转换为米
    zs_interp.push_back(d.z * 1000.0); // 单位转换为米

    // 计算拟合结果
    vector<double> coords = CalculateSatellitePosition(d.toe, param);
    xs_interp_fit.push_back(coords[0]);
    ys_interp_fit.push_back(coords[1]);
    zs_interp_fit.push_back(coords[2]);

    // 下一个插值点
    t_interp = t_interp + interval;
  }

  // 计算内插误差
  vector<double> xs_err, ys_err, zs_err;
  double x_err, y_err, z_err;
  double sum_x_err = 0, sum_y_err = 0, sum_z_err = 0;
  for (size_t i = 0; i < xs_interp.size(); i++) {
    x_err = xs_interp[i] - xs_interp_fit[i];
    y_err = ys_interp[i] - ys_interp_fit[i];
    z_err = zs_interp[i] - zs_interp_fit[i];
    xs_err.push_back(x_err);
    ys_err.push_back(y_err);
    zs_err.push_back(z_err);
    sum_x_err += x_err * x_err;
    sum_y_err += y_err * y_err;
    sum_z_err += z_err * z_err;
  }
  double rms_x_err = sqrt(sum_x_err / xs_interp.size());
  double rms_y_err = sqrt(sum_y_err / ys_interp.size());
  double rms_z_err = sqrt(sum_z_err / zs_interp.size());

  cout << "- Interpolation Error: " << endl;
  for (size_t i = 0; i < xs_interp.size(); i++) {
    cout << " - " << toe_list[i] << " " << xs_err[i] << " " << ys_err[i] << " "
         << zs_err[i] << endl;
  }
  cout << endl;
  cout << " - Interpolation RMS Error: " << rms_x_err << " " << rms_y_err << " "
       << rms_z_err << endl;
  cout << endl;

  // 验证外推精度
  // 对拟合参数左侧外推1小时
  vector<double> toe_list_left;
  vector<double> xs_left, ys_left, zs_left;
  vector<double> xs_left_fit, ys_left_fit, zs_left_fit;
  GNSSTime t_left = start_time - 3600;
  while (t_left < start_time) {
    // 计算外推结果
    PosAndClock d = sp3File.GetInterpolatedData(prn, t_left, 10);
    toe_list_left.push_back(d.toe);
    xs_left.push_back(d.x * 1000.0); // 单位转换为米
    ys_left.push_back(d.y * 1000.0); // 单位转换为米
    zs_left.push_back(d.z * 1000.0); // 单位转换为米

    // 计算拟合结果
    vector<double> coords = CalculateSatellitePosition(d.toe, param);
    xs_left_fit.push_back(coords[0]);
    ys_left_fit.push_back(coords[1]);
    zs_left_fit.push_back(coords[2]);

    // 下一个外推点
    t_left = t_left + interval;
  }

  // 计算外推误差
  vector<double> xs_err_left, ys_err_left, zs_err_left;
  double x_err_left, y_err_left, z_err_left;
  double sum_x_err_left = 0, sum_y_err_left = 0, sum_z_err_left = 0;
  for (size_t i = 0; i < xs_left.size(); i++) {
    x_err_left = xs_left[i] - xs_left_fit[i];
    y_err_left = ys_left[i] - ys_left_fit[i];
    z_err_left = zs_left[i] - zs_left_fit[i];
    xs_err_left.push_back(x_err_left);
    ys_err_left.push_back(y_err_left);
    zs_err_left.push_back(z_err_left);
    sum_x_err_left += x_err_left * x_err_left;
    sum_y_err_left += y_err_left * y_err_left;
    sum_z_err_left += z_err_left * z_err_left;
  }

  double rms_x_err_left = sqrt(sum_x_err_left / xs_left.size());
  double rms_y_err_left = sqrt(sum_y_err_left / ys_left.size());
  double rms_z_err_left = sqrt(sum_z_err_left / zs_left.size());

  cout << "- Left Extrapolation Error: " << endl;
  for (size_t i = 0; i < xs_left.size(); i++) {
    cout << " - " << toe_list_left[i] << " " << xs_err_left[i] << " "
         << ys_err_left[i] << " " << zs_err_left[i] << endl;
  }
  cout << endl;
  cout << " - Left Extrapolation RMS Error: " << rms_x_err_left << " "
       << rms_y_err_left << " " << rms_z_err_left << endl;
  cout << endl;

  // 对拟合参数右侧外推1小时
  vector<double> toe_list_right;
  vector<double> xs_right, ys_right, zs_right;
  vector<double> xs_right_fit, ys_right_fit, zs_right_fit;
  GNSSTime t_right = end_time + 3600;
  while (t_right > end_time) {
    // 计算外推结果
    PosAndClock d = sp3File.GetInterpolatedData(prn, t_right, 10);
    toe_list_right.push_back(d.toe);
    xs_right.push_back(d.x * 1000.0); // 单位转换为米
    ys_right.push_back(d.y * 1000.0); // 单位转换为米
    zs_right.push_back(d.z * 1000.0); // 单位转换为米

    // 计算拟合结果
    vector<double> coords = CalculateSatellitePosition(d.toe, param);
    xs_right_fit.push_back(coords[0]);
    ys_right_fit.push_back(coords[1]);
    zs_right_fit.push_back(coords[2]);

    // 下一个外推点
    t_right = t_right - interval;
  }

  // 计算外推误差
  vector<double> xs_err_right, ys_err_right, zs_err_right;
  double x_err_right, y_err_right, z_err_right;
  double sum_x_err_right = 0, sum_y_err_right = 0, sum_z_err_right = 0;
  for (size_t i = 0; i < xs_right.size(); i++) {
    x_err_right = xs_right[i] - xs_right_fit[i];
    y_err_right = ys_right[i] - ys_right_fit[i];
    z_err_right = zs_right[i] - zs_right_fit[i];
    xs_err_right.push_back(x_err_right);
    ys_err_right.push_back(y_err_right);
    zs_err_right.push_back(z_err_right);
    sum_x_err_right += x_err_right * x_err_right;
    sum_y_err_right += y_err_right * y_err_right;
    sum_z_err_right += z_err_right * z_err_right;
  }

  double rms_x_err_right = sqrt(sum_x_err_right / xs_right.size());
  double rms_y_err_right = sqrt(sum_y_err_right / ys_right.size());
  double rms_z_err_right = sqrt(sum_z_err_right / zs_right.size());

  cout << "- Right Extrapolation Error: " << endl;
  for (size_t i = 0; i < xs_right.size(); i++) {
    cout << " - " << toe_list_right[i] << " " << xs_err_right[i] << " "
         << ys_err_right[i] << " " << zs_err_right[i] << endl;
  }

  cout << endl;

  cout << " - Right Extrapolation RMS Error: " << rms_x_err_right << " "
       << rms_y_err_right << " " << rms_z_err_right << endl;

  cout << endl;
  cout << "================================================================="
       << endl;
  cout << endl;

  // 结果输出到JSON文件
  using json = nlohmann::json;

  json j;
  j["prn"] = prn;
  j["start_time"] = start_time.ToString();
  j["end_time"] = end_time.ToString();
  j["interval"] = interval;
  j["param"] = param.GetParam();
  j["delta_param"] = delta_param.GetParam();

  j["left_extrapolation"]["toe"] = toe_list_left;
  j["left_extrapolation"]["x"] = xs_err_left;
  j["left_extrapolation"]["y"] = ys_err_left;
  j["left_extrapolation"]["z"] = zs_err_left;
  j["left_extrapolation"]["rms_x"] = rms_x_err_left;
  j["left_extrapolation"]["rms_y"] = rms_y_err_left;
  j["left_extrapolation"]["rms_z"] = rms_z_err_left;

  j["right_extrapolation"]["toe"] = toe_list_right;
  j["right_extrapolation"]["x"] = xs_err_right;
  j["right_extrapolation"]["y"] = ys_err_right;
  j["right_extrapolation"]["z"] = zs_err_right;
  j["right_extrapolation"]["rms_x"] = rms_x_err_right;
  j["right_extrapolation"]["rms_y"] = rms_y_err_right;
  j["right_extrapolation"]["rms_z"] = rms_z_err_right;

  j["interpolation"]["toe"] = toe_list;
  j["interpolation"]["x"] = xs_err;
  j["interpolation"]["y"] = ys_err;
  j["interpolation"]["z"] = zs_err;
  j["interpolation"]["rms_x"] = rms_x_err;
  j["interpolation"]["rms_y"] = rms_y_err;
  j["interpolation"]["rms_z"] = rms_z_err;

  ofstream ofs("result.json");
  ofs << j.dump(4) << endl;
  ofs.close();

  return 0;
}
