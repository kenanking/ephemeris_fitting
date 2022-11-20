/*
 * File: main.cc
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 11:49:25 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

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
  param.sqrtA = 5153.661;
  param.e = 0.01201;
  param.omega = 0.94435;
  param.OMEGA = 0.36391;
  param.i0 = 0.98890;
  param.M0 = 2.82278;

  // 拟合相关设置
  FittingOptions options;
  options.maxIter = 150; // 最大迭代次数
  options.tol = 1e-5; // 迭代收敛阈值，前后两次中误差的变化小于该值时
  options.deltaF = 0;        // ΔF采用米级精度
  options.delta = 1e-6;      // 计算toe时刻速度使用的微小量
  options.fitKepler = true; // 是否使用Kepler初轨作为初始值
  options.fixDeltaOrbitParam = false; // 是否固定数值求导的微小量（在迭代中更新）
  options.verbose = true; // 输出迭代过程

  // 设置数值导数计算时的微小量初值
  DeltaOrbitParam delta_param = DeltaOrbitParam();

  // 拟合星历
  cout << "========================== Fitting Orbit ========================="
       << endl;

  FitEphemeris(sp3Data, param, delta_param, options);

  cout << "================================================================="
       << endl;
  cout << endl;

  // 输出拟合结果
  cout << "========================== Fitting Result ========================"
       << endl;
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
  cout << endl;
  cout << "================================================================="
       << endl;
  cout << endl;

  return 0;
}
