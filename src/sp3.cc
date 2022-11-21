/*
 * File: sp3.cc
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Monday, 21st November 2022 10:06:43 am
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */
#include "sp3.h"

#include <fstream>
#include <sstream>

#include "utils.h"

namespace sp3 {
SP3File::SP3File() {}
SP3File::~SP3File() {}

void SP3File::readFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::string line;

  // Read first line
  std::getline(file, line);
  versionIdentifier = line.substr(1, 1).c_str()[0];
  std::string mode = line.substr(2, 1);
  if (mode == "P")
    flag = MODE_FLAG::MODE_P;
  else if (mode == "V")
    flag = MODE_FLAG::MODE_V;
  else
    flag = MODE_FLAG::MODE_UNK;
  int year = stoi(line.substr(3, 4));
  int month = stoi(line.substr(8, 2));
  int day = stoi(line.substr(11, 2));
  int hour = stoi(line.substr(14, 2));
  int minute = stoi(line.substr(17, 2));
  double second = stod(line.substr(20, 11));
  gnsstime = gnsstime::GNSSTime(year, month, day, hour, minute, second,
                                gnsstime::TIME_SYSTEM::GPST);
  epoch = stoi(line.substr(32, 7));
  dataUsed = trim(line.substr(40, 5));
  coordinateSystem = trim(line.substr(46, 5));
  orbitType = trim(line.substr(52, 3));
  agency = trim(line.substr(56, 4));

  // Read second line
  line.clear();
  std::getline(file, line);
  gpsWeek = stoi(line.substr(3, 4));
  gpsSeconds = stod(line.substr(8, 15));
  epochInterval = stof(line.substr(24, 14));
  modifiedJD = stoi(line.substr(39, 5));
  fractionalDay = stod(line.substr(45, 15));

  // Read Satellite IDs
  std::vector<std::string> satNames;
  int lineCount = 0; // line count
  int satCount = 0;  // satellite count

  while (!file.eof()) {
    line.clear();
    std::getline(file, line);
    lineCount++;

    if (lineCount == 1) {
      numberOfSatellites = stoi(line.substr(4, 2));
      satCount = 0;
      int i = 0;
      for (; i < 17 && satCount < numberOfSatellites; i++, satCount++) {
        satNames.push_back(line.substr(9 + 3 * i, 3));
      }
    } else {
      int i = 0;
      for (; i < 17 && satCount < numberOfSatellites; i++, satCount++) {
        satNames.push_back(line.substr(9 + 3 * i, 3));
      }
    }

    if (satCount == numberOfSatellites)
      break;
  }

  // Read Accuracy of Satellite IDs
  lineCount = 0;
  satCount = 0;
  while (!file.eof()) {
    line.clear();
    std::getline(file, line);
    if (line.rfind("++") != 0)
      continue;

    lineCount++;

    if (lineCount == 1)
      satCount = 0;

    int i = 0;
    for (; i < 17 && satCount < numberOfSatellites; i++, satCount++) {
      satInfo.insert(std::pair<std::string, int>(
          satNames.at(satCount), stoi(line.substr(9 + 3 * i, 3))));
    }

    if (satCount == numberOfSatellites)
      break;
  }

  // Move cursor to the beginning of the data
  while (!file.eof()) {
    line.clear();
    std::getline(file, line);
    if (line.rfind("*") == 0)
      break;
  }

  lineCount = 0;
  Epoch epoch;
  double toe = 0;
  int epochCount = 0;
  while (!file.eof()) {
    if (lineCount != 0) {
      line.clear();
      std::getline(file, line);
    }
    lineCount++;

    // If startwith * , then it is a new epoch
    if (line.rfind("*") == 0) {
      if (epochCount > 0) {
        ephemerisData.push_back(epoch);
      }
      epochCount++;
      epoch = Epoch();
      int year = stoi(line.substr(3, 4));
      int month = stoi(line.substr(8, 2));
      int day = stoi(line.substr(11, 2));
      int hour = stoi(line.substr(14, 2));
      int minute = stoi(line.substr(17, 2));
      float second = stof(line.substr(20, 11));
      epoch.time = gnsstime::GNSSTime(year, month, day, hour, minute, second,
                                      gnsstime::TIME_SYSTEM::GPST);
      int week;
      epoch.time.ToGPSTime(&week, &toe);
    }

    // If startwith EOF, then it is the end of the file
    else if (line.rfind("EOF") == 0) {
      ephemerisData.push_back(epoch);
      break;
    }

    else {
      PosAndClock record;
      std::string satID = line.substr(1, 3);
      record.toe = toe;
      record.x = stod(line.substr(4, 14));
      record.y = stod(line.substr(18, 14));
      record.z = stod(line.substr(32, 14));
      record.clock = stod(line.substr(46, 14));

      epoch.data.insert(std::pair<std::string, PosAndClock>(satID, record));
    }
  }
}

std::vector<PosAndClock>
SP3File::GetKNearestNeighbors(const std::string &satName,
                              const gnsstime::GNSSTime &time, int k) {
  std::vector<PosAndClock> result;

  std::vector<double> tks;
  for (auto &epoch : ephemerisData) {
    tks.push_back(abs(epoch.time - time));
  }

  std::vector<int> indices = GetKSmallestIndices(tks, k);
  for (auto &index : indices) {
    result.push_back(ephemerisData.at(index).data.at(satName));
  }

  return result;
}

PosAndClock SP3File::GetInterpolatedData(const std::string &satName,
                                         gnsstime::GNSSTime &time, int k) {
  std::vector<PosAndClock> neighbors = GetKNearestNeighbors(satName, time, k);
  std::vector<double> ts, xs, ys, zs, clocks;
  for (auto &neighbor : neighbors) {
    ts.push_back(neighbor.toe);
    xs.push_back(neighbor.x);
    ys.push_back(neighbor.y);
    zs.push_back(neighbor.z);
    clocks.push_back(neighbor.clock);
  }

  int week;
  double toe;
  time.ToGPSTime(&week, &toe);

  double x = LagrangeInterpolate(ts, xs, toe);
  double y = LagrangeInterpolate(ts, ys, toe);
  double z = LagrangeInterpolate(ts, zs, toe);
  double clock = LagrangeInterpolate(ts, clocks, toe);

  PosAndClock result{x, y, z, clock, toe};
  return result;
}

} // namespace sp3