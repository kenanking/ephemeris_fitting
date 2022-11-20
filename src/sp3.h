/*
 * File: sp3.h
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 2:40:45 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#ifndef SP3_H
#define SP3_H

#include <map>
#include <string>
#include <vector>

#include "gnsstime.h"

namespace sp3 {

enum MODE_FLAG {
  MODE_P = 0x01,   // Position
  MODE_V = 0x02,   // Velocity
  MODE_UNK = 0x04, // Unknown
};

struct PosAndClock {
  double toe;     // time of ephemeris (s)
  double x, y, z; // x, y, z coordinates (km)
  double clock;   // clock time (10^-6 s)
};

struct Epoch {
  gnsstime::GNSSTime time;
  std::map<std::string, PosAndClock> data;
};

class SP3File {
public:
  char versionIdentifier;       // Version Identifier
  MODE_FLAG flag;               // Position/Velocity Mode Flag
  gnsstime::GNSSTime gnsstime;  // Start Year, Month, Day, Hour, Minute, Second
  int epoch;                    // Number of Epochs
  std::string dataUsed;         // Data Used
  std::string coordinateSystem; // Coordinate System
  std::string orbitType;        // Orbit Type
  std::string agency;           // Array
  int gpsWeek;                  // GPS Week
  double gpsSeconds;            // GPS Seconds of Week
  double epochInterval;         // Epoch Interval
  int modifiedJD;               // Modified Julian Day Start
  double fractionalDay;         // Fractional Day
  int numberOfSatellites;       // Number of Satellites
  std::map<std::string, int> satInfo; // Satellite IDs and its Accuracy

  std::vector<Epoch> ephemerisData; // Ephemeris Data (Epochs)

  SP3File();
  ~SP3File();

  void readFile(const std::string &filename);

  std::vector<PosAndClock> GetKNearestNeighbors(const std::string &satName,
                                                const gnsstime::GNSSTime &time,
                                                int k);
};

}; // namespace sp3

#endif // SP3_H