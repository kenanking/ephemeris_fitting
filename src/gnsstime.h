/*
 * File: gnsstime.h
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 2:39:58 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#ifndef GNSSTIME_H
#define GNSSTIME_H

#include <ctime>
#include <iostream>
#include <string>

namespace gnsstime {

const static double GPST_0[] = {1980, 1, 6, 0, 0, 0}; // GPS start time
const static double BDST_0[] = {2006, 1, 1, 0, 0, 0}; // BDS start time

const static int MAX_LEAPS = 64; // Max number of leap seconds table

const static double LEAPS_TABLE[MAX_LEAPS + 1]
                               [7] = // Leap Seconds (Y,M,D,h,m,s,utc-gpst)
    {{2017, 1, 1, 0, 0, 0, -18},
     {2015, 7, 1, 0, 0, 0, -17},
     {2012, 7, 1, 0, 0, 0, -16},
     {2009, 1, 1, 0, 0, 0, -15},
     {2006, 1, 1, 0, 0, 0, -14},
     {1999, 1, 1, 0, 0, 0, -13},
     {1997, 7, 1, 0, 0, 0, -12},
     {1996, 1, 1, 0, 0, 0, -11},
     {1994, 7, 1, 0, 0, 0, -10},
     {1993, 7, 1, 0, 0, 0, -9},
     {1992, 7, 1, 0, 0, 0, -8},
     {1991, 1, 1, 0, 0, 0, -7},
     {1990, 1, 1, 0, 0, 0, -6},
     {1988, 1, 1, 0, 0, 0, -5},
     {1985, 7, 1, 0, 0, 0, -4},
     {1983, 7, 1, 0, 0, 0, -3},
     {1982, 7, 1, 0, 0, 0, -2},
     {1981, 7, 1, 0, 0, 0, -1},
     {0}};

enum TIME_SYSTEM {
  UTC = 0x01,  // UTC time
  GPST = 0x02, // GPS time
  BDST = 0x04, // BDS time
};

void Time2Sec(int year, int month, int day, int hour, int minute, double second,
              time_t *time, double *second_fraction);
void Sec2Time(time_t time, double second_fraction, int *year, int *month,
              int *day, int *hour, int *minute, double *second);

class GNSSTime {

public:
  GNSSTime();
  GNSSTime(int year, int month, int day, int hour, int minute, double second,
           TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  GNSSTime(const double time[6], TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  GNSSTime(time_t time, double second_fraction,
           TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  GNSSTime(const GNSSTime &t);

  ~GNSSTime();

  inline time_t GetSecond() const { return second_; };
  inline double GetSecondFraction() const { return second_fraction_; };

  void SetTime(int year, int month, int day, int hour, int minute,
               double second, TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  void SetTime(time_t time, double second_fraction,
               TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  void GetTime(int *year, int *month, int *day, int *hour, int *minute,
               double *second, TIME_SYSTEM time_system = TIME_SYSTEM::UTC);
  void GetTime(time_t *time, double *second_fraction,
               TIME_SYSTEM time_system = TIME_SYSTEM::UTC);

  void ToGPSTime(int *week, double *seconds_of_week);

  std::string ToString(TIME_SYSTEM time_system = TIME_SYSTEM::UTC);

  // Overload operators
  GNSSTime &operator=(const GNSSTime &t);

  GNSSTime operator+(const double &sec);
  GNSSTime operator-(const double &sec);

  double operator-(const GNSSTime &t);

  bool operator==(const GNSSTime &other);
  bool operator!=(const GNSSTime &other);
  bool operator<(const GNSSTime &other);
  bool operator>(const GNSSTime &other);
  bool operator<=(const GNSSTime &other);
  bool operator>=(const GNSSTime &other);

private:
  time_t second_ = 0; // Time in seconds since 1970-01-01 00:00:00 UTC
  double second_fraction_ = 0.0;  // Second fraction
  TIME_SYSTEM time_system_ = UTC; // Time system

  void GPST2UTC(time_t &time, double &second_fraction);
  void UTC2GPST(time_t &time, double &second_fraction);

  void GPST2BDST(time_t &time, double &second_fraction);
  void BDST2GPST(time_t &time, double &second_fraction);
};

inline GNSSTime TimeAdd(GNSSTime time, double second) {
  // Add seconds to time
  double second_fraction = time.GetSecondFraction() + second;
  int tt = (int)floor(second_fraction);

  return GNSSTime(time.GetSecond() + tt, second_fraction - tt);
}

inline double TimeDiff(GNSSTime t1, GNSSTime t2) {
  // Calculate the difference between two times in seconds
  return difftime(t1.GetSecond(), t2.GetSecond()) + t1.GetSecondFraction() -
         t2.GetSecondFraction();
}

}; // namespace gnsstime

#endif // GNSSTIME_H