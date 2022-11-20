/*
 * File: gnsstime.cc
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 2:39:53 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "gnsstime.h"

namespace gnsstime {

void Time2Sec(int year, int month, int day, int hour, int minute, double second,
              time_t *time, double *second_fraction) {
  // Convert time to seconds since 1970-01-01 00:00:00 UTC

  // Check if the input time is valid
  if (year < 1970 || year > 2099 || month < 1 || month > 12 || day < 1 ||
      day > 31 || hour < 0 || hour > 23 || minute < 0 || minute > 59 ||
      second < 0.0 || second >= 60.0) {
    throw std::invalid_argument("Invalid time");
  }

  const int doy[] = {1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
  int days = (year - 1970) * 365 + (year - 1969) / 4 + doy[month - 1] + day - 2;
  if (year % 4 == 0 && month > 2) {
    // Leap year if year % 4 == 0 and month > 2 in 1970-2099
    days++;
  }
  int sec = (int)floor(second);

  *time = (time_t)days * 86400 + hour * 3600 + minute * 60 + sec;
  *second_fraction = second - sec;
}

void Sec2Time(time_t time, double second_fraction, int *year, int *month,
              int *day, int *hour, int *minute, double *second) {
  // Convert seconds since 1970-01-01 00:00:00 UTC to time

  // Check if the input time is valid
  if (time < 0 || second_fraction < 0.0 || second_fraction >= 1.0) {
    throw std::invalid_argument("Invalid time");
  }

  const int day_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
                              31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
                              31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
                              31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int tmp_days, tmp_sec, tmp_mon, tmp_day;

  tmp_days = (int)(time / 86400);
  tmp_sec = (int)(time % 86400);
  for (tmp_day = tmp_days % 1461, tmp_mon = 0; tmp_mon < 48; tmp_mon++) {
    if (tmp_day < day_in_month[tmp_mon]) {
      break;
    }
    tmp_day -= day_in_month[tmp_mon];
  }

  *year = 1970 + tmp_days / 1461 * 4 + tmp_mon / 12;
  *month = tmp_mon % 12 + 1;
  *day = tmp_day + 1;
  *hour = tmp_sec / 3600;
  *minute = (tmp_sec % 3600) / 60;
  *second = (tmp_sec % 60) + second_fraction;
}

GNSSTime::GNSSTime() {
  // Default constructor
  second_ = 0;
  second_fraction_ = 0.0;
}

GNSSTime::GNSSTime(int year, int month, int day, int hour, int minute,
                   double second, TIME_SYSTEM time_system) {
  // Constructor
  SetTime(year, month, day, hour, minute, second, time_system);
}

GNSSTime::GNSSTime(const double time[6], TIME_SYSTEM time_system) {
  // Constructor
  SetTime(int(time[0]), int(time[1]), int(time[2]), int(time[3]), int(time[4]),
          time[5]);
}

GNSSTime::GNSSTime(time_t time, double second_fraction,
                   TIME_SYSTEM time_system) {
  // Constructor
  second_ = time;
  second_fraction_ = second_fraction;
}

GNSSTime::GNSSTime(const GNSSTime &t) {
  // Copy constructor
  second_ = t.second_;
  second_fraction_ = t.second_fraction_;
}

GNSSTime &GNSSTime::operator=(const GNSSTime &t) {
  // Assignment operator
  second_ = t.second_;
  second_fraction_ = t.second_fraction_;
  return *this;
}

GNSSTime::~GNSSTime() {
  // Destructor
}

void GNSSTime::SetTime(int year, int month, int day, int hour, int minute,
                       double second, TIME_SYSTEM time_system) {
  // Set time
  time_t time;
  double second_fraction;
  Time2Sec(year, month, day, hour, minute, second, &time, &second_fraction);
  SetTime(time, second_fraction, time_system);
}

void GNSSTime::SetTime(time_t time, double second_fraction,
                       TIME_SYSTEM time_system) {
  // Set time
  second_ = time;
  second_fraction_ = second_fraction;

  switch (time_system) {
  case TIME_SYSTEM::UTC:
    break;
  case TIME_SYSTEM::GPST:
    GPST2UTC(second_, second_fraction_);
    break;
  case TIME_SYSTEM::BDST:
    break;
  default:
    break;
  }
}

void GNSSTime::GetTime(int *year, int *month, int *day, int *hour, int *minute,
                       double *second, TIME_SYSTEM time_system) {
  // Get time
  time_t time = second_;
  double second_fraction = second_fraction_;
  GetTime(&time, &second_fraction, time_system);
  Sec2Time(time, second_fraction, year, month, day, hour, minute, second);
}

void GNSSTime::GetTime(time_t *time, double *second_fraction,
                       TIME_SYSTEM time_system) {
  // Get time
  *time = second_;
  *second_fraction = second_fraction_;

  switch (time_system) {
  case TIME_SYSTEM::UTC:
    break;
  case TIME_SYSTEM::GPST:
    UTC2GPST(*time, *second_fraction);
    break;
  case TIME_SYSTEM::BDST:
    break;
  default:
    break;
  }
}

void GNSSTime::UTC2GPST(time_t &time, double &second_fraction) {
  // Convert UTC time to GPS time
  GNSSTime t = GNSSTime(time, second_fraction, TIME_SYSTEM::UTC);

  for (int i = 0; LEAPS_TABLE[i][0] > 0; i++) {
    if (TimeDiff(t, GNSSTime(LEAPS_TABLE[i])) >= 0.0) {
      t = TimeAdd(t, -LEAPS_TABLE[i][6]);
      break;
    }
  }
  time = t.GetSecond();
  second_fraction = t.GetSecondFraction();
}

void GNSSTime::GPST2UTC(time_t &time, double &second_fraction) {
  // Convert GPS time to UTC time
  GNSSTime t = GNSSTime(time, second_fraction, TIME_SYSTEM::GPST);

  for (int i = 0; LEAPS_TABLE[i][0] > 0; i++) {
    t = TimeAdd(t, LEAPS_TABLE[i][6]);
    if (TimeDiff(t, GNSSTime(LEAPS_TABLE[i])) >= 0.0) {
      break;
    }
  }
  time = t.GetSecond();
  second_fraction = t.GetSecondFraction();
}

void GNSSTime::GPST2BDST(time_t &time, double &second_fraction) {
  // Convert GPS time to BDS time
  GNSSTime t = GNSSTime(time, second_fraction, TIME_SYSTEM::GPST);
  t = TimeAdd(t, -14.0);
  time = t.GetSecond();
  second_fraction = t.GetSecondFraction();
}

void GNSSTime::BDST2GPST(time_t &time, double &second_fraction) {
  // Convert BDS time to GPS time
  GNSSTime t = GNSSTime(time, second_fraction, TIME_SYSTEM::BDST);
  t = TimeAdd(t, 14.0);
  time = t.GetSecond();
  second_fraction = t.GetSecondFraction();
}

void GNSSTime::ToGPSTime(int *week, double *seconds_of_week) {
  // Convert time to GPS time
  GNSSTime t = GNSSTime(GPST_0);

  time_t second_tmp = second_;
  double second_fraction_tmp = second_fraction_;
  UTC2GPST(second_tmp, second_fraction_tmp);

  time_t sec = second_tmp - t.GetSecond();
  *week = int(sec / (7 * 24 * 3600));
  *seconds_of_week = (double)(sec - *week * 7 * 24 * 3600) + second_fraction_tmp;
}

std::string GNSSTime::ToString(TIME_SYSTEM time_system) {
  // Convert time to string
  int year, month, day, hour, minute;
  double second;
  GetTime(&year, &month, &day, &hour, &minute, &second, time_system);

  std::ostringstream oss;
  oss << std::setfill('0') << std::setw(4) << year << "-" << std::setw(2)
      << month << "-" << std::setw(2) << day << " " << std::setw(2) << hour
      << ":" << std::setw(2) << minute << ":" << std::setw(10)
      << std::setprecision(7) << std::fixed << second;
  return oss.str();
}

GNSSTime GNSSTime::operator+(const double &sec) {
  // Add time
  GNSSTime time = *this;
  time = TimeAdd(time, sec);
  return time;
}

GNSSTime GNSSTime::operator-(const double &sec) {
  // Subtract time
  GNSSTime time = *this;
  time = TimeAdd(time, -sec);
  return time;
}

double GNSSTime::operator-(const GNSSTime &other) {
  return (this->second_ - other.second_) +
         (this->second_fraction_ - other.second_fraction_);
}

bool GNSSTime::operator==(const GNSSTime &other) {
  if (abs(*this - other) <= 1E-10)
    return true;
  else
    return false;
}

bool GNSSTime::operator!=(const GNSSTime &other) {
  if (abs(*this - other) <= 1E-10)
    return false;
  else
    return true;
}

bool GNSSTime::operator>(const GNSSTime &other) {
  if ((*this - other) > 0)
    return true;
  else
    return false;
}

bool GNSSTime::operator>=(const GNSSTime &other) {
  if ((*this - other) >= 0)
    return true;
  else
    return false;
}

bool GNSSTime::operator<(const GNSSTime &other) {
  if ((*this - other) < 0)
    return true;
  else
    return false;
}

bool GNSSTime::operator<=(const GNSSTime &other) {
  if ((*this - other) <= 0)
    return true;
  else
    return false;
}

} // namespace gnsstime
