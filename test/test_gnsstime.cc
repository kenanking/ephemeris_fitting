/*
 * File: test_gnsstime.cc
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 2:41:34 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#include <cassert>
#include <iostream>

#include "gnsstime.h"

using namespace gnsstime;

void TestTime2Sec() {
  time_t time;
  double second_fraction;
  Time2Sec(2022, 11, 16, 12, 55, 13.123, &time, &second_fraction);
  assert(time == 1668603313);
  assert(abs(second_fraction - 0.123) < 1e-6);
  std::cout << "TestTime2Sec_1 passed" << std::endl;

  Time2Sec(2023, 11, 16, 12, 55, 13.123456, &time, &second_fraction);
  assert(time == 1700139313);
  assert(abs(second_fraction - 0.123456) < 1e-6);
  std::cout << "TestTime2Sec_2 passed" << std::endl;

  Time2Sec(1999, 12, 31, 23, 59, 59.9, &time, &second_fraction);
  assert(time == 946684799);
  assert(abs(second_fraction - 0.9) < 1e-6);
  std::cout << "TestTime2Sec_3 passed" << std::endl;
}

void TestSec2Time() {
  int year, month, day, hour, minute;
  double second;
  Sec2Time(1668603313, 0.123, &year, &month, &day, &hour, &minute, &second);
  assert(year == 2022);
  assert(month == 11);
  assert(day == 16);
  assert(hour == 12);
  assert(minute == 55);
  assert(abs(second - 13.123) < 1e-6);
  std::cout << "TestSec2Time_1 passed" << std::endl;

  Sec2Time(1700139313, 0.123456, &year, &month, &day, &hour, &minute, &second);
  assert(year == 2023);
  assert(month == 11);
  assert(day == 16);
  assert(hour == 12);
  assert(minute == 55);
  assert(abs(second - 13.123456) < 1e-6);
  std::cout << "TestSec2Time_2 passed" << std::endl;

  Sec2Time(946684799, 0.9, &year, &month, &day, &hour, &minute, &second);
  assert(year == 1999);
  assert(month == 12);
  assert(day == 31);
  assert(hour == 23);
  assert(minute == 59);
  assert(abs(second - 59.9) < 1e-6);
  std::cout << "TestSec2Time_3 passed" << std::endl;
}

void TestGNSSTime() {
  GNSSTime gnsstime1(2022, 11, 16, 12, 55, 13.123);
  int year, month, day, hour, minute;
  double second;
  gnsstime1.GetTime(&year, &month, &day, &hour, &minute, &second);
  assert(year == 2022);
  assert(month == 11);
  assert(day == 16);
  assert(hour == 12);
  assert(minute == 55);
  assert(abs(second - 13.123) < 1e-6);

  gnsstime1.GetTime(&year, &month, &day, &hour, &minute, &second,
                    gnsstime::GPST);
  assert(year == 2022);
  assert(month == 11);
  assert(day == 16);
  assert(hour == 12);
  assert(minute == 55);
  assert(abs(second - 31.123) < 1e-6);
}

void TestToGPST() {
  // Ground truth from 

  GNSSTime gt = GNSSTime(2022, 10, 1, 6, 0, 0, TIME_SYSTEM::UTC);
  int week;
  double seconds_of_week;

  gt.ToGPSTime(&week, &seconds_of_week);
  assert(week == 2229);
  assert(abs(seconds_of_week - 540018) < 1e-6);

  gt = GNSSTime(2000, 7, 22, 8, 6, 19, TIME_SYSTEM::UTC);
  gt.ToGPSTime(&week, &seconds_of_week);
  assert(week == 1071);
  assert(abs(seconds_of_week - 547592) < 1e-6);
}

int main() {
  std::cout << "Test gnsstime.cc" << std::endl;

  std::cout << "TestTime2Sec" << std::endl;
  TestTime2Sec();
  std::cout << std::endl;

  std::cout << "TestSec2Time" << std::endl;
  TestSec2Time();
  std::cout << std::endl;

  std::cout << "TestGNSSTime" << std::endl;
  TestGNSSTime();
  std::cout << std::endl;

  std::cout << "TestToGPST" << std::endl;
  TestToGPST();
  std::cout << std::endl;

  return 0;
}
