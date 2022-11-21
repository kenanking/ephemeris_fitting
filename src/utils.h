/*
 * File: utils.h
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Monday, 21st November 2022 9:56:56 am
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#ifndef UTILS_H
#define UTILS_H

#include <queue>
#include <string>
#include <vector>

#include <Eigen/Dense>

inline void ltrimHelper(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));
}

inline void rtrimHelper(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());
}

inline std::string ltrim(std::string s) {
  ltrimHelper(s);
  return s;
}

inline std::string rtrim(std::string s) {
  rtrimHelper(s);
  return s;
}

inline std::string trim(std::string s) {
  ltrimHelper(s);
  rtrimHelper(s);
  return s;
}

inline std::vector<int> GetKSmallestIndices(const std::vector<double> &v,
                                            int k) {
  std::priority_queue<double> pq;
  for (auto d : v) {
    if (pq.size() >= (size_t)k && pq.top() > d) {
      pq.push(d);
      pq.pop();
    } else if (pq.size() < (size_t)k) {
      pq.push(d);
    }
  }

  double smallest = pq.top();

  std::vector<int> indices;
  for (int i = 0; i < v.size(); i++) {
    if (v[i] <= smallest) {
      indices.push_back(i);
    }
  }

  return indices;
}

inline double LagrangeInterpolate(std::vector<double> &x,
                                  std::vector<double> &y, double x0) {
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

#endif // UTILS_H