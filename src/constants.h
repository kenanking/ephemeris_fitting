/*
 * File: constants.h
 * File Created: Wednesday, 16th November 2022 10:19:47 pm
 * Author: Yan Tang (360383464@qq.com)
 * -----
 * Last Modified: Sunday, 20th November 2022 2:39:45 pm
 * Modified By: Yan Tang (360383464@qq.com>)
 * -----
 * Copyright 2022 - 2022 Yan Tang
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

const double light_spd = 2.99792458E+8; // Speed of light in vacuum (m/s)

const double pi = 3.1415926535897932385; // pi -- Pi

const double e = 2.7182818284590452354; // e  -- Euler's Constant

namespace WGS84 {
/**
 * @brief Constant Parameters for WGS-84 reference ellipsoid.
 *
 */

const double GM = 3.9860050E+14; // Earth's gravitational constant (m^3/s^2)
const double OMEGA_E = 7.2921151467E-5; // Angular velocity of the earth (rad/s)
const double a = 6378137.0;             // Semi-major axis of the ellipsoid (m)
const double f = 1 / 298.257223563;     // Flattening of the ellipsoid
} // namespace WGS84

#endif // CONSTANTS_H