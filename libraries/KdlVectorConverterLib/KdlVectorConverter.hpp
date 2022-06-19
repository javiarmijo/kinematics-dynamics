// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#ifndef __KDL_VECTOR_CONVERTER_HPP__
#define __KDL_VECTOR_CONVERTER_HPP__

#include <vector>

#include <kdl/frames.hpp>

/**
 * @ingroup kinematics-dynamics-libraries
 * \defgroup KdlVectorConverterLib
 *
 * @brief Contains classes related to KDL and std::vector classes.
 */

/**
 * @ingroup KdlVectorConverterLib
 * @brief Collection of utilities related to KDL and std::vector classes.
 */
namespace roboticslab::KdlVectorConverter
{

/**
 * @brief Convert from std::vector<double> to KDL::Frame
 *
 * @param x 6-element vector describing a position in cartesian space; first
 * three elements denote translation (meters), last three denote rotation in
 * scaled axis-angle representation (radians).
 *
 * @return Resulting KDL::Frame object.
 */
KDL::Frame vectorToFrame(const std::vector<double> & x);

/**
 * @brief Convert from KDL::Frame to std::vector<double>
 *
 * @param f Input KDL::Frame object.
 *
 * @return Resulting 6-element vector describing a position in cartesian space; first
 * three elements denote translation (meters), last three denote rotation in scaled
 * axis-angle representation (radians).
 */
std::vector<double> frameToVector(const KDL::Frame & f);

/**
 * @brief Convert from std::vector<double> to KDL::Twist
 *
 * @param xdot 6-element vector describing a velocity in cartesian space; first
 * three elements denote translational velocity (meters/second), last three denote
 * angular velocity (radians/second).
 *
 * @return Resulting KDL::Twist object.
 */
KDL::Twist vectorToTwist(const std::vector<double> & xdot);

/**
 * @brief Convert from KDL::Twist to std::vector<double>
 *
 * @param t Input KDL::Twist object
 *
 * @return Resulting 6-element vector describing a velocity in cartesian space; first
 * three elements denote translational velocity (meters/second), last three denote
 * angular velocity (radians/second).
 */
std::vector<double> twistToVector(const KDL::Twist & t);

} // namespace roboticslab::KdlVectorConverter

#endif // __KDL_VECTOR_CONVERTER_HPP__
