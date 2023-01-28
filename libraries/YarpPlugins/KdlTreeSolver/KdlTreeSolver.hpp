// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#ifndef __KDL_TREE_SOLVER_HPP__
#define __KDL_TREE_SOLVER_HPP__

#include <map>
#include <string>
#include <vector>

#include <yarp/dev/DeviceDriver.h>

#include <kdl/tree.hpp>
#include <kdl/treefksolver.hpp>
#include <kdl/treeiksolver.hpp>
#include <kdl/treeidsolver.hpp>

#include "ICartesianSolver.h"

namespace roboticslab
{

/**
 * @ingroup YarpPlugins
 * @defgroup KdlTreeSolver
 *
 * @brief Contains roboticslab::KdlSolver.
 */

/**
 * @ingroup KdlTreeSolver
 * @brief The KdlTreeSolver class implements ICartesianSolver.
 */
class KdlTreeSolver : public yarp::dev::DeviceDriver,
                      public ICartesianSolver
{
public:
    KdlTreeSolver() : fkSolverPos(nullptr),
                      ikSolverPos(nullptr),
                      ikSolverVel(nullptr),
                      idSolver(nullptr)
    {}

    // -- ICartesianSolver declarations. Implementation in ICartesianSolverImpl.cpp --

    // Get number of joints for which the solver has been configured.
    int getNumJoints() override;

    // Get number of TCPs for which the solver has been configured.
    int getNumTcps() override;

    // Append an additional link.
    bool appendLink(const std::vector<double> & x) override;

    // Restore original kinematic chain.
    bool restoreOriginalChain() override;

    // Change reference frame.
    bool changeOrigin(const std::vector<double> & x_old_obj, const std::vector<double> & x_new_old, std::vector<double> & x_new_obj) override;

    // Perform forward kinematics.
    bool fwdKin(const std::vector<double> & q, std::vector<double> & x) override;

    // Obtain difference between supplied pose inputs.
    bool poseDiff(const std::vector<double> & xLhs, const std::vector<double> & xRhs, std::vector<double> & xOut) override;

    // Perform inverse kinematics.
    bool invKin(const std::vector<double> & xd, const std::vector<double> & qGuess, std::vector<double> & q, const reference_frame frame) override;

    // Perform differential inverse kinematics.
    bool diffInvKin(const std::vector<double> & q, const std::vector<double> & xdot, std::vector<double> & qdot, const reference_frame frame) override;

    // Perform inverse dynamics.
    bool invDyn(const std::vector<double> & q, std::vector<double> & t) override;

    // Perform inverse dynamics.
    bool invDyn(const std::vector<double> & q, const std::vector<double> & qdot, const std::vector<double> & qdotdot, const std::vector<double> & ftip, std::vector<double> & t) override;

    // -------- DeviceDriver declarations. Implementation in IDeviceImpl.cpp --------

    bool open(yarp::os::Searchable & config) override;

    bool close() override;

protected:
    std::vector<std::string> endpoints;
    std::map<std::string, std::string> mergedEndpoints;
    KDL::Tree tree;
    KDL::TreeFkSolverPos * fkSolverPos;
    KDL::TreeIkSolverPos * ikSolverPos;
    KDL::TreeIkSolverVel * ikSolverVel;
    KDL::TreeIdSolver * idSolver;
};

} // namespace roboticslab

#endif // __KDL_TREE_SOLVER_HPP__
