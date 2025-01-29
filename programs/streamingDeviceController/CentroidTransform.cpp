#include "CentroidTransform.hpp"

#include <cmath>

#include <yarp/os/LogStream.h>
#include <yarp/os/Time.h>

#include <kdl/frames.hpp>

#include <KdlVectorConverter.hpp>

#include "LogComponent.hpp"

using namespace roboticslab;

constexpr auto ROT_FACTOR = 0.1;

CentroidTransform::CentroidTransform()
    : streamingDevice(nullptr),
      permanenceTime(0.0)
{}

bool CentroidTransform::setTcpToCameraRotation(yarp::os::Bottle * b)
{
    if (b->size() != 3)
    {
        yCWarning(SDC) << "Bottle size must equal 3, was:" << b->size();
        return false;
    }

    double roll = b->get(0).asFloat64() * M_PI / 180.0;
    double pitch = b->get(1).asFloat64() * M_PI / 180.0;
    double yaw = b->get(2).asFloat64() * M_PI / 180.0;

    yCInfo(SDC) << "centroidFrameRPY [rad]:" << roll << pitch << yaw;

    rot_tcp_camera = KDL::Rotation::RPY(roll, pitch, yaw);

    return true;
}

bool CentroidTransform::acceptBottle(yarp::os::Bottle * b)
{
    if (b)
    {
        if (b->size() != 2)
        {
            yCWarning(SDC) << "Malformed input bottle, size" << b->size() << "(expected 2)";
            return false;
        }

        lastBottle = *b;
        lastAcquisition.update();
        return true;
    }

    return yarp::os::Time::now() - lastAcquisition.getTime() <= permanenceTime;
}

bool CentroidTransform::processStoredBottle() const
{
    // object centroids scaled to fit into [-1, 1] range
    double cx = lastBottle.get(0).asFloat64(); // points right
    double cy = lastBottle.get(1).asFloat64(); // points down

    std::vector<double> x;

    if (!streamingDevice->iCartesianControl->stat(x))
    {
        yCWarning(SDC) << "stat failed";
        return false;
    }

    KDL::Frame H_base_tcp = KdlVectorConverter::vectorToFrame(x);

    // express camera's z axis (points "forward") in base frame
    KDL::Vector v_base = H_base_tcp.M * rot_tcp_camera * KDL::Vector(0, 0, 1);
    KDL::Frame H_base_target = KdlVectorConverter::vectorToFrame(streamingDevice->data);

    double norm = KDL::dot(H_base_target.p, v_base);

    if (norm <= 0.0)
    {
        // no action if we move away from the target (negative TCP's z axis)
        return false;
    }

    // project target vector into TCP's z axis, refer result to base frame
    H_base_target.p = v_base * norm;

    // find axis along which to rotate (in TCP frame) given pixel coords
    KDL::Vector coords(cx, cy, 0);
    KDL::Vector tcp_axis = KDL::Rotation::RotZ(KDL::PI_2) * coords;
    KDL::Vector base_axis = H_base_tcp.M * rot_tcp_camera * tcp_axis;

    // rotate towards the target in base frame
    H_base_target.M = KDL::Rotation::Rot(base_axis, coords.Norm() * ROT_FACTOR);

    // apply changes to input transform
    streamingDevice->data = KdlVectorConverter::frameToVector(H_base_target);

    return true;
}
