#include "WiimoteSensorDevice.hpp"

#include <cmath>
#include <algorithm>

#include <yarp/os/LogStream.h>
#include <yarp/sig/Vector.h>

#include "LogComponent.hpp"

using namespace roboticslab;

WiimoteSensorDevice::WiimoteSensorDevice(yarp::os::Searchable & config, bool usingPose)
    : StreamingDevice(config),
      iAnalogSensor(nullptr),
      mode(NONE),
      usingPose(usingPose),
      step(0.0)
{
    data.resize(3);  // already called by base constructor
    buffer.resize(5);
    step = config.check("step", yarp::os::Value(DEFAULT_STEP), "").asFloat64();
}

bool WiimoteSensorDevice::acquireInterfaces()
{
    bool ok = true;

    if (!yarp::dev::PolyDriver::view(iAnalogSensor))
    {
        yCWarning(SDC) << "Could not view iAnalogSensor";
        ok = false;
    }

    return ok;
}

bool WiimoteSensorDevice::initialize(bool usingStreamingPreset)
{
    if (usingPose && step <= 0.0)
    {
        yCWarning(SDC) << "Invalid step:" << step;
        return false;
    }

    if (usingStreamingPreset)
    {
        int cmd = usingPose ? VOCAB_CC_POSE : VOCAB_CC_TWIST;

        if (!iCartesianControl->setParameter(VOCAB_CC_CONFIG_STREAMING_CMD, cmd))
        {
            yCWarning(SDC) << "Unable to preset streaming command";
            return false;
        }
    }

    if (!iCartesianControl->setParameter(VOCAB_CC_CONFIG_FRAME, ICartesianSolver::TCP_FRAME))
    {
        yCWarning(SDC) << "Unable to set TCP reference frame";
        return false;
    }

    return true;
}

bool WiimoteSensorDevice::acquireData()
{
    yarp::sig::Vector data;
    iAnalogSensor->read(data);

    yCDebug(SDC) << data.toString(4, 1);

    if (data.size() != 5)
    {
        yCWarning(SDC) << "Invalid data size:" << data.size();
        return false;
    }

    for (int i = 0; i < data.size(); i++)
    {
        buffer[i] = data[i];
    }

    return true;
}

bool WiimoteSensorDevice::transformData(double scaling)
{
    bool buttonA = buffer[2] == 1.0;
    bool buttonB = buffer[3] == 1.0;
    bool yawActive = buffer[4] == 1.0;

    if (buttonA && buttonB)
    {
        mode = ROT;
    }
    else if (buttonA)
    {
        mode = FWD;
    }
    else if (buttonB)
    {
        mode = BKWD;
    }
    else
    {
        mode = NONE;
        return true;
    }

    data[1] = -buffer[1] / scaling;

    if (yawActive)
    {
        data[0] = 0.0;
        data[2] = buffer[0] / scaling;
    }
    else
    {
        data[0] = buffer[0] / scaling;
        data[2] = 0.0;
    }

    return true;
}

bool WiimoteSensorDevice::hasValidMovementData() const
{
    return mode != NONE;
}

void WiimoteSensorDevice::sendMovementCommand(double timestamp)
{
    std::vector<double> xdot(6, 0.0);
    std::copy(data.begin(), data.end(), xdot.begin() + 3);

    switch (mode)
    {
    case FWD:
        xdot[2] = step;
        break;
    case BKWD:
        xdot[2] = -step;
        break;
    case ROT:
        xdot[2] = 0.0; // for the sake of explicitness
        break;
    default:
        return;
    }

    if (usingPose)
    {
        iCartesianControl->pose(xdot);
    }
    else
    {
        iCartesianControl->twist(xdot);
    }
}

void WiimoteSensorDevice::stopMotion()
{
    if (!usingPose)
    {
        std::vector<double> zeros(6, 0.0);
        iCartesianControl->twist(zeros);
    }
}
