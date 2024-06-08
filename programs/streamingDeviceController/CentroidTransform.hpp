#ifndef __CENTROID_TRANSFORM_HPP__
#define __CENTROID_TRANSFORM_HPP__

#include <yarp/os/Bottle.h>
#include <yarp/os/Stamp.h>

#include <kdl/frames.hpp>

#include "StreamingDevice.hpp"

namespace roboticslab
{

class StreamingDevice;

/**
 * @ingroup streamingDeviceController
 *
 * @see @cite eona2020icarsc
 */
class CentroidTransform
{
public:
    //! Constructor
    CentroidTransform();

    //! Register handle to device
    void registerStreamingDevice(StreamingDevice * streamingDevice)
    { this->streamingDevice = streamingDevice; }

    //! Set TCP to camera frame
    bool setTcpToCameraRotation(yarp::os::Bottle * b);

    //! Set new permanence time
    void setPermanenceTime(double permanenceTime)
    { this->permanenceTime = permanenceTime; }

    //! Register or dismiss incoming bottle
    bool acceptBottle(yarp::os::Bottle * b);

    //! Process last stored bottle
    bool processStoredBottle() const;

private:
    StreamingDevice * streamingDevice;
    double permanenceTime;
    yarp::os::Bottle lastBottle;
    yarp::os::Stamp lastAcquisition;
    KDL::Rotation rot_tcp_camera;
};

} // namespace roboticslab

#endif // __CENTROID_TRANSFORM_HPP__
