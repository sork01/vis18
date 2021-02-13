// *neode.onsave* setgo g++ -I../ -I../../../inviwo/include -I../../../inviwo/modules -I../../../inviwo/modules/_generated -I../../../inviwo/ext -fsyntax-only integrator.h
/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labtopo/labtopomoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <functional>

namespace inviwo {

class IVW_MODULE_LABTOPO_API Integrator {
// Friends
// Types
public:
    struct StreamlineOptions
    {
        bool enableStepLimit = false;
        bool enableArcLengthLimit = false;
        bool enableLowVelocityLimit = false;
        bool enableBoundaryStop = false;
        bool enableFieldNormalization = false;
        
        float stepSize = 0.0f;
        int maxSteps = 0;
        float maxArcLength = 0.0f;
        float minVelocity = 0.0f;
    };
// Construction / Deconstruction
public:
    Integrator();
    virtual ~Integrator() = default;

// Methods
public:
    static vec2 RK4(const Volume* vol, const vec2& position, const float stepsize);
    static vec2 RK4Normalized(const Volume* vol, const vec2& position, const float stepsize);
    
    static void streamline(const Volume* vol, size3_t dims, const vec2& startPoint,
        int maxSteps, float stepSize, std::vector<vec2>& linePoints);
    
    static void streamline(const Volume* vol, size3_t dims, const vec2& startPoint,
        const Integrator::StreamlineOptions& opts, std::vector<vec2>& linePoints);
    
    // TODO: Build on the last assignment by either copying the integrator code
    // here and in the respective .cpp or include the header from that
    // assignment with #include <lablic/integrator.h> in the files
    // where it is needed.
    // You may want to consider adding a helper function that computes an entire streamline
    // if you have not done so for the last assignments already.
};

}  // namespace inviwo
