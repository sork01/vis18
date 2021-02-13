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

#include <lablic/lablicmoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <functional>

namespace inviwo {

class IVW_MODULE_LABLIC_API Integrator {
// Friends
// Types
public:
// Construction / Deconstruction
public:
    Integrator();
    virtual ~Integrator() = default;

// Methods
public:

    // TODO: Build on the last assignment by copying the integrator code
    // here and in the respective .cpp. Mark a small change in the vector field sampling
    // (change the parameter VolumeRAM* to Volume*, i.e., vr to vol.get()).
    // You may want to consider adding a helper function that computes an entire streamline
    // if you have not done so for the last assignment already.

};

}  // namespace inviwo
