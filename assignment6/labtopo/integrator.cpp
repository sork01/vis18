// *neode.onsave* setgo g++ -I../ -I../../../inviwo/include -I../../../inviwo/modules -I../../../inviwo/modules/_generated -I../../../inviwo/ext -fsyntax-only integrator.cpp
/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labtopo/integrator.h>
#include <labtopo/interpolator.h>

namespace inviwo {

Integrator::Integrator() {}

vec2 Integrator::RK4(const Volume* vol, const vec2& position, const float stepsize)
{
    vec2 v1 = Interpolator::sampleFromField(vol, position);
    vec2 v2 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v1);
    vec2 v3 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v2);
    vec2 v4 = Interpolator::sampleFromField(vol, position + stepsize*v3);

    return position + stepsize * ((v1/6.0f) + (v2/3.0f) + (v3/3.0f) + (v4/6.0f));
}

vec2 Integrator::RK4Normalized(const Volume* vol, const vec2& position, const float stepsize)
{
    vec2 v1 = Interpolator::sampleFromField(vol, position);
    vec2 v2 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v1);
    vec2 v3 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v2);
    vec2 v4 = Interpolator::sampleFromField(vol, position + stepsize*v3);
    
    if (glm::length(v1) < 0.0001f || glm::length(v2) < 0.0001f || glm::length(v3) < 0.0001f || glm::length(v4) < 0.0001f)
    {
        return vec2(0,0);
    }
    
    v1 = glm::normalize(v1);
    v2 = glm::normalize(v2);
    v3 = glm::normalize(v3);
    v4 = glm::normalize(v4);
    
    return position + stepsize * ((v1/6.0f) + (v2/3.0f) + (v3/3.0f) + (v4/6.0f));
}

void Integrator::streamline(const Volume* vol, size3_t dims, const vec2& startPoint,
    int maxSteps, float stepSize, std::vector<vec2>& linePoints)
{
    Integrator::StreamlineOptions opts;
    opts.enableStepLimit = true;
    opts.maxSteps = maxSteps;
    opts.stepSize = stepSize;
    
    streamline(vol, dims, startPoint, opts, linePoints);
}

void Integrator::streamline(const Volume* vol, size3_t dims, const vec2& startPoint,
    const Integrator::StreamlineOptions& opts, std::vector<vec2>& linePoints)
{
    float arcLength = 0.0f;
    
    vec2 oldPosition = startPoint;
    vec2 newPosition = oldPosition;
    
    linePoints.push_back(startPoint);
    
    for (int i = 0;; ++i)
    {
        if (opts.enableStepLimit && i >= opts.maxSteps)
        {
            break;
        }
        
        // std::cerr << "streamline with stepsize " << opts.stepSize << " at step " << i << " position " << oldPosition << "\n";
        
        if (opts.enableFieldNormalization)
        {
            newPosition = RK4Normalized(vol, oldPosition, opts.stepSize);
        }
        else
        {
            newPosition = RK4(vol, oldPosition, opts.stepSize);
        }
        
        // std::cerr << "streamline at step " << i << " moved to " << newPosition << "\n";
        
        linePoints.push_back(newPosition);
        
        if (opts.enableArcLengthLimit)
        {
            arcLength += glm::length(newPosition - oldPosition);
            
            if (arcLength >= opts.maxArcLength)
            {
                break;
            }
        }
        
        if (opts.enableBoundaryStop && (newPosition.x < 0 || newPosition.x > (dims.x - 1) || newPosition.y < 0 || newPosition.y > (dims.y - 1)))
        {
            break;
        }
        
        if (opts.enableLowVelocityLimit && glm::length(Interpolator::sampleFromField(vol, newPosition)) <= opts.minVelocity)
        {
            break;
        }
        
        oldPosition = newPosition;
    }
}

// TODO: Implementation of the functions defined in the header file integrator.h

}  // namespace inviwo
