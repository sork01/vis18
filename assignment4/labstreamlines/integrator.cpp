/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo
{

Integrator::Integrator()
{
}

vec2 Integrator::sampleFromField(const VolumeRAM* vr, size3_t dims, const vec2& position)
{
    // Sampled outside the domain!
    if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 || position[1] > dims[1] - 1)
    {
        return vec2(0, 0);
    }

    int x0 = int(position[0]);
    int y0 = int(position[1]);

    // Leads to accessing only inside the volume
    // Coefficients computation takes care of using the correct values
    if (x0 == dims[0] - 1)
    {
        x0--;
    }
    if (y0 == dims[1] - 1)
    {
        y0--;
    }

    auto f00 = vr->getAsDVec2(size3_t(x0, y0, 0));
    auto f10 = vr->getAsDVec2(size3_t(x0 + 1, y0, 0));
    auto f01 = vr->getAsDVec2(size3_t(x0, y0 + 1, 0));
    auto f11 = vr->getAsDVec2(size3_t(x0 + 1, y0 + 1, 0));

    float x = position[0] - x0;
    float y = position[1] - y0;

    vec2 f;

    for (int i = 0; i < 2; i++)
    {
        f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) + f11[i] * x * y;
    }

    return f;
}

vec2 Integrator::Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, const float stepsize)
{
  return position + stepsize * sampleFromField(vr, dims, position);
}

vec2 Integrator::RK4(const VolumeRAM* vr, size3_t dims, const vec2& position, const float stepsize)
{
    vec2 v1 = sampleFromField(vr, dims, position);
    vec2 v2 = sampleFromField(vr, dims, position + (stepsize/2.0f)*v1);
    vec2 v3 = sampleFromField(vr, dims, position + (stepsize/2.0f)*v2);
    vec2 v4 = sampleFromField(vr, dims, position + stepsize*v3);

    return position + stepsize * ((v1/6.0f) + (v2/3.0f) + (v3/3.0f) + (v4/6.0f));
}

vec2 Integrator::RK4Normalized(const VolumeRAM* vr, size3_t dims, const vec2& position, const float stepsize)
{
    vec2 v1 = sampleFromField(vr, dims, position);
    vec2 v2 = sampleFromField(vr, dims, position + (stepsize/2.0f)*v1);
    vec2 v3 = sampleFromField(vr, dims, position + (stepsize/2.0f)*v2);
    vec2 v4 = sampleFromField(vr, dims, position + stepsize*v3);
    
    if (glm::length(v1) < 0.0001f || glm::length(v2) < 0.0001f || glm::length(v3) < 0.0001f || glm::length(v4) < 0.0001f)
    {
        return position;
    }
    
    v1 = glm::normalize(v1);
    v2 = glm::normalize(v2);
    v3 = glm::normalize(v3);
    v4 = glm::normalize(v4);
    
    return position + stepsize * ((v1/6.0f) + (v2/3.0f) + (v3/3.0f) + (v4/6.0f));
}
// TODO: Implement a single integration step here

//vec2 Integrator::Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, ...)
//{
// Access the vector field with sampleFromField(vr, dims, ...)
//}

//vec2 Integrator::RK4(const VolumeRAM* vr, size3_t dims, const vec2& position, ...)
//{
//
//}

} // namespace
