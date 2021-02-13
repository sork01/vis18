// *neode.onsave* setgo g++ -I../ -I../../../inviwo/include -I../../../inviwo/modules -I../../../inviwo/modules/_generated -I../../../inviwo/ext -fsyntax-only topology.cpp
/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labtopo/integrator.h>
#include <labtopo/interpolator.h>
#include <labtopo/topology.h>

#define SEPARATRIX_INTEGRATION_DIRECTION_AUTO 0
#define SEPARATRIX_INTEGRATION_DIRECTION_POSITIVE 1
#define SEPARATRIX_INTEGRATION_DIRECTION_NEGATIVE 2

namespace inviwo
{

static float costheta(vec2 a, vec2 b)
{
    return glm::dot(a, b) / (glm::length(a) * glm::length(b));
}

static float angle(vec2 a, vec2 b)
{
    return acos(costheta(a, b));
}

const vec4 Topology::ColorsCP[7] =
    {
        vec4(1, 1, 0, 1),  // Saddle             (Yellow)
        vec4(0, 0, 1, 1),  // AttractingNode     (Blue)
        vec4(1, 0, 0, 1),  // RepellingNode      (Red)
        vec4(0.5, 0, 1, 1),// AttractingFocus    (Purple)
        vec4(1, 0.5, 0, 1),// RepellingFocus     (Orange)
        vec4(0, 1, 0, 1),  // Center             (Green)
        vec4(1, 1, 1, 1)   // Unknown            (White)
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",  // Class identifier
    "Vector Field Topology",// Display name
    "KTH Lab",              // Category
    CodeState::Experimental,// Code state
    Tags::None,             // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const
{
    return processorInfo_;
}

Topology::Topology()
    : Processor(), outMesh("meshOut"), inData("inData")
    , propClassifyCriticalPoints("classifyCriticalPoints", "Classify Critical Points", false)
    , propDrawSeparatrices("drawSeparatrices", "Draw Separatrices", false)
    , propSeparatrixIntegrationDirection("separatrixIntegrationDirection", "Separatrix Integration Direction")
    , propSeparatrixOffset("separatrixOffset", "Separatrix Offset", 0.01f, 0.005f, 5.0f, 0.001f)
    , propSeparatrixStepsize("separatrixStepsize", "Separatrix Stepsize", 0.01f, 0.005f, 5.0f, 0.001f)
    , propSeparatrixStepLimit("separatrixStepLimit", "Separatrix Step Limit", 1000, 1, 10000, 1)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment (optional));
// propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);

    // TODO: Register additional properties
    // addProperty(propertyName);
    
    addProperty(propClassifyCriticalPoints);
    addProperty(propDrawSeparatrices);
    
    propSeparatrixIntegrationDirection.addOption("auto", "Automatic", SEPARATRIX_INTEGRATION_DIRECTION_AUTO);
    propSeparatrixIntegrationDirection.addOption("startPositive", "Positive then negative", SEPARATRIX_INTEGRATION_DIRECTION_POSITIVE);
    propSeparatrixIntegrationDirection.addOption("startNegative", "Negative then positive", SEPARATRIX_INTEGRATION_DIRECTION_NEGATIVE);
    addProperty(propSeparatrixIntegrationDirection);
    addProperty(propSeparatrixOffset);
    addProperty(propSeparatrixStepsize);
    addProperty(propSeparatrixStepLimit);
}

void Topology::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VolumeRAM* vr = vol->getRepresentation<VolumeRAM>();
    uvec3 dims = vr->getDimensions();


    auto base = vol->getBasis();
    vec2 cellSize(base[0][0] / dims[0], base[1][1] / dims[1]);
    
    // std::cerr << "basis is " << base[0][0] << "," << base[0][1] << "; " << base[1][0] << "," << base[1][1] << "\n";
    // std::cerr << "cellsize is " << cellSize.x << " by " << cellSize.y << "\n";
    
    // Initialize mesh, vertices and index buffers for the two streamlines and the
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer),
    // or use several index buffers with connectivity type adjacency.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.
    // You can use your previous integration code (copy it over or call it from <lablic/integrator.h>).

    // Looping through all values in the vector field.
    // for (int x = 0; x < dims[0]; ++x)
    // {
        // for (int y = 0; y < dims[1]; ++y)
        // {
            // dvec2 vectorValue = vr->getAsDVec2(uvec3(x, y, 0));
            // std::cerr << "voldata[" << x << "][" << y << "] = " << "(" << vectorValue[0] << "," << vectorValue[1] << ")\n";
        // }
    // }
    
    std::vector<vec2> zeros;
    findAllZeros(vol.get(), vec2(dims.x - 1, dims.y - 1), 0.00001f, zeros);
    
    std::cerr << "number of zeros: " << zeros.size() << "\n";
    
    std::cerr << "====================================================================\n";
    
    for (vec2 point : zeros)
    {
        mat2 jacobian = Interpolator::sampleJacobian(vol.get(), point);
        util::EigenResult eigens = util::eigenAnalysis(jacobian);
        
        TypeCP type = classifyCriticalPoint(eigens);
        
        vec4 color;
        
        if (propClassifyCriticalPoints.get())
        {
            color = ColorsCP[int(type)];
        }
        else
        {
            color = vec4(1,0,0,1);
        }
        
        indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(point.x / (dims.x - 1), point.y / (dims.y - 1), 0),
                            vec3(0), vec3(0), color});
        
        if (propDrawSeparatrices.get() && type == TypeCP::Saddle)
        {
            // Draw separatrices
            std::vector<vec2> linePoints;
            
            float offset = propSeparatrixOffset.get();
            
            Integrator::StreamlineOptions opts;
            opts.enableStepLimit = true;
            opts.maxSteps = propSeparatrixStepLimit.get();
            opts.enableBoundaryStop = true;
            opts.enableLowVelocityLimit = true;
            opts.minVelocity = 0.001f;
            
            if (propSeparatrixIntegrationDirection.get() == SEPARATRIX_INTEGRATION_DIRECTION_AUTO)
            {
                opts.stepSize = propSeparatrixStepsize.get() * (glm::dot(eigens.eigenvectors[0], Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[0])) > 0.0f ? 1 : -1);
            }
            else
            {
                opts.stepSize = propSeparatrixStepsize.get();
                
                if (propSeparatrixIntegrationDirection.get() == SEPARATRIX_INTEGRATION_DIRECTION_NEGATIVE)
                {
                    opts.stepSize = -opts.stepSize;
                }
            }
            
            // opts.enableFieldNormalization = true;
            
            std::cerr << "our point is " << point << "\n";
            std::cerr << "the first eigenvector is " << eigens.eigenvectors[0] << "\n";
            // std::cerr << "the first eigenvector (other way) is " << -eigens.eigenvectors[0] << "\n";
            std::cerr << "the second eigenvector is " << eigens.eigenvectors[1] << "\n";
            // std::cerr << "the second eigenvector (other way) is " << -eigens.eigenvectors[1] << "\n";
            
            std::cerr << "moving a little bit on the first eigenvector takes us to " << (point + offset*eigens.eigenvectors[0]) << "\n";
            std::cerr << "  sampling the field here gives us " << Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[0]) << "\n";
            std::cerr << "  the dot between the eigenvector and the field vector is " << glm::dot(eigens.eigenvectors[0], Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[0])) << "\n";
            std::cerr << "moving a little bit on the second eigenvector takes us to " << (point + offset*eigens.eigenvectors[1]) << "\n";
            std::cerr << "  sampling the field here gives us " << Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[1]) << "\n";
            std::cerr << "  the dot between the eigenvector and the field vector is " << glm::dot(eigens.eigenvectors[1], Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[1])) << "\n";
            
            // std::cerr << "moving the other way on the first eigenvector takes us to " << (point + offset*(-eigens.eigenvectors[0])) << "\n";
            // std::cerr << "moving a little bit on the second eigenvector takes us to " << (point + offset*eigens.eigenvectors[1]) << "\n";
            // std::cerr << "moving the other way on the second eigenvector takes us to " << (point + offset*(-eigens.eigenvectors[1])) << "\n";
            
            Integrator::streamline(vol.get(), dims, point + offset*eigens.eigenvectors[0], opts, linePoints);
            // std::cerr << "we gots " << linePoints.size() << " linePoints\n";
            
            for (int i = 0; i < (linePoints.size() - 1); ++i)
            {
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i].x / (dims.x - 1), linePoints[i].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
                                
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i+1].x / (dims.x - 1), linePoints[i+1].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
            }
            
            linePoints.clear();
            
            Integrator::streamline(vol.get(), dims, point + offset*(-eigens.eigenvectors[0]), opts, linePoints);
            // std::cerr << "we gots " << linePoints.size() << " linePoints\n";
            
            for (int i = 0; i < (linePoints.size() - 1); ++i)
            {
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i].x / (dims.x - 1), linePoints[i].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
                                
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i+1].x / (dims.x - 1), linePoints[i+1].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
            }
            
            if (propSeparatrixIntegrationDirection.get() == SEPARATRIX_INTEGRATION_DIRECTION_AUTO)
            {
                opts.stepSize = propSeparatrixStepsize.get() * (glm::dot(eigens.eigenvectors[1], Interpolator::sampleFromField(vol.get(), point + offset*eigens.eigenvectors[1])) > 0.0f ? 1 : -1);
            }
            else
            {
                opts.stepSize = -propSeparatrixStepsize.get();
                
                if (propSeparatrixIntegrationDirection.get() == SEPARATRIX_INTEGRATION_DIRECTION_NEGATIVE)
                {
                    opts.stepSize = -opts.stepSize;
                }
            }
            
            linePoints.clear();
            
            Integrator::streamline(vol.get(), dims, point + offset*eigens.eigenvectors[1], opts, linePoints);
            // std::cerr << "we gots " << linePoints.size() << " linePoints\n";
            
            for (int i = 0; i < (linePoints.size() - 1); ++i)
            {
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i].x / (dims.x - 1), linePoints[i].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
                                
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i+1].x / (dims.x - 1), linePoints[i+1].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
            }
            
            linePoints.clear();
            
            Integrator::streamline(vol.get(), dims, point + offset*(-eigens.eigenvectors[1]), opts, linePoints);
            // std::cerr << "we gots " << linePoints.size() << " linePoints\n";
            
            for (int i = 0; i < (linePoints.size() - 1); ++i)
            {
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i].x / (dims.x - 1), linePoints[i].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
                                
                indexBufferSeparatrices->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(linePoints[i+1].x / (dims.x - 1), linePoints[i+1].y / (dims.y - 1), 0),
                                    vec3(0), vec3(0), vec4(1,1,1,1)});
            }
        }
    }
    
    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}



Topology::TypeCP Topology::classifyCriticalPoint(util::EigenResult eigens)
{
    float R1 = eigens.eigenvaluesRe[0];
    float R2 = eigens.eigenvaluesRe[1];
    float I1 = eigens.eigenvaluesIm[0];
    float I2 = eigens.eigenvaluesIm[1];
    
    if (((R1 < 0 && R2 > 0) || (R2 < 0 && R1 > 0)) && I1 == 0 && I2 == 0)
    {
        return TypeCP::Saddle;
    }
    else if (R1 > 0 && R2 > 0 && I1 == 0 && I2 == 0)
    {
        return TypeCP::RepellingNode;
    }
    else if (R1 > 0 && R1 == R2 && I1 != 0 && I1 == -I2)
    {
        return TypeCP::RepellingFocus;
    }
    else if (R1 == 0 && R2 == 0 && I1 != 0 && I1 == -I2)
    {
        return TypeCP::Center;
    }
    else if (R1 < 0 && R1 == R2 && I1 != 0 && I1 == -I2)
    {
        return TypeCP::AttractingFocus;
    }
    else if (R1 < 0 && R2 < 0 && I1 == 0 && I2 == 0)
    {
        return TypeCP::AttractingNode;
    }
    else
    {
        return TypeCP::Unknown;
    }
}



bool Topology::possibleZero(const Volume* vol, const TopoBounds& bounds)
{
    std::vector<vec2> corners;
    bounds.getCornerPoints(corners);
    
    vec2 f00 = Interpolator::sampleFromField(vol, corners[0]);
    vec2 f10 = Interpolator::sampleFromField(vol, corners[1]);
    vec2 f01 = Interpolator::sampleFromField(vol, corners[2]);
    vec2 f11 = Interpolator::sampleFromField(vol, corners[3]);
    
    vec2 s00(f00[0] > 0 ? 1 : 0, f00[1] > 0 ? 1 : 0);
    vec2 s10(f10[0] > 0 ? 1 : 0, f10[1] > 0 ? 1 : 0);
    vec2 s01(f01[0] > 0 ? 1 : 0, f01[1] > 0 ? 1 : 0);
    vec2 s11(f11[0] > 0 ? 1 : 0, f11[1] > 0 ? 1 : 0);
    
    vec2 vz(1,1);
    
    return s00 + s10 == vz || s00 + s01 == vz  || s00 + s11 == vz || s10 + s01 == vz || s10 + s11 == vz || s01 + s11 == vz;
}



void Topology::findAllZeros(const Volume* vol, const vec2& subdivs, float threshold, std::vector<vec2>& zeros)
{
    const VolumeRAM* vr = vol->getRepresentation<VolumeRAM>();
    uvec3 dims = vr->getDimensions();
    
    vec2 subsize((dims.x - 1)  / subdivs.x, (dims.y - 1) / subdivs.y);
    std::vector<TopoBounds> subBounds;
    
    float xoffset = 0.0f;
    float yoffset = 0.0f;
    
    for (int y = 0; y < subdivs.y; ++y)
    {
        xoffset = 0.0f;
        
        for (int x = 0; x < subdivs.x; ++x)
        {
            subBounds.push_back(TopoBounds(vec2(xoffset, yoffset), subsize));
            xoffset += subsize.x;
        }
        
        yoffset += subsize.y;
    }
    
    for (TopoBounds bounds : subBounds)
    {
        vec2 result;
        
        if (findZero(vol, bounds, 0.0001f, result))
        {
            zeros.push_back(result);
        }
    }
}

bool Topology::findZero(const Volume* vol, const TopoBounds& bounds, float threshold, vec2& result)
{
    vec2 halfsize = bounds.size / 2.0f;
    
    if (possibleZero(vol, bounds))
    {
        if (bounds.size.x <= threshold)
        {
            result = bounds.getCenter();
            return true;
        }
        else
        {
            if (findZero(vol, TopoBounds(bounds.bottomLeft, halfsize), threshold, result))
            {
                return true;
            }
            else if (findZero(vol, TopoBounds(bounds.bottomLeft + vec2(halfsize.x, 0), halfsize), threshold, result))
            {
                return true;
            }
            else if (findZero(vol, TopoBounds(bounds.bottomLeft + vec2(0, halfsize.y), halfsize), threshold, result))
            {
                return true;
            }
            else if (findZero(vol, TopoBounds(bounds.bottomLeft + vec2(halfsize.x, halfsize.y), halfsize), threshold, result))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
}

}// namespace
