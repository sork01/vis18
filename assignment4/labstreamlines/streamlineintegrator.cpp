/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/streamlineintegrator.h>
#include <labstreamlines/integrator.h>
#include <inviwo/core/util/utilities.h>
#include <inviwo/core/interaction/events/mouseevent.h>

#include <random>
#include <cmath>

#define INTEGRATION_DIRECTION_RIGHT 0
#define INTEGRATION_DIRECTION_LEFT 1
#define INTEGRATION_DIRECTION_BOTH 2

#define MULTIMODE_RANDOM 0
#define MULTIMODE_UNIFORM 1
#define MULTIMODE_BY_MAGNITUDE 2

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , outMesh("meshOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , propSeedMode("seedMode", "Seeds")
    // TODO: Initialize additional properties
    // propertyName("propertyIdentifier", "Display Name of the Propery",
    // default value (optional), minimum value (optional), maximum value (optional), increment
    // (optional)); propertyIdentifier cannot have spaces
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move)
    , propMultiMode("multiMode", "Multipoint Mode")
    , propMultiNumPoints("multiNumPoints", "Number of points", 4, 4, 1000, 1)
    , propMultiGridX("multiGridX", "Grid X", 4, 2, 1000, 1)
    , propMultiGridY("multiGridY", "Grid Y", 4, 2, 1000, 1)
    , propDrawPoints("drawPoints", "Draw points", true)
    , propIntegrationDirection("integrationDirection", "Direction")
    , propIntegrationStepSize("integrationStepSize", "Step size", 0.15f, 0.001f, 5.0f, 0.001f)
    , propNormalizeField("normalizeField", "Normalize field", false)
    , propLimitIntegrationSteps("limitIntegrationSteps", "Limit steps", true)
    , propMaxIntegrationSteps("maxIntegrationSteps", "Maximum number of steps", 30, 1, 1000, 1)
    , propLimitArcLength("limitArcLength", "Limit arc length", false)
    , propMaxArcLength("maxArcLength", "Maximum arc length", 1.0f, 0.0f, 10.0f, 0.01f) // TODO: scale arc length either to fit image [0,1] or to fit data [min,max]
    , propStopAtDomainBoundary("stopAtDomainBoundary", "Stop at domain boundary", true)
    , propUseVelocityThreshold("useVelocityThreshold", "Use velocity threshold", true)
    , propVelocityThreshold("velocityThreshold", "Velocity threshold", 0.001f, 0.001f, 2.0f, 0.001f)
{
    // Register Ports
    addPort(inData);
    addPort(outMesh);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    propMultiMode.addOption("random", "Random", MULTIMODE_RANDOM);
    propMultiMode.addOption("uniform", "Uniform", MULTIMODE_UNIFORM);
    propMultiMode.addOption("byMagnitude", "By Magnitude", MULTIMODE_BY_MAGNITUDE);
    addProperty(propMultiMode);
    addProperty(propMultiNumPoints);
    addProperty(propMultiGridX);
    addProperty(propMultiGridY);

    addProperty(propDrawPoints);

    propIntegrationDirection.addOption("right", "Right", INTEGRATION_DIRECTION_RIGHT);
    propIntegrationDirection.addOption("left", "Left", INTEGRATION_DIRECTION_LEFT);
    propIntegrationDirection.addOption("both", "Both", INTEGRATION_DIRECTION_BOTH);
    addProperty(propIntegrationDirection);

    addProperty(propIntegrationStepSize);
    addProperty(propNormalizeField);
    addProperty(propLimitIntegrationSteps);
    addProperty(propMaxIntegrationSteps);
    addProperty(propLimitArcLength);
    addProperty(propMaxArcLength);
    addProperty(propStopAtDomainBoundary);
    addProperty(propUseVelocityThreshold);
    addProperty(propVelocityThreshold);

    // TODO: Register additional properties
    // addProperty(propertyName);

    // You can hide and show properties for a single seed and hide properties for multiple seeds (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
            // util::hide(...)
        } else {
            util::hide(propStartPoint, mouseMoveStart);
            // util::show(...)
        }
    });

}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    // Handle mouse interaction only if we
    // are in the mode with a single point
    if (propSeedMode.get() == 1) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();
    // Denormalize to volume dimensions
    mousePos.x *= dims.x - 1;
    mousePos.y *= dims.y - 1;
    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vr = vol->getRepresentation<VolumeRAM>();
    dims = vol->getDimensions();
    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMaxValue(vec2(dims.x - 1, dims.y - 1));

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0)
    {
        drawStreamLine(vr, mesh, vertices, propStartPoint.get());
    }
    else
    {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        
        if (propMultiMode.get() == MULTIMODE_RANDOM)
        {
            // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
            std::random_device randDev;
            std::mt19937 mersenneGenerator(randDev());
            std::uniform_real_distribution<> randX(0.0, dims.x - 1);
            std::uniform_real_distribution<> randY(0.0, dims.y - 1);
            
            for (int n = 0; n < propMultiNumPoints.get(); ++n)
            {
                drawStreamLine(vr, mesh, vertices, vec2(randX(mersenneGenerator), randY(mersenneGenerator)));
            }
        }
        else if (propMultiMode.get() == MULTIMODE_UNIFORM)
        {
            // auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
            
            int numPointsX = propMultiGridX.get() + 2;
            int numPointsY = propMultiGridY.get() + 2;
            // int numdrawn = 0;
            
            for (int x = 1; x < (numPointsX - 1); ++x)
            {
                for (int y = 1; y < (numPointsY - 1); ++y)
                {
                    float xf = (float(x)/float(numPointsX - 1)) * (dims.x - 1);
                    float yf = (float(y)/float(numPointsY - 1)) * (dims.y - 1);
                    
                    drawStreamLine(vr, mesh, vertices, vec2(xf, yf));
                    // ++numdrawn;
                    
                    // indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                    // vertices.push_back({vec3(xf / (dims.x - 1), yf / (dims.y - 1), 0),
                                        // vec3(0), vec3(0), vec4(0, 0, 0, 1)});
                }
            }
            
            // LogProcessorInfo("Uniform drew " << numdrawn << " streamlines");
        }
        else if (propMultiMode.get() == MULTIMODE_BY_MAGNITUDE)
        {
            int numRequestedPoints = propMultiNumPoints.get();
            
            std::vector<float> cellMagnitude;
            
            int gridX = 4;
            int gridY = 4;
            
            int numCells = gridX * gridY;
            
            float magnitudeSum = 0.0f;
            
            for (int c = 0; c < numCells; ++c)
            {
                float magn = findLargestMagnitudeUniform(vr, dims, gridX, gridY, c, 4, 4);
                magnitudeSum += magn;
                
                cellMagnitude.push_back(magn);
                
                LogProcessorInfo("Cell " << c << " has largest magnitude " << magn);
            }
            
            int numSeeded = 0;
            
            for (int c = 0; c < numCells; ++c)
            {
                cellMagnitude[c] /= magnitudeSum;
                // cellMagnitude[c] = 1;
                
                LogProcessorInfo("Seeding " << int(numRequestedPoints * cellMagnitude[c]) << " points in cell " << c);
                
                numSeeded += int(numRequestedPoints * cellMagnitude[c]);
                
                seedNStreamlinesRandom(vr, mesh, vertices, gridX, gridY, c, int(numRequestedPoints * cellMagnitude[c]));
            }
            
            LogProcessorInfo("Seeded a total of " << numSeeded << " points");
        }
    }
    
    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

/*
inline static double offsetInLinearInterpolation(double f0, double f1, double val)
{
    // linear interpolation s(x) of a function f(x) where x in [0,1] and f(0) = f0, f(1) = f1
    //   s(x) = (1-x)f0 + xf1
    //   s(x) = f0 - xf0 + xf1
    //   s(x) = f0 + x(f1 - f0)
    // solving s(x) = C
    //   f0 + x(f1 - f0) = C
    //   x(f1 - f0) = C - f0
    //   x = (C - f0)/(f1 - f0)
    
    return (val - f0)/(f1 - f0);
}
*/

// divide the data into <gridX> by <gridY> cells, scan <cell (0-based)> using <n> seeds
float StreamlineIntegrator::findLargestMagnitudeUniform(const VolumeRAM* vr, size3_t dims, int gridX, int gridY, int cell, int nX, int nY)
{
    // e.g gridX=gridY=4
    // cell=4
    // boundary width = 1 / gridX
    // boundary height = 1 / gridY
    // boundary left = (cell % gridX) / gridX
    // boundary bottom = (cell / gridY) / gridY
    
    float boundaryWidth = 1.0f / float(gridX);
    float boundaryHeight = 1.0f / float(gridY);
    float boundaryLeft = float(cell % gridX) / float(gridX);
    float boundaryRight = boundaryLeft + boundaryWidth;
    float boundaryBottom = float(int(float(cell) / float(gridY))) / float(gridY);
    float boundaryTop = boundaryBottom + boundaryHeight;
    
    int numPointsX = nX + 2;
    int numPointsY = nY + 2;
    float largestMagnitude = 0.0f;
    
    // for x in range(1, numPointsX - 1):
    for (int x = 1; x < (numPointsX - 1); ++x)
    {
        // for y in range(1, numPointsY - 1):
        for (int y = 1; y < (numPointsX - 1); ++y)
        {
            float xoffset = float(x)/float(numPointsX - 1);
            float yoffset = float(y)/float(numPointsY - 1);
            
            float px = boundaryLeft - xoffset*boundaryLeft + xoffset*boundaryRight;
            float py = boundaryBottom - yoffset*boundaryBottom + yoffset*boundaryTop;
            
            float magn = glm::length(Integrator::sampleFromField(vr, dims, vec2(px * (dims.x - 1), py * (dims.y - 1))));
            
            if (magn > largestMagnitude)
            {
                largestMagnitude = magn;
            }
        }
    }
    
    return largestMagnitude;
}

// divide the data into <gridX> by <gridY> cells, scan <cell (0-based)> using <n> seeds
void StreamlineIntegrator::seedNStreamlinesRandom(
    const VolumeRAM* vr, std::shared_ptr<BasicMesh>& mesh,
    std::vector<BasicMesh::Vertex>& vertices, int gridX, int gridY, int cell, int n)
{
    float boundaryWidth = 1.0f / float(gridX);
    float boundaryHeight = 1.0f / float(gridY);
    float boundaryLeft = float(cell % gridX) / float(gridX);
    float boundaryRight = boundaryLeft + boundaryWidth;
    float boundaryBottom = float(int(float(cell) / float(gridY))) / float(gridY);
    float boundaryTop = boundaryBottom + boundaryHeight;
    
    // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device randDev;
    std::mt19937 mersenneGenerator(randDev());
    std::uniform_real_distribution<> randX(boundaryLeft * (dims.x - 1), boundaryRight * (dims.x - 1));
    std::uniform_real_distribution<> randY(boundaryBottom * (dims.y - 1), boundaryTop * (dims.y - 1));
    
    for (int i = 0; i < n; ++i)
    {
        drawStreamLine(vr, mesh, vertices, vec2(randX(mersenneGenerator), randY(mersenneGenerator)));
    }
}

void StreamlineIntegrator::drawStreamLine(const VolumeRAM* vr, std::shared_ptr<BasicMesh>& mesh,
                        std::vector<BasicMesh::Vertex>& vertices,
                        const vec2& startPoint)
{
    if (propIntegrationDirection.get() == INTEGRATION_DIRECTION_LEFT)
    {
        drawStreamLineReal(vr, mesh, vertices, startPoint, -1);
    }
    else if (propIntegrationDirection.get() == INTEGRATION_DIRECTION_RIGHT)
    {
        drawStreamLineReal(vr, mesh, vertices, startPoint, 1);
    }
    else if (propIntegrationDirection.get() == INTEGRATION_DIRECTION_BOTH)
    {
        drawStreamLineReal(vr, mesh, vertices, startPoint, -1);
        drawStreamLineReal(vr, mesh, vertices, startPoint, 1);
    }
}

void StreamlineIntegrator::drawStreamLineReal(const VolumeRAM* vr, std::shared_ptr<BasicMesh>& mesh,
                                          std::vector<BasicMesh::Vertex>& vertices,
                                          const vec2& startPoint, const int direction)
{
    // LogProcessorInfo("drawStreamLine(" << startPoint.x << ", " << startPoint.y << ")");
    vec2 (*integrator)(const VolumeRAM*, size3_t, const vec2&, const float) = Integrator::RK4;

    if (propNormalizeField.get())
    {
        // (c) Allow integration in the direction field (normalized vector field)
        integrator = Integrator::RK4Normalized;
    }

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);

    // Draw start point
    if (propDrawPoints.get())
    {
        indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
                            vec3(0), vec3(0), vec4(0, 0, 0, 1)});
    }

    indexBufferLine->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
                        vec3(0), vec3(0), vec4(0.4f, 0.4f, 0.4f, 1)});

    bool limitIntegrationSteps = propLimitIntegrationSteps.get();
    int maxIntegrationSteps = propMaxIntegrationSteps.get();

    bool limitArcLength = propLimitArcLength.get();
    float arcLength = 0.0f;
    float arcLengthMax = propMaxArcLength.get();

    bool drawPoints = propDrawPoints.get();
    bool stopAtBoundary = propStopAtDomainBoundary.get();

    bool useVelocityThreshold = propUseVelocityThreshold.get();
    float velocityThreshold = propVelocityThreshold.get();

    vec2 oldPosition = startPoint;
    vec2 newPosition;

    // (b) Allow different step sizes
    float stepSize = propIntegrationStepSize.get();
    int infoStepsTaken = 0;
    
    /*
    if (propIntegrationDirection.get() != INTEGRATION_DIRECTION_RIGHT)
    {
        // (a) Allow integration in forward and backward direction
        stepSize = -stepSize;
    }
    */
    
    stepSize *= direction;

    for (int i = 0;; ++i)
    {
        if (limitIntegrationSteps && i >= maxIntegrationSteps)
        {
            // (d) Stop the integration after a certain number of steps
            break;
        }

        // static vec2 Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, const float stepsize)
        newPosition = integrator(vr, dims, oldPosition, stepSize);

        indexBufferLine->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(newPosition.x / (dims.x - 1), newPosition.y / (dims.y - 1), 0), vec3(0), vec3(0), vec4(0.4f, 0.4f, 0.4f, 1)});

        if (drawPoints)
        {
            indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
            vertices.push_back({vec3(newPosition.x / (dims.x - 1), newPosition.y / (dims.y - 1), 0), vec3(0), vec3(0), vec4(0, 0, 1, 1)});
        }

        ++infoStepsTaken;

        arcLength += glm::length(newPosition - oldPosition);

        if (limitArcLength && arcLength >= arcLengthMax)
        {
            // (e) Stop the integration after a certain arc length of the stream line
            // LogProcessorInfo("Arc length limit reached");
            break;
        }

        if (stopAtBoundary && (newPosition.x < 0 || newPosition.x > (dims.x - 1)
            || newPosition.y < 0 || newPosition.y > (dims.y - 1)))
        {
            // (f) Stop the integration at the boundary of the domain
            // TODO: clip the last line exactly at the boundary?
            // LogProcessorInfo("Boundary stop (position: " << newPosition.x << ", " << newPosition.y);
            break;
        }

        if (useVelocityThreshold && glm::length(Integrator::sampleFromField(vr, dims, newPosition)) <= velocityThreshold)
        {
            // (g) Stop the integration at zeros of the vector field
            // (h) Stop the integration when the velocity becomes too slow
            // LogProcessorInfo("Velocity threshold stop (position: " << newPosition.x << ", " << newPosition.y);
            break;
        }

        oldPosition = newPosition;
    }

    // LogProcessorInfo("Steps taken: " << infoStepsTaken);
}

}  // namespace inviwo
