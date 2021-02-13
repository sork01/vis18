/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, ...
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_
{
    "org.inviwo.MarchingSquares",      // Class identifier
    "Marching Squares",                // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,           // Code state
    Tags::None,                        // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const
{
    return processorInfo_;
}


MarchingSquares::MarchingSquares()
    :Processor()
    , inData("volumeIn")
    , meshOut("meshOut")
    , propDebug("debug", "Debug")
    , propShowGrid("showGrid", "Show Grid")
    , propDeciderType("deciderType", "Decider Type")
    , propGauss("gauss", "Use Gauss Smoothing")
    , propSigma("sigma", "Sigma Value")
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f),
        vec4(0.0f), vec4(1.0f), vec4(0.1f),
        InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f),
        vec4(0.0f), vec4(1.0f), vec4(0.1f),
        InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
{
    // Register ports
    addPort(inData);
    addPort(meshOut);
	
    // Register properties
    addProperty(propDebug);
    addProperty(propShowGrid);
    addProperty(propGridColor);
	
    addProperty(propDeciderType);
    propDeciderType.addOption("midpoint", "Mid Point", 0);
    propDeciderType.addOption("asymptotic", "Asymptotic", 1);

    addProperty(propGauss);
    addProperty(propSigma);
    addProperty(propMultiple);
    
    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clearPoints();
    propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propNumContours, propIsoTransferFunc);

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]()
    {
        if (propShowGrid.get())
        {
            util::show(propGridColor);
        }
        else
        {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]()
    {
        if (propMultiple.get() == 0)
        {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        }
        else
        {
            util::hide(propIsoValue);
            util::show(propIsoColor, propNumContours);
            
            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            // util::hide(propIsoValue, propIsoColor);
            // util::show(propNumContours, propIsoTransferFunc);
        }
    });

}

void MarchingSquares::process()
{
    if (!inData.hasData()) {
	    return;
    }
    
    bool debug = propDebug.get();
    bool gauss = propGauss.get();
    
    propSigma.setMinValue(0.05);
    propSigma.setMaxValue(2);
    
    // This results in a shared pointer to a volume
    auto vol = inData.getData(); // auto = std::shared_ptr<const Volume>

    // Extract the minimum and maximum value from the input data
    double minValue = vol->dataMap_.valueRange[0];
    double maxValue = vol->dataMap_.valueRange[1];

    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from
    
    // Retreive data in a form that we can access it
    const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
    const size3_t dims = vol->getDimensions();

    // Initialize mesh and vertices
    auto mesh = std::make_shared<BasicMesh>(); // auto = std::shared_ptr<BasicMesh>
    std::vector<BasicMesh::Vertex> vertices;

    // Values within the input data are accessed by the function below
    // It's input is the VolumeRAM from above, the dimensions of the volume
    // and the indeces i and j of the position to be accessed where

    Volume volSmoothed;
    VolumeRAM* vrSmoothed;
    
    if (propGauss.get())
    {
        int i, j;
        float kernel[5][5];
        float sigma = propSigma.get();
        // Make a new Volume for the smoothed isolines
        volSmoothed = Volume(vol->getDimensions(), vol->getDataFormat());
        vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
        
	    // Making the 5x5 Gauss kernel
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                kernel[i][j] = 1.0/(2.0*M_PI*pow(sigma,2)) * exp(-(pow(i-2,2) + pow(j-2,2))/(2.0*pow(sigma,2)));
            }
        }
        // Calculate the sum of the kernel for normalizing
        float sum = 0;
        for (i = 0; i < 5; i++)
        {
            for (j = 0; j < 5; j++) {
                sum += kernel[i][j];
            }
        }
        // Normalizing the kernel
        for (i = 0; i < 5; i++)
        {
            for (j = 0; j < 5; j++) {
                kernel[i][j] = kernel[i][j]/sum;
            }
        }
        
	    // Applying kernel to input data.
        sum = 0;
        float res = 0;
        for (int x = 0; x < dims.x; x++) 
        {
            for (int y = 0; y < dims.y; y++) 
            {
                for (int kx = 0; kx < 5; kx++) 
                {
                    for (int ky = 0; ky < 5; ky++) 
                    {
                        res += kernel[kx][ky] * getInputValue(vr, dims, std::min((int)dims.x -1, std::max(x + kx - 2, 0)),std::min((int)dims.y -1, std::max(y + ky - 2, 0)));
                    }
                }
                vrSmoothed->setFromDouble(vec3(x,y,0), res);
                res = 0;
            }
        }
        
        // Update the data value range
        minValue = getInputValue(vrSmoothed, dims, 0, 0);
        maxValue = getInputValue(vrSmoothed, dims, 0, 0);
        
        for (i = 0; i < dims.x; ++i)
        {
            for (j = 0; j < dims.y; ++j)
            {
                double val = getInputValue(vrSmoothed, dims, i, j);
                
                if (val < minValue)
                {
                    minValue = val;
                }
                
                if (val > maxValue)
                {
                    maxValue = val;
                }
            }
        }
    }
    
    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue << ".");
    
    // float valueat00 = getInputValue(vr, dims, 0, 0);
    // LogProcessorInfo("Value at (0,0) is: " << valueat00);
    // You can assume that dims.z = 1 and do not need to consider others cases

    // TODO (Bonus) Gaussian filter
    // Our input is const, but you need to compute smoothed data and write it somewhere
    // Create an editable volume like this:
    // Volume volSmoothed(vol->getDimensions(), vol->getDataFormat());
    // auto vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
    // Values can be set with
    // vrSmoothed->setFromDouble(vec3(i,j,0), value);
    // getting values works with an editable volume as well
    // getInputValue(vrSmoothed, dims, 0, 0);

    // Grid

    // Properties are accessed with propertyName.get() 
    if (propShowGrid.get())
    {
        // TODO: Add grid lines of the given color 

        // The function drawLineSegments creates two vertices at the specified positions, 
        // that are placed into the Vertex vector defining our mesh. 
        // An index buffer specifies which of those vertices should be grouped into to make up lines/trianges/quads.
        // Here two vertices make up a line segment.
        auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None); // auto = IndexBufferRAM*

        // Draw a line segment from v1 to v2 with a color, the coordinates in the final 
        // image range from 0 to 1 for both x and y
        int x_dim = dims.x - 1;
        int y_dim = dims.y - 1;

        // Horizontal
        for (int i = 0; i <= x_dim; i++) {
            float point = i / (x_dim * 1.0);
            vec2 h1 = vec2(0, point);  // x-coordinate
            vec2 h2 = vec2(1, point);  // y-coordinate
            drawLineSegment(h1, h2, propGridColor.get(), indexBufferGrid,
                            vertices);  // draw a line from x to y
        }
        // Vertical
        for (int j = 0; j <= y_dim; j++) {
            float point = j / (y_dim * 1.0);
            vec2 h1 = vec2(point, 0);  // x-coordinate
            vec2 h2 = vec2(point, 1);  // y-coordinate
            drawLineSegment(h1, h2, propGridColor.get(), indexBufferGrid,
                            vertices);  // draw a line from x to y
        }
    }

    // Iso contours

    if (propMultiple.get() == 0)
    {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue) 
        // and color it with the specified color (propIsoColor)
        IndexBufferRAM* buf = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        if (propGauss.get()) {
            drawIsoLine(vrSmoothed, buf, vertices, propIsoValue.get(), propIsoColor.get());
            
        }
        else {
        drawIsoLine(vr, buf, vertices, propIsoValue.get(), propIsoColor.get());
        }
    }
    else
    {
        // We want numContours=1 to draw the middle isovalue and numContours=2 to draw 0.33 and 0.66
        // so for numContours=1 we use 3 points by splitting the range into (numpoints - 1) equally
        // sized ranges and then we draw isolines at the values of all points except the first and last
        
        IndexBufferRAM* buf;
        int numPoints = propNumContours.get() + 2;
        double step = (maxValue - minValue) / double(numPoints - 1);
        
        for (int i = 1; i < (numPoints - 1); ++i)
        {
            if (debug) { LogProcessorInfo("Drawing isovalue " << float(minValue + i * step) << " ("
                << offsetInLinearInterpolation(minValue, maxValue, float(minValue + i * step))
                << " in range)"); }
            
            buf = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
            if (propGauss.get()) 
            {
                drawIsoLine(vrSmoothed, buf, vertices, minValue + i * step, propIsoColor.get());
            }
            else 
            {
                drawIsoLine(vr, buf, vertices, minValue + i * step, propIsoColor.get());
            }
        }
        
        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshOut.setData(mesh);

}

double MarchingSquares::getInputValue(const VolumeRAM* data, const size3_t dims, 
    const size_t i, const size_t j) {
    // Check if the indices are withing the dimensions of the volume
    if (i < dims.x && j < dims.y) {
        return data->getAsDouble(size3_t(i, j, 0));
    } else {
        LogProcessorError(
            "Attempting to access data outside the boundaries of the volume, value is set to 0");
        return 0;
    }
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}


void MarchingSquares::drawIsoLine(const VolumeRAM* vr, IndexBufferRAM* indexBuffer,
                                  std::vector<BasicMesh::Vertex>& vertices,
                                  const double isoValue, const vec4& color)
{
    bool debug = propDebug.get();
    
    // step 1: mark each grid vertex (i,j) with + or - based on whether value(i,j) is larger or smaller
            // than the requested isovalue c
    
    const size3_t dims = vr->getDimensions();
    
    std::vector<std::vector<int>> vertexSigns = std::vector<std::vector<int>>();
    
    for (unsigned int i = 0; i < dims.x; ++i)
    {
        vertexSigns.push_back(std::vector<int>(dims.y));
    }
    
    for (unsigned int x = 0; x < dims.x; ++x)
    {
        for (unsigned int y = 0; y < dims.y; ++y)
        {
            if (debug) { LogProcessorInfo("Data[" << x << "," << y << "] = " << vr->getAsDouble(size3_t(x, y, 0))); }
            
            if (vr->getAsDouble(size3_t(x, y, 0)) >= isoValue)
            {
                vertexSigns[x][y] = 1;
            }
            else
            {
                vertexSigns[x][y] = 0;
            }
            
            if (debug) { LogProcessorInfo("vertexSigns(" << x << ", " << y << ") = " << vertexSigns[x][y]); }
        }
    }
    
    // step 2: find the intersected grid cells, i.e cells where not all four vertices have the same sign
    
    // there are (dims.x - 1) by (dims.y - 1) cells
    // once a cell (i,j) has been picked, the vertices are 
    //  v0 = (i,j)
    //  v1 = (i+1,j)
    //  v2 = (i,j+1)
    //  v3 = (i+1,j+1)
    
    for (unsigned int i = 0; i < (dims.x - 1); ++i)
    {
        for (unsigned int j = 0; j < (dims.y - 1); ++j)
        {
            // v3 ----- v2
             // |        |
             // |        |
             // |        |
            // v0 ----- v1
            
            // cell corner vertices
            vec2 v[4];
            
            v[0] = vec2(i, j);
            v[1] = vec2(i+1, j);
            v[2] = vec2(i+1, j+1);
            v[3] = vec2(i, j+1);
            
            int signs = vertexSigns[v[0].x][v[0].y] + vertexSigns[v[1].x][v[1].y]
                        + vertexSigns[v[2].x][v[2].y] + vertexSigns[v[3].x][v[3].y];
            
            // If all signs are the same this cell is not intersected
            if (signs == 0 || signs == 4)
            {
                continue;
            }
            
            if (debug) { LogProcessorInfo("Grid cell " << i << ", " << j << " is intersected."); }
            if (debug) { LogProcessorInfo("Grid cell vertices: (" << v[0].x << "," << v[0].y << "), (" << v[1].x << "," << v[1].y << "), (" << v[2].x << "," << v[2].y << "), (" << v[3].x << "," << v[3].y << ")"); }
            
            //         e2
            //    v3 ----- v2
             //    |        |
             // e3 |        | e1
             //    |        |
            //    v0 ----- v1
            //         e0
            
            int edgesum[4] = {
                vertexSigns[v[0].x][v[0].y] + vertexSigns[v[1].x][v[1].y],
                vertexSigns[v[1].x][v[1].y] + vertexSigns[v[2].x][v[2].y],
                vertexSigns[v[2].x][v[2].y] + vertexSigns[v[3].x][v[3].y],
                vertexSigns[v[3].x][v[3].y] + vertexSigns[v[0].x][v[0].y]};
            
            if (debug) { LogProcessorInfo("Edge sums: " << edgesum[0] << "," << edgesum[1] << "," << edgesum[2] << "," << edgesum[3]); }
            
            // data values at cell corners
            double val[4] = {
                vr->getAsDouble(size3_t(v[0].x, v[0].y, 0)),
                vr->getAsDouble(size3_t(v[1].x, v[1].y, 0)),
                vr->getAsDouble(size3_t(v[2].x, v[2].y, 0)),
                vr->getAsDouble(size3_t(v[3].x, v[3].y, 0))};
            
            if (debug) { LogProcessorInfo("Data values: " << val[0] << "," << val[1] << "," << val[2] << "," << val[3]); }
            
    // step 3: for each intersected grid cell, determine the type of intersection from the three core
            // cases: (1) a single vertex/corner has sign different to the three others,
                   // (2) two adjacent (sharing an edge) vertices have differing signs to the two others
                   // (3) two non-adjacent (diagonal) vertices have differing signs to the two others
            
            if (signs == 1 || signs == 3)
            {
                // case 1) starting from the differing vertex, find the point of intersection on the x-axis-
                        // aligned edge adjacent to the vertex by using linear interpolation, then do the
                        // same for the point of intersection on the y-axis-aligned edge. draw a line 
                        // connecting the two points of intersection.
                
                //  1 ------ 0        and any rotation
                 // |        |
                 // |        |
                 // |        |
                //  0 ------ 0
                
                //      or
                
                //  0 ------ 1        and any rotation
                 // |        |
                 // |        |
                 // |        |
                //  1 ------ 1
                
                //         e2
                //    v3 ----- v2
                 //    |        |
                 // e3 |        | e1
                 //    |        |
                //    v0 ----- v1
                //         e0
                
                // in both cases there are exactly two edges with sum=1
                
                for (int n = 0; n < 4; ++n)
                {
                    if (edgesum[n] + edgesum[(n + 1) % 4] == 2) // e0 and e1, e1 and e2, ...
                    {
                        // the differing vertex is v[(n+1) % 4]
                        
                        if (debug) { LogProcessorInfo("Found the differing vertex v[" << ((n + 1) % 4) << "] between edges e[" << n << "] and e[" << ((n + 1) % 4) << "]"); }
                        
                        // the two edges to interpolate on are (v[n], v[(n+1) % 4]) and (v[(n+1) % 4], v[(n+2) % 4])
                        
                        float p1offset = offsetInLinearInterpolation(val[n], val[(n+1) % 4], isoValue);
                        float p2offset = offsetInLinearInterpolation(val[(n+1) % 4], val[(n+2) % 4], isoValue);
                        
                        if (debug) { LogProcessorInfo("Offsets: " << p1offset << ", " << p2offset); }
                        
                        vec2 p1 = v[n] + p1offset * (v[(n+1) % 4] - v[n]);
                        vec2 p2 = v[(n+1) % 4] + p2offset * (v[(n+2) % 4] - v[(n+1) % 4]);
                        
                        if (debug) { LogProcessorInfo("Point 1: " << p1.x << "," << p1.y); }
                        if (debug) { LogProcessorInfo("Point 2: " << p2.x << "," << p2.y); }
                        
                        p1.x /= float(dims.x - 1);
                        p1.y /= float(dims.y - 1);
                        
                        p2.x /= float(dims.x - 1);
                        p2.y /= float(dims.y - 1);
                        
                        drawLineSegment(p1, p2, color, indexBuffer, vertices);
                    }
                }
            }
            else if (signs == 2 && (edgesum[0] == 2 || edgesum[1] == 2 || edgesum[2] == 2 || edgesum[3] == 2))
            {
                // case 2) by linear interpolation find the two points of intersection on both of the
                        // edges each adjacent to a pair of differing vertices. draw a line connecting the
                        // points.
                
                //  1 ------ 0        and any rotation
                 // |        |
                 // |        |
                 // |        |
                //  1 ------ 0
                
                //         e2
                //    v3 ----- v2
                 //    |        |
                 // e3 |        | e1
                 //    |        |
                //    v0 ----- v1
                //         e0
                
                // exactly one edge has sum=2
                
                for (int n = 0; n < 4; ++n)
                {
                    if (edgesum[n] == 2)
                    {
                        // the two intersected edges are (v[n], v[(n-1) % 4]) and (v[(n+1) % 4], v[(n+2) % 4])
                        
                        if (debug) { LogProcessorInfo("Found the edge e[" << n << "]"); }
                        
                        float p1offset = offsetInLinearInterpolation(val[n], val[(n+4-1) % 4], isoValue);
                        float p2offset = offsetInLinearInterpolation(val[(n+1) % 4], val[(n+2) % 4], isoValue);
                        
                        if (debug) { LogProcessorInfo("Offsets: " << p1offset << ", " << p2offset); }
                        
                        vec2 p1 = v[n] + p1offset * (v[(n+4-1) % 4] - v[n]);
                        vec2 p2 = v[(n+1) % 4] + p2offset * (v[(n+2) % 4] - v[(n+1) % 4]);
                        
                        if (debug) { LogProcessorInfo("Point 1: " << p1.x << "," << p1.y); }
                        if (debug) { LogProcessorInfo("Point 2: " << p2.x << "," << p2.y); }
                        
                        p1.x /= float(dims.x - 1);
                        p1.y /= float(dims.y - 1);
                        
                        p2.x /= float(dims.x - 1);
                        p2.y /= float(dims.y - 1);
                        
                        drawLineSegment(p1, p2, color, indexBuffer, vertices);
                    }
                }
            }
            else
            {
                // case 3) by linear interpolation on each of the four edges find the four points of
                        // intersection, then use an ambiguity resolution strategy to decide which pairs
                        // of points should be connected.
                
                //  1 ------ 0        and any rotation
                 // |        |
                 // |        |
                 // |        |
                //  0 ------ 1
                
                //         e2
                //    v3 ----- v2
                 //    |        |
                 // e3 |        | e1
                 //    |        |
                //    v0 ----- v1
                //         e0
                
                // all four edges are intersected so we don't need to search for edges (yay)
                
                float p01offset = offsetInLinearInterpolation(val[0], val[1], isoValue);
                float p12offset = offsetInLinearInterpolation(val[1], val[2], isoValue);
                float p23offset = offsetInLinearInterpolation(val[2], val[3], isoValue);
                float p30offset = offsetInLinearInterpolation(val[3], val[0], isoValue);
                
                // vec2 p01 = v[0] + p01offset * (v[1] - v[0]);
                // vec2 p12 = v[1] + p12offset * (v[2] - v[1]);
                // vec2 p23 = v[2] + p23offset * (v[3] - v[2]);
                // vec2 p30 = v[3] + p30offset * (v[0] - v[3]);
                
                vec2 p[4] = {
                    v[0] + p01offset * (v[1] - v[0]), // p01
                    v[1] + p12offset * (v[2] - v[1]), // p12
                    v[2] + p23offset * (v[3] - v[2]), // p23
                    v[3] + p30offset * (v[0] - v[3])}; // p30
                
                p[0].x /= float(dims.x - 1);
                p[0].y /= float(dims.y - 1);
                p[1].x /= float(dims.x - 1);
                p[1].y /= float(dims.y - 1);
                p[2].x /= float(dims.x - 1);
                p[2].y /= float(dims.y - 1);
                p[3].x /= float(dims.x - 1);
                p[3].y /= float(dims.y - 1);
                
                int whichCase = 0;
                
                // midpoint decider (bad strategy)
                if (propDeciderType.get() == 0)
                {
                    if (.25 * (val[0] + val[1] + val[2] + val[3]) >= isoValue)
                    {
                        whichCase = (vertexSigns[v[0].x][v[0].y] == 0) ? 1 : 2;
                    }
                    else
                    {
                        whichCase = (vertexSigns[v[0].x][v[0].y] == 1) ? 1 : 2;
                    }
                    
                    if (whichCase == 1)
                    {
                        // upper case from slides
                        // connect p01 and p30, connect p12 and p23
                        drawLineSegment(p[0], p[3], color, indexBuffer, vertices);
                        drawLineSegment(p[1], p[2], color, indexBuffer, vertices);
                    }
                    else if (whichCase == 2)
                    {
                        // lower case from slides
                        // connect p01 and p12, connect p23 and p30
                        drawLineSegment(p[0], p[1], color, indexBuffer, vertices);
                        drawLineSegment(p[2], p[3], color, indexBuffer, vertices);
                    }
                }
                else if (propDeciderType.get() == 1) // asymptotic decider
                {
                    // sort the intersection points by either x- or y-coordinate
                    vec2 order[4];
                    
                    for (int i = 0; i < 4; ++i)
                    {
                        order[i] = p[i];
                    }
                    
                    for (int i = 0; i < 3; ++i)
                    {
                        for (int j = i + 1; j < 4; ++j)
                        {
                            if (order[j].x < order[i].x)
                            {
                                std::swap(order[i], order[j]);
                            }
                        }
                    }
                    
                    if (debug)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            LogProcessorInfo("Ordered point " << i << " has x=" << order[i].x);
                        }
                    }
                    
                    // from the sorted list, connect the first two and the last two
                    
                    drawLineSegment(order[0], order[1], color, indexBuffer, vertices);
                    drawLineSegment(order[2], order[3], color, indexBuffer, vertices);
                }
            }
        }
    }
}

inline double MarchingSquares::offsetInLinearInterpolation(double f0, double f1, double val)
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

} // namespace
