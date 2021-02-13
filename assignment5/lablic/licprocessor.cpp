/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/licprocessor.h>
#include <lablic/integrator.h>
#include <lablic/interpolator.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <inviwo/core/datastructures/transferfunction.h>
#include <ctime>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , propRender("renderEnable", "Render", false)
    , propFilterLength("filterLength", "Filter length", 20, 10, 200, 2)
    , propPreview("preview", "Preview mode", true)
    , propStepSize("stepSize", "Step size", 0.005f, 0.001f, 1.0f, 0.001f)
    , propColorize("colorize", "Colorize", false)
    , propStepSizeSuggested("stepSizeSuggested", "Suggested step size", 0.0f, 0.0f, 100.0f, 0.001f)
    , propCont("cont", "Contrast enhancement")

// TODO: Register additional properties

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    addProperty(propRender);
    addProperty(propPreview);
    addProperty(propColorize);
    addProperty(propFilterLength);
    addProperty(propStepSizeSuggested);
    addProperty(propStepSize);
    addProperty(propCont);

}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    vectorFieldDims_ = vol->getDimensions();
    auto vr = vol->getRepresentation<VolumeRAM>();

    // An accessible form of on image is retrieved analogous to a volume
    auto tex = noiseTexIn_.getData();
    texDims_ = tex->getDimensions();
    auto tr = tex->getRepresentation<ImageRAM>();
    
    float autoStepX = float(vectorFieldDims_.x) / float(texDims_.x);
    float autoStepY = float(vectorFieldDims_.y) / float(texDims_.y);
    
    float autoStep = autoStepX < autoStepY ? autoStepX: autoStepY;
    std::cerr << "autostep = " << autoStep;
    
    propStepSizeSuggested.setMinValue(autoStep);
    propStepSizeSuggested.setMaxValue(autoStep);
    // propStepSizeSuggested.setReadOnly(true);
    
    if (!propRender.get())
    {
        return;
    }
    
    // Prepare the output, it has the same dimensions and datatype as the output
    // and an editable version is retrieved analogous to a volume
    // auto outImage = tex->clone();
    auto outImage = std::make_shared<Image>();
    outImage->setDimensions(size2_t(texDims_.x, texDims_.y));
    
    auto outLayer = outImage->getColorLayer();
    outLayer->setDataFormat(DataVec4UInt8::get());
    
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();

    // To access the image at a floating point position, you can call
    //      Interpolator::sampleFromGrayscaleImage(tr, somePos)
    
    int filterLen = propFilterLength.get();
    
    LogProcessorInfo("Rendering started");
    float timeRenderStarted = clock() / CLOCKS_PER_SEC;
    
    int pixelStep = propPreview.get() ? 2 : 1;
    int prevPct = -1;
    
    for (auto y = 0; y < texDims_.y; y += pixelStep)
    {
        float timeElapsed = (clock() / CLOCKS_PER_SEC) - timeRenderStarted;
        int pct = int(100.0f*float(y)/float(texDims_.y));
        
        if (pct != prevPct)
        {
            std::cerr << pct << "% (" << int((timeElapsed/(float(y)/float(texDims_.y)))-timeElapsed) << " seconds remain)\n";
            prevPct = pct;
        }
        
        for (auto x = 0; x < texDims_.x; x += pixelStep)
        {
            vec2 fieldPos = squareCoordsConvert(
                vec2(x, y), vec2(0, 0), vec2(texDims_.x - 1, texDims_.y - 1),
                vec2(0, 0), vec2(vectorFieldDims_.x - 1, vectorFieldDims_.y - 1));
            
            std::vector<vec2> lineForward;
            std::vector<vec2> lineBack;
            
            // TODO: auto stepsize = datawidth / imagewidth ? data per pixel
            
            streamlineIntegrator(vol.get(), vectorFieldDims_, fieldPos, filterLen/2, propStepSize.get(), lineForward);
            streamlineIntegrator(vol.get(), vectorFieldDims_, fieldPos, filterLen/2, -propStepSize.get(), lineBack);
            
            double sum = 0.0;
            
            for (int i = 0; i < filterLen/2; ++i)
            {
                vec2 texPosFwd = squareCoordsConvert(
                    lineForward[i], vec2(0, 0), vec2(vectorFieldDims_.x - 1, vectorFieldDims_.y - 1),
                    vec2(0, 0), vec2(texDims_.x - 1, texDims_.y - 1));
                
                vec2 texPosBwd = squareCoordsConvert(
                    lineBack[i], vec2(0, 0), vec2(vectorFieldDims_.x - 1, vectorFieldDims_.y - 1),
                    vec2(0, 0), vec2(texDims_.x - 1, texDims_.y - 1));
                
                sum += Interpolator::sampleFromGrayscaleImage(tr, texPosFwd);
                sum += Interpolator::sampleFromGrayscaleImage(tr, texPosBwd);
            }
            
            int val = int(sum/double(filterLen));
            lr->setFromDVec4(size2_t(x, y), dvec4(val, val, val, 255));
            
            if (pixelStep == 2)
            {
                lr->setFromDVec4(size2_t(x+1, y), dvec4(val, val, val, 255));
                lr->setFromDVec4(size2_t(x, y+1), dvec4(val, val, val, 255));
                lr->setFromDVec4(size2_t(x+1, y+1), dvec4(val, val, val, 255));
            }
        }
    }
    
    float timeRenderEnded = clock() / CLOCKS_PER_SEC;
    
    LogProcessorInfo("Rendered in " << timeRenderEnded - timeRenderStarted << " seconds");
    
    float scale = 0.2f;
    if (propCont.get())
    {
        for (int j = 0; j < texDims_.y; j++) {
            for (int i = 0; i < texDims_.x; i++) {
                int val = int(lr->getAsDVec4(size2_t(i, j))[0]);
                
                if (val < 127)
                {
                    val = val * (1.0f - scale);
                }
                else 
                {
                    val = val * (1.0f + scale);
                    val = glm::clamp(val, 0, 255);
                }
                lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));

            }
        }
    }
    
    if (propColorize.get())
    {
        // Allocate some memory to store magnitude of each pixel
        float** mag = new float*[texDims_.x];
        for (auto i = 0; i < texDims_.x; i++){
            mag[i] = new float[texDims_.y];
        }
        
        //Calculate magnitude
        float minMag = 0, maxMag = 0, _mag = 0;
        for (auto j = 0; j < texDims_.y; j++){
            for (auto i = 0; i < texDims_.x; i++){
                
                vec2 fieldPos = squareCoordsConvert(
                    vec2(i, j), vec2(0, 0), vec2(texDims_.x - 1, texDims_.y - 1),
                    vec2(0, 0), vec2(vectorFieldDims_.x - 1, vectorFieldDims_.y - 1));
                
                // _mag = glm::length(Interpolator::sampleFromField(vol.get(), vec2(i, j)));
                _mag = glm::length(Interpolator::sampleFromField(vol.get(), fieldPos));
                mag[i][j] = _mag;
                if (_mag < minMag)
                    minMag = _mag;
                if (_mag > maxMag)
                    maxMag = _mag;
            }
        }
        
        TransferFunction tf = TransferFunction();
        tf.addPoint(0.0f, vec4(0.0f, 1.0f, 0.0f, 1.0f));
        tf.addPoint(0.5f, vec4(1.0f, 1.0f, 0.0f, 1.0f));
        tf.addPoint(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f));
        
        // double green[2] = {-105.0/255.0, 150.0/255.0};
        // double yellow[2] = {0.0/255.0, 255.0/255.0};
        // double red[2] = {100.0/255.0, 355.0/255.0};
        
        //Normalize magnitude (min-max scaling) and set intensity of color
        for (auto j = 0; j < texDims_.y; j++){
            for (auto i = 0; i < texDims_.x; i++){
                
                double val = lr->getAsDVec4(size2_t(i, j))[0];
                
                // double amt_green = 0.0;
                // double amt_yellow = 0.0;
                // double amt_red = 0.0;
                
                mag[i][j] = (mag[i][j] - minMag) / (maxMag - minMag);
                
                // double magVal = mag[i][j];
                
                // if (magVal > green[0] && magVal < green[1])
                // {
                    // amt_green = 1.0 - std::abs(magVal - (green[0] + (green[1] - green[0])/2.0))/((green[1] - green[0])/2.0);
                // }
                
                // if (magVal > yellow[0] && magVal < yellow[1])
                // {
                    // amt_yellow = 1.0 - std::abs(magVal - (yellow[0] + (yellow[1] - yellow[0])/2.0))/((yellow[1] - yellow[0])/2.0);
                // }
                
                // if (magVal > red[0] && magVal < red[1])
                // {
                    // amt_red = 1.0 - std::abs(magVal - (red[0] + (red[1] - red[0])/2.0))/((red[1] - red[0])/2.0);
                // }
                
                // double rgbMagn = std::sqrt(amt_green*amt_green + amt_yellow*amt_yellow + amt_red*amt_red);
                
                // if (std::abs(rgbMagn) < 0.0001)
                // {
                    // amt_green = amt_yellow = amt_red = 0;
                // }
                // else
                // {
                    // amt_green /= rgbMagn;
                    // amt_yellow /= rgbMagn;
                    // amt_red /= rgbMagn;
                // }
                
                vec4 colors = tf.sample(mag[i][j]);
                
                lr->setFromDVec4(size2_t(i, j), dvec4(val*colors[0], val*colors[1], val*colors[2], 255));
                
                // int Red, Green, Blue;
                // Red = int(val * mag[i][j]);
                // Green = 0;
                // Blue = int(abs(1 - val) * mag[i][j]);
                // lr->setFromDVec4(size2_t(i, j), dvec4(Red, Green, Blue, 255));
            }
        }
        
        for (auto i = 0; i < texDims_.x; i++){
            delete mag[i];
        }
        delete mag;
        mag = nullptr;
    }
    
    licOut_.setData(outImage);
}

// TODO: Move to integrator
// vec2 LICProcessor::RK4_(const Volume* vol, size3_t dims, const vec2& position, const float stepsize)
// {
    // vec2 v1 = Interpolator::sampleFromField(vol, position);
    // vec2 v2 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v1);
    // vec2 v3 = Interpolator::sampleFromField(vol, position + (stepsize/2.0f)*v2);
    // vec2 v4 = Interpolator::sampleFromField(vol, position + stepsize*v3);

    // return position + stepsize * ((v1/6.0f) + (v2/3.0f) + (v3/3.0f) + (v4/6.0f));
// }

vec2 LICProcessor::RK4(const Volume* vol, size3_t dims, const vec2& position, const float stepsize)
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

// TODO: Move to integrator?
void LICProcessor::streamlineIntegrator(const Volume* vol, size3_t dims, const vec2& startPoint,
    int stepLimit, const float stepSize, std::vector<vec2>& linePoints)
{
    // we will call method twice per point to integrate in "both directions"
    // integrate stream line up to a certain arc length from (x,y)
    // return vector of points on the streamline
    
    // float arcLength = 0.0f;
    
    vec2 oldPosition = startPoint;
    vec2 newPosition = oldPosition;
    
    linePoints.push_back(startPoint);
    
    for (int i = 0;; ++i)
    {
        if (i >= stepLimit)
        {
            break;
        }
        
        newPosition = RK4(vol, dims, oldPosition, stepSize);
        linePoints.push_back(newPosition);
        
        // arcLength += glm::length(newPosition - oldPosition);
        
        // if (arcLength >= arcLengthLimit)
        // {
            // break;
        // }
        
        /*
        // Boundary stop
        if (newPosition.x < 0 || newPosition.x > (dims.x - 1) || newPosition.y < 0 || newPosition.y > (dims.y - 1))
        {
            break;
        }
        
        // Velocity threshold stop
        if (glm::length(sampleFromField(vr, dims, newPosition)) <= velocityThreshold)
        {
            break;
        }
        */
        
        oldPosition = newPosition;
    }
}

vec2 LICProcessor::squareCoordsConvert(vec2 coord_in, vec2 min_in, vec2 max_in, vec2 min_out, vec2 max_out)
{
    vec2 offset_in = vec2();
    offset_in.x = float(coord_in.x)/float(max_in.x - min_in.x);
    offset_in.y = float(coord_in.y)/float(max_in.y - min_in.y);

    vec2 coord_out = vec2();

    // lerp
    coord_out.x = (1.0f - offset_in.x) * min_out.x + offset_in.x * max_out.x;
    coord_out.y = (1.0f - offset_in.y) * min_out.y + offset_in.y * max_out.y;

    return coord_out;
}

}  // namespace inviwo
