/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:36
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/noisetexturegenerator.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo NoiseTextureGenerator::processorInfo_{
    "org.inviwo.NoiseTextureGenerator",  // Class identifier
    "Noise Texture Generator",           // Display name
    "KTH Labs",                                  // Category
    CodeState::Experimental,             // Code state
    Tags::None,                          // Tags
};

const ProcessorInfo NoiseTextureGenerator::getProcessorInfo() const { return processorInfo_; }

NoiseTextureGenerator::NoiseTextureGenerator()
    : Processor()
    , texOut_("texOut")
    , texSize_("texSize", "Texture Size", vec2(512, 512), vec2(1, 1), vec2(2048, 2048), vec2(1, 1))
    , propBW("showBW", "Black & White")
    , propSeed("Seed", "Seed", 1, 1, 2147483647, 1)


    
    // TODO: Register additional properties

{
    // Register ports
    addPort(texOut_);

    // Register properties
    addProperty(texSize_);
    addProperty(propBW);
    addProperty(propSeed);
    


    // TODO: Register additional properties

}

void NoiseTextureGenerator::process() {
    // The output of the generation process is an Image
    auto outImage = std::make_shared<Image>();

    // In Inviwo, images consist of multiple layers, but only the
    // color layer is relevant for us
    auto outLayer = outImage->getColorLayer();

    outImage->setDimensions(size2_t(texSize_.get().x, texSize_.get().y));
    // With the data format DataVec4UInt8 values for RGB-alpha range between 0 and 255
    outLayer->setDataFormat(DataVec4UInt8::get());

    // Just like we did with the volume in other assignments we need to retrieve an editable
    // representation of the object we want to modify (here a layer)
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();


    std::srand(propSeed.get());

    for (int j = 0; j < texSize_.get().y; j++) {
        for (int i = 0; i < texSize_.get().x; i++) {
            
            int val = (std::rand()% 255) + 1;
            
            if (propBW) 
            {
                if (val < 128)
                {
                    val = 0;
                }
                else 
                {
                    val = 255;
                }
            }
            
            lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));
        }
    }
    
    texOut_.setData(outImage);
}

}  // namespace inviwo
