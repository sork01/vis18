// *neode.onsave* setgo g++ -I../ -I../../../inviwo/include -I../../../inviwo/modules -I../../../inviwo/modules/_generated -I../../../inviwo/ext -fsyntax-only topology.h
/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/ports/meshport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/eventproperty.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/properties/optionproperty.h>
#include <labtopo/utils/gradients.h>
#include <labtopo/labtopomoduledefine.h>
#include <labtopo/topobounds.h>

namespace inviwo
{

/** \docpage{org.inviwo.Topology, Vector Field Topology}
    ![](org.inviwo.Topology.png?classIdentifier=org.inviwo.Topology)

    Generate the topological skeleton of a vector field.

    ### Inports
      * __data__ The input here is 2-dimensional vector field (with vectors of
      two components thus two values within each voxel) but it is represented
      by a 3-dimensional volume.
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
    

    ### Outports
      * __outMesh__ The output mesh contains points and linesegments that make up
      the topological skeleton.
*/
class IVW_MODULE_LABTOPO_API Topology : public Processor
{
public:
    // All possible first order critical points in 2D vector fields.
    enum class TypeCP
    {
        Saddle = 0,
        AttractingNode = 1,
        RepellingNode = 2,
        AttractingFocus = 3,
        RepellingFocus = 4,
        Center = 5,
        Unknown = 6
    };

    // Colors according to the TypeCP enum.
    static const vec4 ColorsCP[7];

    // Construction / Deconstruction
  public:
    Topology();
    virtual ~Topology() = default;

    // Methods
  public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

  protected:
    // Our main computation function
    virtual void process() override;
    bool possibleZero(const Volume* vol, const TopoBounds& bounds);
    void findAllZeros(const Volume* vol, const vec2& subdivs, float threshold, std::vector<vec2>& zeros);
    bool findZero(const Volume* vol, const TopoBounds& bounds, float threshold, vec2& result);
    TypeCP classifyCriticalPoint(util::EigenResult eigens);

    //TODO: You may want to declare additional functions here, e.g., extractCriticalPoints.

    // Ports
  public:
    // Input data
    VolumeInport inData;

    // Output mesh
    MeshOutport outMesh;

// Properties
public:
    BoolProperty propClassifyCriticalPoints;
    BoolProperty propDrawSeparatrices;
    TemplateOptionProperty<int> propSeparatrixIntegrationDirection;
    FloatProperty propSeparatrixStepsize;
    IntProperty propSeparatrixStepLimit;
    FloatProperty propSeparatrixOffset;
};

}// namespace inviwo
