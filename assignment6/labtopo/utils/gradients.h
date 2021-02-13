/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <labtopo/labtopomoduledefine.h>
#include <labtopo/utils/tnt/jama_eig.h>

namespace inviwo {
namespace util {

struct EigenResult
{
    mat2 eigenvectors;
    vec2 eigenvaluesRe;
    vec2 eigenvaluesIm;
};

EigenResult eigenAnalysis(mat2 matrix);

}// namespace
}
