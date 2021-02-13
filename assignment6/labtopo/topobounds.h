#pragma once
#include <inviwo/core/common/inviwo.h>

class TopoBounds
{
public:
    TopoBounds();
    TopoBounds(inviwo::vec2 _bottomLeft, inviwo::vec2 _size);
    TopoBounds(const TopoBounds& other);
    virtual ~TopoBounds() = default;
    void getCornerPoints(std::vector<inviwo::vec2>& points) const;
    inviwo::vec2 getCenter() const;
    std::string toString() const;
    
    inviwo::vec2 bottomLeft;
    inviwo::vec2 size;
};
