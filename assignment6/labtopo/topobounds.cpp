// *neode.onsave* setgo g++ -I/home/lur/Desktop/inviwo/ext -I/home/lur/Desktop/inviwo/modules/_generated -I/home/lur/Desktop/inviwo/include topobounds.cpp -o /tmp/topobounds && /tmp/topobounds
#include "topobounds.h"

TopoBounds::TopoBounds() : bottomLeft(0,0), size(0,0)
{
}

TopoBounds::TopoBounds(inviwo::vec2 _bottomLeft, inviwo::vec2 _size) : bottomLeft(_bottomLeft), size(_size)
{
}

TopoBounds::TopoBounds(const TopoBounds& other) : bottomLeft(other.bottomLeft), size(other.size)
{
}

void TopoBounds::getCornerPoints(std::vector<inviwo::vec2>& points) const
{
    points.emplace_back(inviwo::vec2(bottomLeft));
    points.emplace_back(inviwo::vec2(bottomLeft) + inviwo::vec2(size.x, 0));
    points.emplace_back(inviwo::vec2(bottomLeft) + inviwo::vec2(0, size.y));
    points.emplace_back(inviwo::vec2(bottomLeft) + inviwo::vec2(size.x, size.y));
}

inviwo::vec2 TopoBounds::getCenter() const
{
    return bottomLeft + inviwo::vec2(size.x/2.0f, size.y/2.0f);
}

std::string TopoBounds::toString() const
{
    std::ostringstream ss;
    ss << "TopoBounds(bottomLeft=" << bottomLeft << ", size=" << size << ")";
    return ss.str();
}

// int main(int argc, char** argv)
// {
    // TopoBounds b(inviwo::vec2(0,0), inviwo::vec2(10,20));
    
    // std::cout << b.toString() << "(center: " << b.getCenter() << ")\n";
// }
