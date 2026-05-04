#ifndef ENGINE_NO_H
#define ENGINE_NO_H

#include <vector>

namespace Engine {

class No {
public:
    int id;
    double x, y;
    std::vector<int> gdlGlobais;

    No(int id, double x, double y, std::vector<int> gdls);
};

} // namespace Engine

#endif // ENGINE_NO_H
