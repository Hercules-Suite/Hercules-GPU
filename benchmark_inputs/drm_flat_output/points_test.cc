#include "../../quake/forward/drm_output_points.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>

template<class Function> void rejects(Function function) {
    bool failed = false;
    try { function(); } catch (const std::exception&) { failed = true; }
    assert(failed);
}

int main(int argc, char** argv) {
    std::istringstream csv("latitude,longitude,depth\r\n34,-118,0\n35,-117,2\n34,-118,0 # repeat\n");
    const auto points = drm_output::read_points(csv);
    assert(points.size() == 3);
    assert((points[0] == drm_output::Point{{34, -118, 0}}));
    assert((points[1] == drm_output::Point{{35, -117, 2}}));
    assert(points[0] == points[2]); // Input order and duplicate rows are intentional.
    std::istringstream space("# comment\nlatitude longitude depth\r\n34 -118 0\n");
    assert(drm_output::read_points(space).size() == 1);
    std::istringstream headerless("34 -118 0\n");
    assert(drm_output::read_points(headerless)[0] == points[0]);
    for (const char* bad : {"", "latitude longitude depth\n", "34 -118", "34 -118 0 1",
                            "nan -118 0", "34 inf 0", "91 0 0", "0 181 0", "0 0 -1",
                            "1a 0 0", "latitude longitude depth\nlatitude longitude depth"}) {
        rejects([bad]() { std::istringstream input(bad); drm_output::read_points(input); });
    }
    const double boundary[5] = {2, 3, 4, 6, 2};
    const auto box = drm_output::box_points(boundary, 1, 10, 20);
    assert(box.size() == 116); // Inner: 32; outer: 84.
    assert((box[0] == drm_output::Point{{12, 23, 0}}));
    assert((box[1] == drm_output::Point{{12, 23, 1}}));
    assert((box[2] == drm_output::Point{{12, 23, 2}}));
    assert((box[3] == drm_output::Point{{12, 24, 0}}));
    assert((box[32] == drm_output::Point{{11, 22, 0}}));
    const std::set<drm_output::Point> unique(box.begin(), box.end());
    assert(unique.size() == box.size()); // No duplicated edges/corners or shared surfaces.
    assert(unique.count(drm_output::Point{{13, 24, 0}}) == 0); // No top lid.
    assert(unique.count(drm_output::Point{{13, 24, 2}}) == 1); // Inner bottom.
    assert(unique.count(drm_output::Point{{13, 24, 3}}) == 1); // Outer bottom.
    for (size_t i = 0; i < box.size(); ++i) {
        const auto& p = box[i];
        if (i < 32) assert(p[0] == 12 || p[0] == 14 || p[1] == 23 || p[1] == 26 || p[2] == 2);
        else assert(p[0] == 11 || p[0] == 15 || p[1] == 22 || p[1] == 27 || p[2] == 3);
    }
    rejects([&]() { drm_output::box_points(boundary, 0, 0, 0); });
    rejects([&]() { drm_output::box_points(boundary, 0.7, 0, 0); });
    rejects([&]() { drm_output::box_points(boundary, 1e-10, 0, 0); });
    rejects([&]() { drm_output::box_points(boundary, 1, INFINITY, 0); });
    const double huge[5] = {0, 0, 100000, 100000, 1};
    rejects([&]() { drm_output::box_points(huge, 1, 0, 0); });
    if (argc == 2) {
        std::ifstream input(argv[1]);
        auto supplied = drm_output::read_points(input);
        std::cout << "Supplied file: " << supplied.size() << " points\n";
    }
    std::cout << "DRM geometry and input regressions passed\n";
}
