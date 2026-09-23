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
    // Unequal lengths catch swapped face axes; each face is a complete grid.
    const drm_output::BoxSpec spec{{10, 20, 4, 6, 3, 1, 0}};
    const auto box = drm_output::box_faces(spec);
    const int shapes[10][2] = {{7,4},{7,4},{5,4},{5,4},{5,7},
                               {5,3},{5,3},{3,3},{3,3},{3,5}};
    // Exterior: 28+28+20+20+35=131; interior:15+15+9+9+15=63.
    assert(box.points.size() == 194);
    size_t offset = 0;
    for (int face = 0; face < 10; ++face) {
        const auto& grid = box.faces[face];
        assert(grid.rows == shapes[face][0] && grid.columns == shapes[face][1]);
        assert(grid.offset == offset);
        const int side = face / 5, kind = face % 5;
        const double xmin = 8 + side, xmax = 12 - side;
        const double ymin = 17 + side, ymax = 23 - side;
        const double depth = 3 - side;
        for (int row = 0; row < grid.rows; ++row)
            for (int col = 0; col < grid.columns; ++col) {
                const auto& point = box.points[offset + row*grid.columns + col];
                if (kind < 2) {
                    assert(point[0] == (kind == 0 ? xmin : xmax));
                    assert(point[1] == ymin + row && point[2] == col);
                } else if (kind < 4) {
                    assert(point[0] == xmin + row);
                    assert(point[1] == (kind == 2 ? ymin : ymax) && point[2] == col);
                } else {
                    assert(point[0] == xmin + row && point[1] == ymin + col && point[2] == depth);
                }
            }
        offset += grid.rows * grid.columns;
    }
    assert(offset == box.points.size());
    const std::set<drm_output::Point> unique(box.points.begin(), box.points.end());
    assert(unique.size() < box.points.size()); // Face edges/corners intentionally repeat.
    assert(unique.count(drm_output::Point{{10,20,0}}) == 0); // No top lid.
    assert(unique.count(drm_output::Point{{10,20,3}}) == 1);
    assert(unique.count(drm_output::Point{{10,20,2}}) == 1);
    for (double strike : {90.0, -90.0, 37.0, 397.0}) {
        auto rotated_spec = spec;
        rotated_spec[6] = strike;
        const auto rotated = drm_output::box_faces(rotated_spec);
        const double a = std::remainder(strike, 360.0) * std::acos(-1.0) / 180;
        for (size_t i = 0; i < box.points.size(); ++i) {
            const double x = box.points[i][0] - 10, y = box.points[i][1] - 20;
            const auto& p = rotated.points[i];
            assert(std::abs((p[0]-10)*std::cos(a) + (p[1]-20)*std::sin(a) - x) < 1e-12);
            assert(std::abs(-(p[0]-10)*std::sin(a) + (p[1]-20)*std::cos(a) - y) < 1e-12);
            assert(p[2] == box.points[i][2]);
        }
    }
    const auto example = drm_output::box_faces(drm_output::BoxSpec{{0,0,100,100,50,1,0}});
    assert(example.points.size() == 60406);
    assert(example.faces[0].rows == 101 && example.faces[0].columns == 51);
    assert(example.faces[5].rows == 99 && example.faces[5].columns == 50);
    assert((example.points[0] == drm_output::Point{{-50,-50,0}}));
    assert((example.points[example.faces[5].offset] == drm_output::Point{{-49,-49,0}}));
    for (double angle : {0.0, 90.0, 180.0, 270.0}) {
        const auto cardinal = drm_output::box_faces(drm_output::BoxSpec{{50,50,100,100,50,1,angle}});
        for (const auto& p : cardinal.points) assert(p[0] >= 0 && p[1] >= 0);
    }
    const auto fractional = drm_output::box_faces(drm_output::BoxSpec{{0,0,0.3,0.5,0.4,0.1,0}});
    for (int side = 0; side < 2; ++side) {
        const auto& front = fractional.faces[side*5];
        const auto& bottom = fractional.faces[side*5+4];
        // The front-bottom edge has identical coordinates in both face arrays.
        for (int row = 0; row < front.rows; ++row)
            assert(fractional.points[front.offset + row*front.columns + front.columns-1] ==
                   fractional.points[bottom.offset + row]);
    }
    for (auto bad : {drm_output::BoxSpec{{0,0,4,6,3,0,0}},
                     drm_output::BoxSpec{{0,0,2,6,3,1,0}},
                     drm_output::BoxSpec{{0,0,4,2,3,1,0}},
                     drm_output::BoxSpec{{0,0,4,6,1,1,0}},
                     drm_output::BoxSpec{{0,0,4,6,3,0.7,0}},
                     drm_output::BoxSpec{{0,0,100000,100000,3,1,0}},
                     drm_output::BoxSpec{{0,0,4,6,3,1,INFINITY}}})
        rejects([&]() { drm_output::box_faces(bad); });
    std::istringstream table("number_output_drm_boxes = 2\noutput_drm_boxes =\n"
        "# center_x center_y lx ly depth h strike\n10 20 4 6 3 1 0\n\n30 40 6 4 4 1 37 # rotated\nnext_key = 1\n");
    const auto specs = drm_output::read_boxes(table, "output_drm_boxes", 2);
    assert(specs.size() == 2 && specs[0] == spec && specs[1][6] == 37);
    for (const char* bad : {"", "output_drm_boxes =\n", "output_drm_boxes = 0 0 4 6 3 1 0",
        "output_drm_boxes =\n0 0 4 6 3 1", "output_drm_boxes =\n0 0 4 6 3 1 0 8",
        "output_drm_boxes =\n0 0 4 6 3 1 nan", "output_drm_boxes =\n0 0 4 6 3 1 0\n0 0 4 6 3 1 0"})
        rejects([bad]() { std::istringstream in(bad); drm_output::read_boxes(in, "output_drm_boxes", 1); });
    if (argc == 2) {
        std::ifstream input(argv[1]);
        auto supplied = drm_output::read_points(input);
        std::cout << "Supplied file: " << supplied.size() << " points\n";
    }
    std::cout << "DRM geometry and input regressions passed\n";
}
