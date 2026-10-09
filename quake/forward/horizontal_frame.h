/* Horizontal coordinates in metres; azimuth clockwise from the parent +X
 * (UTM grid north for geographic frames). Z and depth are unchanged. */
#ifndef HORIZONTAL_FRAME_H
#define HORIZONTAL_FRAME_H

#include <math.h>

struct horizontal_frame_t {
    double c, s;
};

/* Evaluate trigonometry once at setup. Snap cardinal rotations exactly so
 * boundary coordinates do not drift outside their grid through roundoff. */
static inline horizontal_frame_t horizontal_frame(double degrees) {
    double angle = remainder(degrees, 360.0);
    if (angle == 0.0) return {1.0, 0.0};
    if (angle == 90.0) return {0.0, 1.0};
    if (angle == -90.0) return {0.0, -1.0};
    if (fabs(angle) == 180.0) return {-1.0, 0.0};
    angle *= 3.14159265358979323846 / 180.0;
    return {cos(angle), sin(angle)};
}

static inline void horizontal_to_local(horizontal_frame_t frame,
                                        double north, double east,
                                        double* x, double* y) {
    *x = frame.c * north + frame.s * east;
    *y = -frame.s * north + frame.c * east;
}

static inline void horizontal_to_parent(horizontal_frame_t frame,
                                         double x, double y,
                                         double* north, double* east) {
    *north = frame.c * x - frame.s * y;
    *east = frame.s * x + frame.c * y;
}

#endif
