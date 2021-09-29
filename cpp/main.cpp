#include <iostream>
#include "HobbyCurve.h"

void hobby_ctrl_points(std::vector<std::pair<double, double>> points,
                       double tension,
                       bool cyclic,
                       double begin_curl,
                       double end_curl);

int main() {
    std::vector<std::pair<double, double>> points{
            {0.,  0.,},
            {10., 15.},
            {20., 0.},
            {10., -10.}
    };
    hobby_ctrl_points(points, 3, false, 1, 1);

    return 0;
}


void hobby_ctrl_points(std::vector<std::pair<double, double>> points,
                       double tension,
                       bool cyclic,
                       double begin_curl,
                       double end_curl) {
    HobbyCurve curve(points, tension, cyclic, begin_curl, end_curl);
    std::vector<std::pair<double, double>> ctrl_pts = curve.get_ctrl_pts();
    int i = 0;
    while (i < ctrl_pts.size()) {
        double x_1 = ctrl_pts[i].first;
        double y_1 = ctrl_pts[i].second;
        double x_2 = ctrl_pts[i + 1].first;
        double y_2 = ctrl_pts[i + 1].second;
        std::cout << x_1 << " " << y_1 << " and " << x_2 << " " << y_2 << std::endl;
        i += 2;
    }
}
