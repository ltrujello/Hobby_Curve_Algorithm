#ifndef CPP_HOBBYCURVE_H
#define CPP_HOBBYCURVE_H
#include <vector>
#include "HobbyPoint.h"

class HobbyCurve {
public:
    // Constructor
    HobbyCurve(std::vector<std::pair<double,double>> &points, double tension, bool cyclic, double begin_curl, double end_curl);
    // Destructor
    ~HobbyCurve();
    // Algorithm procedures
    std::vector<HobbyPoint> get_points();
    std::vector<std::pair<double, double>> get_ctrl_pts();
    void calculate_d_vals();
    void calculate_psi_vals();
    void calculate_theta_vals();
    void calculate_phi_vals();
    void calculate_ctrl_pts();

private:
    std::vector<HobbyPoint> points;
    std::vector<std::pair<double, double>> ctrl_pts;
    bool cyclic;
    double begin_curl;
    double end_curl;
    int num_points;
};

double velocity(double theta, double phi);


#endif //CPP_HOBBYCURVE_H
