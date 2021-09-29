#include <Eigen>
#include <cmath>
#include <vector>
#include "HobbyCurve.h"
#include "HobbyPoint.h"

// Constructor
HobbyCurve::HobbyCurve(std::vector<std::pair<double, double>> &input_points,
                       double tension,
                       bool cyclic,
                       double begin_curl,
                       double end_curl) :
                       cyclic(cyclic),
                       begin_curl(begin_curl),
                       end_curl(end_curl) {
    if (input_points.size() < 2)
        throw std::runtime_error("Hobby Algorithm needs more thant two points");
    for (auto &point : input_points) {
        HobbyPoint new_point{
                {point.first, point.second},
                1.0 / tension,
                1.0 / tension,
        };
        points.emplace_back(new_point);
    }
    num_points = points.size();
}

// Destructor
HobbyCurve::~HobbyCurve() = default;

// Getter
std::vector<HobbyPoint> HobbyCurve::get_points() {
    return points;
}

// Main method
std::vector<std::pair<double,double>> HobbyCurve::get_ctrl_pts() {
    calculate_d_vals();
    calculate_psi_vals();
    calculate_theta_vals();
    calculate_phi_vals();
    calculate_ctrl_pts();
    return ctrl_pts;
}

void HobbyCurve::calculate_d_vals() {
    /* Calculates the pairwise distances between the points. */
    // Skip last point if path is non-cyclic
    int end = cyclic ? num_points : num_points - 1;
    for (int i = 0; i < end; i++) {
        HobbyPoint &z_i = points[i];
        const HobbyPoint z_j = points[(i + 1) % num_points];
        z_i.d_val = std::abs((z_i.cmplx - z_j.cmplx));
    }
}

void HobbyCurve::calculate_psi_vals() {
    /* Calculates the psi values by subtracting pairwise phases. */
    // Skip first and last point if path is non-cyclic
    int start = cyclic ? 0 : 1;
    int end = cyclic ? num_points : num_points - 1;
    for (int i = start; i < end; i++) {
        const HobbyPoint &z_h = (i == 0 ? points[num_points - 1] : points[i - 1]);
        HobbyPoint &z_i = points[i];
        const HobbyPoint &z_j = points[(i + 1) % num_points];
        const std::complex<double> polygonal_turn = (z_j.cmplx - z_i.cmplx) / (z_i.cmplx - z_h.cmplx);
        z_i.psi = std::atan2(polygonal_turn.imag(), polygonal_turn.real());
    }
}

void HobbyCurve::calculate_theta_vals() {
    /* Calculates the theta values by creating a linear system whose solutions are the values. */
    Eigen::VectorXd A(num_points); // Inappropriate names, but they mirror Knuth's notation.
    Eigen::VectorXd B(num_points);
    Eigen::VectorXd C(num_points);
    Eigen::VectorXd D(num_points);
    Eigen::VectorXd R(num_points);

    // Calculate the entries of the five vectors.
    // Skip first and last point if path is non-cyclic.
    int start = cyclic ? 0 : 1;
    int end = cyclic ? num_points : num_points - 1;
    for (int i = start; i < end; i++) {
        const HobbyPoint &z_h = (i == 0 ? points[num_points - 1] : points[i - 1]);
        const HobbyPoint &z_i = points[i];
        const HobbyPoint &z_j = points[(i + 1) % num_points];

        A(i) = z_h.alpha / (std::pow(z_i.beta, 2) * z_h.d_val);
        B(i) = (3 - z_h.alpha) / (std::pow(z_i.beta, 2) * z_h.d_val);
        C(i) = (3 - z_j.beta) / (std::pow(z_i.alpha, 2) * z_i.d_val);
        D(i) = z_j.beta / (std::pow(z_i.alpha, 2) * z_i.d_val);
        R(i) = -B(i) * z_i.psi - D(i) * z_j.psi;
    }

    // Set up matrix M such that the soln. Mx = R are the theta values.
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(num_points, num_points);
    for (int i = start; i < end; i++) {
        // Fill i-th row of M
        i == 0 ? M(i, num_points - 1) = A(i) : M(i, i - 1) = A(i);
        M(i, i) = B(i) + C(i);
        M(i, (i + 1) % num_points) = D(i);
    }

    // Special formulas for first and last rows of M with non-cyclic paths.
    if (not cyclic) {
        // First row of M
        double alpha_0 = points[0].alpha;
        double beta_1 = points[1].beta;
        double xi_0 = (std::pow(alpha_0, 2) * begin_curl) / std::pow(beta_1, 2);
        M(0, 0) = alpha_0 * xi_0 + 3 - beta_1;
        M(0, 1) = (3 - alpha_0) * xi_0 + beta_1;
        R(0) = -((3 - alpha_0) * xi_0 + beta_1) * points[1].psi;
        // Last row of M
        double alpha_n_1 = points[num_points - 2].alpha;
        double beta_n = points[num_points - 1].beta;
        double xi_n = (std::pow(beta_n, 2) * end_curl) / std::pow(alpha_n_1, 2);
        M(num_points - 1, num_points - 2) = (3 - beta_n) * xi_n + alpha_n_1;
        M(num_points - 1, num_points - 1) = (beta_n * xi_n + 3 - alpha_n_1);
        R(num_points - 1) = 0;
    }
    // Solve for theta values
    Eigen::VectorXd thetas = M.partialPivLu().solve(R);
    for (int i = 0; i < num_points; i++) {
        points[i].theta = thetas(i);
    }
}

void HobbyCurve::calculate_phi_vals() {
    /* Calculates the phi_k values via the relationship theta_k + phi_k + psi_k = 0. */
    for (auto &point : points)
        point.phi = -(point.psi + point.theta);
}

void HobbyCurve::calculate_ctrl_pts() {
    /* Calculates the Bezier control points from z_i to z_{i+1}. */
    // Skip last point if path is non-cyclic
    int end = cyclic ? num_points : num_points - 1;
    for (int i = 0; i < end; i++) {
        const HobbyPoint z_i = points[i];
        const HobbyPoint z_j = points[(i + 1) % num_points];
        double rho_coefficient = z_i.alpha * velocity(z_i.theta, z_j.phi);
        double sigma_coefficient = z_j.beta * velocity(z_j.phi, z_i.theta);

        const std::complex<double> ctrl_point_a =
                z_i.cmplx + (1. / 3.) * rho_coefficient * std::polar(1., z_i.theta) * (z_j.cmplx - z_i.cmplx);
        const std::complex<double> ctrl_point_b =
                z_j.cmplx - (1. / 3.) * sigma_coefficient * std::polar(1., -z_j.phi) * (z_j.cmplx - z_i.cmplx);
        ctrl_pts.emplace_back(std::make_pair(ctrl_point_a.real(), ctrl_point_a.imag()));
        ctrl_pts.emplace_back(std::make_pair(ctrl_point_b.real(), ctrl_point_b.imag()));
    }
}

double velocity(double theta, double phi) {
    /* Metafont's velocity function. */
    double numerator = 2 + std::sqrt(2.) * (std::sin(theta) - (1. / 16.) * std::sin(phi)) *
                           (std::sin(phi) - (1. / 16) * std::sin(theta)) * (
                                   std::cos(theta) - std::cos(phi));
    double denominator = (1 + (1. / 2.) * (std::sqrt(5.) - 1) * std::cos(theta) +
                          (1. / 2.) * (3 - std::sqrt(5.)) * std::cos(phi));
    return numerator / denominator;
}