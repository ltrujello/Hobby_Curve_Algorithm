import numpy as np
import cmath

"""
Implementation of John Hobby's Bezier curve algorithm in Python.
The algorithm is very simple and efficient. Details are on page 112, 113 in Knuth's METAFONT: The Program.
"""


class HobbyPoint(complex):
    """A class for associating numerical quantities from Hobby's algorithm with points that appear on a Hobby curve.
    We subclass `complex` to perform complex arithmetic with points on the Hobby curve, as required in the algorithm."""

    def __new__(cls, x: float, y: float, tension: float) -> 'HobbyPoint':
        return super().__new__(cls, x, y)

    def __init__(self, x: float, y: float, tension: float) -> None:
        self.x = x
        self.y = y
        # In what follows, we use Knuth's notation in our variable names.
        self.alpha = 1 / tension
        self.beta = 1 / tension
        self.d_val = 0  # Distance between this point and next.
        self.theta = 0  # Angle of polygonal line from this point to next.
        self.phi = 0  # Offset angle.
        self.psi = 0  # Another offset angle.

    def __repr__(self) -> str:
        return f"{(self.x, self.y)}"


class HobbyCurve:
    """A class for calculating the control points required to draw a Hobby curve."""

    def __init__(self, points: list[tuple], tension: float = 1, cyclic: bool = False, begin_curl: float = 1,
                 end_curl: float = 1) -> None:
        self.points = [HobbyPoint(*point, tension) for point in points]
        self.ctrl_pts = []
        self.is_cyclic = cyclic
        self.begin_curl = begin_curl
        self.end_curl = end_curl
        self.num_points = len(points)

    def get_crl_points(self) -> list[tuple]:
        """Calculates and returns all of the control points of the Hobby curve."""
        self.calculate_d_vals()
        self.calculate_psi_vals()
        self.calculate_theta_vals()
        self.calculate_phi_vals()
        self.ctrl_pts = self.calculate_ctrl_pts()
        return self.ctrl_pts

    def calculate_d_vals(self) -> None:
        """Calculates the pairwise distances between the points."""
        # Skip last point if path is non-cyclic
        point_inds = range(self.num_points) if self.is_cyclic else range(self.num_points - 1)
        for i in point_inds:
            z_i = self.points[i % self.num_points]
            z_j = self.points[(i + 1) % self.num_points]
            z_i.d_val = abs(z_i - z_j)

    def calculate_psi_vals(self) -> None:
        """Calculates the psi values by subtracting pairwise phases."""
        # Skip first and last point if path is non-cyclic
        point_inds = range(self.num_points) if self.is_cyclic else range(1, self.num_points - 1)
        for i in point_inds:
            z_h = self.points[i - 1]
            z_i = self.points[i]
            z_j = self.points[(i + 1) % len(self.points)]
            polygonal_turn = (z_j - z_i) / (z_i - z_h)
            z_i.psi = np.arctan2(polygonal_turn.imag, polygonal_turn.real)

    def calculate_theta_vals(self) -> None:
        """Calculates the theta values by creating a linear system whose solutions are the values."""
        A = np.zeros(self.num_points)  # Inappropriate names, but they mirror Knuth's notation.
        B = np.zeros(self.num_points)
        C = np.zeros(self.num_points)
        D = np.zeros(self.num_points)
        R = np.zeros(self.num_points)

        # Calculate the entries of the five vectors.
        # Skip first and last point if path is non-cyclic.
        point_ind = range(self.num_points) if self.is_cyclic else range(1, self.num_points - 1)
        for i in point_ind:
            z_h = self.points[i - 1]
            z_i = self.points[i]
            z_j = self.points[(i + 1) % len(self.points)]

            A[i] = z_h.alpha / (z_i.beta ** 2 * z_h.d_val)
            B[i] = (3 - z_h.alpha) / (z_i.beta ** 2 * z_h.d_val)
            C[i] = (3 - z_j.beta) / (z_i.alpha ** 2 * z_i.d_val)
            D[i] = z_j.beta / (z_i.alpha ** 2 * z_i.d_val)
            R[i] = -B[i] * z_i.psi - D[i] * z_j.psi

        # Set up matrix M such that the soln. Mx = R are the theta values.
        M = np.zeros((self.num_points, self.num_points))
        for i in range(self.num_points):
            # Fill i-th row of M
            M[i][i - 1] = A[i]
            M[i][i] = B[i] + C[i]
            M[i][(i + 1) % self.num_points] = D[i]

        # Special formulas for first and last rows of M with non-cyclic paths.
        if not self.is_cyclic:
            # First row of M
            alpha_0 = self.points[0].alpha
            beta_1 = self.points[1].beta
            xi_0 = (alpha_0 ** 2 * self.begin_curl) / beta_1 ** 2
            M[0][0] = alpha_0 * xi_0 + 3 - beta_1
            M[0][1] = (3 - alpha_0) * xi_0 + beta_1
            R[0] = -((3 - alpha_0) * xi_0 + beta_1) * self.points[1].psi
            # Last row of M
            alpha_n_1 = self.points[-2].alpha
            beta_n = self.points[-1].beta
            xi_n = (beta_n ** 2 * self.end_curl) / alpha_n_1 ** 2
            M[-1][-2] = (3 - beta_n) * xi_n + alpha_n_1
            M[-1][-1] = (beta_n * xi_n + 3 - alpha_n_1)
            R[-1] = 0

        # Solve for theta values.
        thetas = np.linalg.solve(M, R)
        for i, point in enumerate(self.points):
            point.theta = thetas[i]

    def calculate_phi_vals(self) -> None:
        """Calculates the phi_k values via the relationship theta_k + phi_k + psi_k = 0."""
        for point in self.points:
            point.phi = - (point.psi + point.theta)

    def calculate_ctrl_pts(self) -> list[tuple]:
        """Calculates the Bezier control points from z_i to z_{i+1}."""
        ctrl_pts = []
        # Skip last point if path is non-cyclic
        point_inds = range(self.num_points) if self.is_cyclic else range(self.num_points - 1)
        for i in point_inds:
            z_i = self.points[i]
            z_j = self.points[(i + 1) % self.num_points]
            rho_coefficient = z_i.alpha * velocity(z_i.theta, z_j.phi)
            sigma_coefficient = z_j.beta * velocity(z_j.phi, z_i.theta)

            ctrl_pt_a = z_i + (1 / 3) * rho_coefficient * cmath.exp(complex(0, z_i.theta)) * (z_j - z_i)
            ctrl_pt_b = z_j - (1 / 3) * sigma_coefficient * cmath.exp(complex(0, -z_j.phi)) * (z_j - z_i)
            ctrl_pts.append((ctrl_pt_a.real, ctrl_pt_a.imag))
            ctrl_pts.append((ctrl_pt_b.real, ctrl_pt_b.imag))
        return ctrl_pts

    def __repr__(self) -> str:
        cartesian_points = [(point.real, point.imag) for point in self.points]
        return repr(cartesian_points)


def hobby_ctrl_points(points: list[tuple], tension: float = 1, cyclic: bool = False, begin_curl: float = 1,
                      end_curl: float = 1) -> list[tuple]:
    """Calculates all cubic Bezier control points, based on John Hobby's algorithm, and pretty prints them."""
    curve = HobbyCurve(points, tension=tension, cyclic=cyclic, begin_curl=begin_curl, end_curl=end_curl)
    ctrl_points = curve.get_crl_points()

    # Calculate whitespace padding for pretty print.
    max_pad = 0
    for ctrl_point in ctrl_points:
        x, y = ctrl_point
        # Calculate number of digits in x, y before decimal, and take the max for nice padding.
        padding = max(1 if abs(x) <= 0.1 else int(np.ceil(np.log10(abs(x)))) + 1,
                      1 if abs(y) <= 0.1 else int(np.ceil(np.log10(abs(y)))) + 1)
        if max_pad < padding:
            max_pad = padding

    # Pretty print control points.
    precision = 10
    space = precision + max_pad + 1  # +1 for negative sign
    i = 0
    while i < len(ctrl_points) - 1:
        x_1, y_1 = ctrl_points[i]
        x_2, y_2 = ctrl_points[i + 1]
        print(f"({x_1:<{space}.{precision}f}, {y_1:<{space}.{precision}f}) "
              f"and "
              f"({x_2:<{space}.{precision}f}, {y_2:<{space}.{precision}f})")
        i += 2
    return ctrl_points


def velocity(theta: float, phi: float) -> float:
    """Metafont's velocity function."""
    numerator = 2 + np.sqrt(2) * (np.sin(theta) - (1 / 16) * np.sin(phi)) * (np.sin(phi) - (1 / 16) * np.sin(theta)) * (
            np.cos(theta) - np.cos(phi))
    denominator = (1 + (1 / 2) * (np.sqrt(5) - 1) * np.cos(theta) + (1 / 2) * (3 - np.sqrt(5)) * np.cos(phi))
    return numerator / denominator
