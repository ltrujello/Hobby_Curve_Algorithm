import numpy as np
from hobby import HobbyCurve
from python2_hobby import mp_to_tikz

"""Tests to compare the control points generated by the old python code and the new code I wrote. 

Points are selected uniformly from [10, -10]. For each such random selection, we  
independently and dependently specify and vary the available options and numerical parameters such as the 
path being cyclic, the tension, and the curl.
"""

# Run `pytest test.py`.

def test_cyclic():
    """Test algorithms against each other for different points with cyclic specified."""
    run_and_compare_algorithms(cyclic=True)


def test_cyclic_and_tension():
    """Test algorithms against each other with cyclic, tension specified."""
    for tension in np.arange(0.5, 5.5, 0.5):
        run_and_compare_algorithms(tension=tension, cyclic=True)


def test_noncyclic():
    """Test algorithms against each other with noncyclic specified."""
    run_and_compare_algorithms(cyclic=False)


def test_noncyclic_and_tension():
    """Test algorithms against each other with noncyclic, tension specified."""
    for tension in np.arange(0.5, 5.5, 0.5):
        run_and_compare_algorithms(tension=tension, cyclic=False)


def test_noncyclic_and_tension_and_curls():
    """Test algorithms against each other with noncyclic, tension, curls specified."""
    for curl in np.arange(0.5, 5.5, 0.5):
        for tension in np.arange(0.5, 5.5, 0.5):
            run_and_compare_algorithms(tension=tension, cyclic=False, curl=curl)


def generate_points(num_points):
    """Generate a random number of points from -10 to 10."""
    points = []
    for _ in range(num_points):
        x = np.random.uniform(-10, 10)
        y = np.random.uniform(-10, 10)
        points.append((x, y))
    return points


def run_and_compare_algorithms(tension=1, cyclic=True, curl=1):
    """Wrapper for testing running the algorithms against each other."""
    for num_points in range(10, 21):
        input_points = generate_points(num_points)
        old_ctrl_points = mp_to_tikz(
            " .. ".join([str(point) for point in input_points]) + f" {'.. cycle' if cyclic else ''} ",
            options=f"tension ={tension}, curl ={curl}")
        new_ctrl_points = HobbyCurve(input_points, cyclic=cyclic, tension=tension, begin_curl=curl,
                                     end_curl=curl).get_crl_points()
        compare_ctrl_points(old_ctrl_points, new_ctrl_points)


def compare_ctrl_points(old_ctrl_pts, new_ctrl_pts):
    """Observe if there is significant distance between the control points calculated via
    the old code and the new code."""
    max_error = 1e-8  # If there's an implementation error, it will be very off. Minor differences are due to libraries.
    assert len(old_ctrl_pts) == len(new_ctrl_pts), f"{len(old_ctrl_pts)}, {len(new_ctrl_pts)}"

    for i in range(len(old_ctrl_pts)):
        x_1, y_1 = old_ctrl_pts[i]
        x_2, y_2 = new_ctrl_pts[i]
        x_error, y_error = abs(x_1 - x_2), abs(y_1 - y_2)
        # Main test
        assert x_error < max_error, x_error
        assert y_error < max_error, y_error
