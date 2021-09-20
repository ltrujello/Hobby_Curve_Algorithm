# Hobby Curve Algorithm
This is a Python implementation of John Hobby's curve drawing algorithm. The algorithm
takes in a sequence of points and calculates the control points which describe Bezier curves 
that interpolate the sequence very nicely. 

## Usage
My implementation is `hobby.py`, so use that. `python2_hobby.py` is an older working version I used to compare against.

The main and only function you need to use is
```python
def hobby_ctrl_points(points: list[tuple], tension: float = 1, cyclic: bool = True, begin_curl: float = 1,
                      end_curl: float = 1) -> None:
    """Calculates all cubic Bezier control points, based on John Hobby's algorithm, and pretty prints them."""
```
Supply your points that you want to interpolate in the `points` argument and specify any additional parameters you want, 
although sensible default values are used.

This function writes to `sys.std.out` (your terminal, if you run it in a terminal) and pretty prints the calculated 
control points. By "pretty prints", I mean it aligns the calculated tuples up in columns for easy viewing.

## Example
```python
>>> from hobby import hobby_ctrl_points
>>> points = [(0, 0), (10, 10), (20, 0), (10, -10)]
>>> hobby_ctrl_points(points, cyclic=True, tension=3)
(0.0000000000  , 1.8409491661  ) and (8.1590508339  , 10.0000000000 )
(11.8409491661 , 10.0000000000 ) and (20.0000000000 , 1.8409491661  )
(20.0000000000 , -1.8409491661 ) and (11.8409491661 , -10.0000000000)
(8.1590508339  , -10.0000000000) and (0.0000000000  , -1.8409491661 )
```


## Motivation
The existing Python code `python2_hobby.py` is very bad Python. It's very messy code, overcomplicated, and 
it is actually implemented incorrectly in two places that luckily cancel each 
other out (Lines 192, 193, and function `f` on 183). It makes a good object oriented abstraction, but instead of writing class methods, 
it creates independent functions that accept class instances which change the attributes, which is just 
an error-prone class method with extra steps. So I felt the need to write a clean and modern version. 

## Testing
Run `pytest test.py` to see if the tests pass. The tests are thorough and I've run them enough times to be confident that
there isn't a bug.
