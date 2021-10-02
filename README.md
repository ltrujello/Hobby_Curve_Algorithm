# Hobby Curve Algorithm
This contains separate Python, Javascript, and C++ implementations of John Hobby's curve drawing algorithm. The algorithm takes in a list of points and calculates the control points of a sequence of Bezier splines which pass through the given points. The resulting shape is very nice and suitable for mathematical drawings.

I only meant to create a Python implementation, but I needed one in Javascript to run in a browser so I wrote one in Javascript. 
Then I thought I could use a C++ wrapped implementation for Javascript. So now all three are here.

## Usage
You probably want to use the Python implementation.

My implementation is `hobby.py`, so use that. `python2_hobby.py` is an older working version I found on 
TeX StackExchange that I used to compare mine against.

The main and only function you need to use is
```python
def hobby_ctrl_points(points: list[tuple], 
                      tension: float=1, 
                      cyclic: bool=True, 
                      begin_curl: float=1,
                      end_curl: float=1, 
                      debug: bool=False) -> list[tuple]:
    """Calculates all cubic Bezier control points, based on John Hobby's algorithm, and pretty prints them."""
```
* **Input:** Supply the points that you want to interpolate in the `points` argument and specify any additional parameters you want. As you can tell, default values are also provided.

* **Output:** This function returns a list of the control points that were calculated and writes to `sys.std.out` (your terminal, if you run it in a terminal) and pretty prints the control points. By "pretty prints", I mean it aligns the calculated tuples up in columns for easy viewing.

* **Debug Mode:** Running this with `debug=True` will print out all of the parameters associated with each point that are used in the algorithm.
This is useful for debugging if you are writing your own implementation. Because this algorithm relies on independent quantities to be 
calculated correctly, if you have an error, this setting makes it easier for you to pinpoint which quantity you are calculating incorrectly and so you 
can catch bugs faster. 

## Example
```python
>>> from hobby import hobby_ctrl_points
>>> points = [(0, 0), (10, 10), (20, 0), (10, -10)]
>>> hobby_ctrl_points(points, tension=3)
(0.0000000000  , 1.8409491661  ) and (8.1590508339  , 10.0000000000 )
(11.8409491661 , 10.0000000000 ) and (20.0000000000 , 1.8409491661  )
(20.0000000000 , -1.8409491661 ) and (11.8409491661 , -10.0000000000)
(8.1590508339  , -10.0000000000) and (0.0000000000  , -1.8409491661 )
```
The Javascript and C++ implementations work similarly.

## Motivation
I originally made this to supersede an existing Python implementation `python2_hobby.py` because it is very bad Python. It's very messy code, overcomplicated, and  it is actually implemented incorrectly in two places that luckily cancel each 
other out (Lines 192, 193, and function `f` on 183). So I felt the need to write a clean and modern version. 

## Testing
Run `pytest test.py` to see if the tests pass. The tests are thorough and I've run them enough times to be confident that
there isn't a bug.
