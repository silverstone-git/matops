# matops
A Simple Matrix Manipulation Library in Python which uses the better column-number,row-number convention 


<b>USAGE</b>

*Windows users can just run the executable in the folder to use it as a program instead

If Python interpreter is installed (on any os), open it up and type:

>>> import matops
>>> matops.run()


<b>USE AS A LIBRARY:</b>

>>> from matops import *
>>> printmat(inverse([[1, 2, 3], [4, 5, 6], [7, 8, 4]]))
to find the inverse of -->
___________
| 1  4  7 |
| 2  5  8 |
| 3  6  4 |
-----------

The program has a roundmat function and other parameters for representability purposes
>>> printmat(roundmat(inverse([[1, 2, 3], [0+1e-18, 9, 10], [0-1e-17, 0-1e-12, 90]]), roundto = 7)), neatness = 10)
--------------------------------------------
1.0          -0.0         0.0
-0.2222222   0.1111111    0.0
-0.008642    -0.0123457   0.0111111
--------------------------------------------

>>> printmat(roundmat(coremult(gauss_inv([[0, 43, 0, 11], [12, 0, 48, 0], [99, 49, 23, 0], [45, 0, 10, 50]]), [[0, 43, 0, 11], [12, 0, 48, 0], [99, 49, 23, 0], [45, 0, 10, 50]])))
-----------------------------------
1.0    0.0   0.0   -0.0
0.0    1.0   0.0   0.0
0.0    0.0   1.0   -0.0
-0.0   0.0   0.0   1.0
-----------------------------------

