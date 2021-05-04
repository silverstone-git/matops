# matops
A Simple Matrix Manipulation Library in Python which uses the better column-number,row-number convention 


<b>USAGE</b>

If Python interpreter is installed (on any os), open it up and type:

> import matops<br>
> matops.run()


<b>USE AS A LIBRARY:</b>

> from matops import *<br>
> printmat(inverse([[1, 2, 3], [4, 5, 6], [7, 8, 4]]))<br>

To find the inverse of :

----------
> 1  4  7 <br>
> 2  5  8 <br>
> 3  6  4 
----------

The program has a roundmat function and other parameters for representability purposes
> printmat(roundmat(inverse([[1, 2, 3], [0+1e-18, 9, 10], [0-1e-17, 0-1e-12, 90]]), roundto = 7)), neatness = 10)<br>

>1.0          -0.0             0.0<br>
>-0.2222222   0.1111111        0.0<br>
>-0.008642    -0.0123457       0.0111111<br>


> printmat(roundmat(coremult(gauss_inv([[0, 43, 0, 11], [12, 0, 48, 0], [99, 49, 23, 0], [45, 0, 10, 50]]), [[0, 43, 0, 11], [12, 0, 48, 0], [99, 49, 23, 0], [45, 0, 10, 50]])))<br>

>1.0    0.0   0.0   -0.0<br>
>0.0    1.0   0.0   0.0<br>
>0.0    0.0   1.0   -0.0<br>
>-0.0   0.0   0.0   1.0


