#!/usr/bin/python3
'''
___USAGE___
>>> import matops
>>> help(matops)

example:
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
'''

"""
___DESCRIPTION___
Author => Aryan Sidhwani, github.com/silverstone-git
Matrix Functions
Functions to solve various problems in Linear Algebra, or any problem converted to that form
The whole module uses [column number], [row number] convention in the same order  
If you are seeing details of an internal process in cmd, the function being used is still in development

Project Started: 9th April 2020
First Stable: 2nd May 2020
Last Updated: 2nd July 2021

#NoImportsGangGang
"""

pi = 3.141592653589793238462643383279502884197169399375105820974944


#FUNCTIONS

def mod(x):
    """
    Returns the absolute value of input number x
    """
    if x < 0:
        return -x
    else:
        return x


def atan(z, n = 1000):
    """
    Arc Tangent or Inverse Tangent Function, z is the usual value input, n is number of iterations, or accuracy
    """
    s = 0
    for ct in range(n+1):
        #print(ct)
        try:
            s += ((-1) ** ct) * ((z ** (2*ct + 1)) / (2*ct + 1))   #formula derived from 1/1-x G.P. expansion transformed to 1/1+x2 and integrate both sides
        except:
            print("Overflow error was encountered")
            break
        
    return s


def factorial(n):
    """
    Computes and Returns the factorial of an integer
    """
    if n == 0:
        return 1
    try:
        n = int(n)
    except:
        print("Factorial of floats, etc are not yet developed, planning on using gamma function..")
        return None
    pr = 1
    for ct in range(1, n + 1):
        pr *= ct
    return pr


def sin(x, n=85):
    """
    Computes and returns the sine of an angle given in radians using taylor approximation
    """
    x = x % (2 * pi) 
    taylor = 0
    for ct in range(n):
        taylor += (( (-1) ** ct ) * (x ** ((2 * ct) + 1))) / (factorial((2 * ct) + 1))
    return taylor


def cos(x, tol = 1e-12):
    """
    Returns the cosine of an angle given in radians using sine function
    """
    sinsquare = sin(x) ** 2
    if 1 - sinsquare >= 0:
        return (1 - (sin(x) ** 2)) ** 0.5
    elif 1 - sinsquare > -tol:
        return 0
    else:
        print("sin function diverged, hqhq cos x is 1 anyway")
        return 1


def tan(x):
    """
    Returns the tangent of an angle given in radians using the sine and cosine
    """
    return sin(x)/cos(x)


def read():
    """
    No arguments
    Returns a List of lists of Column Vectors, used as a Matrix
    Call it a Second Order Composite List
    """
    print("\n")
    print("While entering a Matrix, if a value is mistyped, type 404e2 to repeat that particular column vector")
    print("If 404e2 is typed at the 1st value, last column vector will be ignored and shall be entered again")
    print("\n")
    n = 404e2
    rep = True
    rep1 = False
    while n == 404e2 or rep:
        try:
            m = int(input("Enter the number of rows of matrix: "))
            n = int(input("Enter the number of columns of matrix: "))
        except:
            print("Please Enter integers")
            rep1 = True
        really = input("Proceed? ")
        if really == "n" or really == "no" or really == "No" or really == "NO" or really == "0":
            rep = True
        elif not rep1:
            rep = False
    mat = []
    ccol = 1
    while ccol <= n:
        cols = []
        for crow in range(1,m+1):
            try:
                x = float(eval(input("Enter the {0} th entry of the vertical vector: ".format(crow))))
            except:
                print("Please enter a decimal or an integer")
                x = 404e2
            if x == 404e2 and crow != 1:
                ccol -= 1
                print("Repeating last column vector")
                break
            elif x == 404e2 and crow == 1:
                mat.pop()
                ccol -=2
                print("Repeating last column vector")
                break
            cols.append(x)
        if x != 404e2:
            mat.append(cols)
        if ccol == n and crow == m:
            really = input("Proceed? ")
            if really == "n" or really == "N" or really == "No" or really == "no" or really == "0":
                mat.pop()
                ccol -= 1
                print("Repeating last column vector")
        ccol += 1
    print("Entered Matrix:")
    printmat(mat)
    return mat


def det(mat):
    """
    A Composite List of Column Vectors as Argument
    Returns the determinant of a Transformation
    """
    if len(mat) == 2 and len(mat[0]) == 2:
        return ((mat[0][0]*mat[1][1])-(mat[1][0]*mat[0][1]))
    else:
        s = 0
        for a in range(0,len(mat)):
            matmod = []
            for b1 in range(0,a):
                matmod.append(mat[b1][1:])
            for b2 in range(a+1,len(mat)):
                matmod.append(mat[b2][1:])

            # recursion used
            s += ((-1)**a)*(mat[a][0])*(det(matmod))
        return s


def refine(mat):
    '''
    A Composite List of Column Vectors as Argument
    Returns the Matrix with Zero Vectors removed
    '''
    for a in range(len(mat)):
        if len(a) == 0:
            mat.pop(a)
    return mat


def trace(m1):
    '''
    A Composite List of Column Vectors as Argument
    Returns the Trace of a Matrix, i.e. sum of L2R Diagonal Terms
    '''
    if len(m1) != len(m1[0]):
        print("Not a square matrix!")
        return None
    else:
        traceres = 0
        for diagelementindex in range(len(m1)):
            traceres += m1[diagelementindex][diagelementindex]
        return traceres


def negatefn(listofcfs):
    '''
    List of numbers defining the coefficients of a polynomial as argument
    Equivalent to the transformation x --> -x, to check for odd or even function
    '''
    if (len(listofcfs) % 2) == 0:
        for index in range(len(listofcfs)):
            listofcfs[index] *= (-1)**(index+1)
    else:
        for index in range(len(listofcfs)):
            listofcfs[index] *= (-1)**(index)
    return listofcfs


def v_norm(vector, p = 2):
    '''
    One Vector Coefficient List argument and an optional p argument (default value 2)
    Known as "|v|p" (p in subscript) or "the norm" in maths
    Returns the calculated norm of the vector
    '''
    try:
        _ = vector[1]
        vector = [vector]
    except:
        pass
    s = 0
    for index in range(len(vector[0])):
        s += (vector[0][index])**p
    return s**(1/p)


def sgnchangefn(listofcfs):
    '''
    Used by descartes_ros() to calculate number of roots
    A list of coefficients of a polynomial as argument
    returns the number of sign changes 
    '''
    def poscheck(coeff):
        #checks if the number if positive
        if coeff > 0:
            pos = True
        else:
            pos = False
        return pos
    sgnchange = 0
    for index in range(1,len(listofcfs)):
        if poscheck(listofcfs[index]) != poscheck(listofcfs[index-1]):
            sgnchange += 1
    return sgnchange


def descartes_ros(listofcfs):
    '''
    Descartes Rule of Signs
    Returns the number of roots of a polynomial, given;
    The list of coefficients of a polynomial
    '''
    return (sgnchangefn(listofcfs) + sgnchangefn(negatefn(listofcfs)))


def solveqn(listofcfs, acc = 3000):
    '''
    One List Argument
    Tries to return the roots of the polynomial equation defined by list of coefficients*
    *Highest Order Coefficient at the 0th Index and the Constant at Last one
    '''
    if listofcfs[0] == 0:
        while True:
            if listofcfs[0] == 0:
                listofcfs.pop(0)
            else:
                break
        print("reduced list: ", listofcfs)
    if len(listofcfs) == 2:
        return -listofcfs[1]/listofcfs[0]
    elif len(listofcfs) == 3:
        return [(((((listofcfs[1]**2)-(4*listofcfs[0]*listofcfs[2]))**(1/2))-listofcfs[1])/(2*listofcfs[0])), ((-(((listofcfs[1]**2)-(4*listofcfs[0]*listofcfs[2]))**(1/2))-listofcfs[1])/(2*listofcfs[0]))]
    if listofcfs[0] != 1:
        a = listofcfs[0]
        for icoeff in range(len(listofcfs)):
            listofcfs[icoeff] /= a
    mat = []
    n = len(listofcfs) - 1
    for icf in range(n):
        if icf == 0 :
            vector = ([0] * (n-1) ) + [-listofcfs[-1] ]
        else:
            vector = ([0] * (icf-1) ) + [1] + ([0] * (n - icf -1)) + [-listofcfs[-icf-1] ]
        #print("vector is: ", vector)
        mat.append(vector)
    return eigenvalues(mat, acc = acc)


def coremult(mat1,mat2):
    '''
    Two Composite Lists of Column Vectors as Arguments
    PreMultiplies the First Matrix to the Second and Returns the Result
    '''
    if len(mat1) != len(mat2[0]):
        print("\n\n")
        print("ERROR: Number of columns of first matrix should be equal to Number of rows of second to Multiply")
        return None
    matmultres = []
    for a in mat2:
        matmultcol = []
        for b in range(len(mat1[0])):
            s = 0
            for c in range(len(mat2[0])):
                s += a[c] * mat1[c][b]
            matmultcol.append(s)
        matmultres.append(matmultcol)
    return(matmultres)


def multiply(mat1a, mat2a):
    '''
    Two Second Order Composite List Arguments
    Returns Multiplied Matrix and asks to multiply further
    '''
    multrep = None
    while multrep != "N" or "n" or "no" or "No" or "NO":
        mat1a = coremult(mat1a,mat2a)
        if mat1a == None:
            return None
        print("Resultant Matrix:")
        printmat(mat1a)
        print("\n\n")
        multrep = input("Do you want to Multiply the result with another matrix??\nY or N >> ")
        if multrep == "n" or multrep == "No" or multrep == "NO" or multrep == "no" or multrep == "N":
            break
        print("\n\n")
        mat2a = read()
    return mat1a


def powermat(m1, exp):
    '''
    A Second Order Composite List of Column Vectors Argument 
    and an integer Argument
    Returns the Matrix raised to the exp-th power
    '''
    multres = m1
    if exp < -1:
        return powermat(inverse(m1), -exp)
    elif exp == 0:
        return idmat(1)
    else:
        for times in range(exp-1):
            multres = coremult(multres, m1)
        return multres


def gauss_inv(thematrix):
    '''
    One Second Order Composite List of Column Vectors Argument 
    step by step converts a Matrix into Reduced Row Echelon Form
    '''

    n = len(thematrix)

    gmat = zeromat(n)
    #print("zero matrix: ", gmat)
    for icol in range(n):
        for irow in range(n):
            gmat[icol][irow] = thematrix[icol][irow]

    #print("edited zero matrix: ", gmat)
    #print("original matrix: ", thematrix)
    
    # inverse calculation starts with the identity matrix
    # and then same transformations as main matrix are done to it
    inv = idmat(len(gmat))
    pivot_found = False
    irow = 0
    icol = 0
    while irow < n and icol < n:
        skipstep = False
        if 'bool' in str(type(pivot_found)):
            if gmat[icol][irow] != 0:

                # finding the first nonzero entry in column
                inv = rowscalmult(inv, irow, 1/gmat[icol][irow])
                gmat = rowscalmult(gmat, irow, 1/gmat[icol][irow])

                pivot_found = irow

                #print("gmat after finding pivot ==> ")
                #printmat(gmat)
                irow += 1
                if irow  == n:
                    break

            else:
                inv = rowswap(inv, irow, irow+1)
                gmat = rowswap(gmat, irow, irow+1)
                skipstep = True
        
        if not skipstep:
            # r ---> r + scalar * pivotrow
            inv = rowmultadd(inv, irow, pivot_found, -gmat[icol][irow])
            gmat = rowmultadd(gmat, irow, pivot_found, -gmat[icol][irow])
            #print("gmat after row scaling and adding ==> ")
            #printmat(gmat)
            irow += 1

        if irow == n:
            icol += 1
            irow = pivot_found + 1
            pivot_found = False


    #print("Front Gaussian Elimination done")
    #print("Starting Gauss-Jordan")
    icol = n-1
    irow = n-2
    upno = 1
    while True:
        inv = rowmultadd(inv, irow, irow + upno, -gmat[icol][irow])
        gmat = rowmultadd(gmat, irow, irow + upno, -gmat[icol][irow])

        #print("gmat after row scaling and adding ==> ")
        #printmat(gmat)

        upno += 1
        irow -= 1

        if irow == -1:
            icol -= 1
            irow = upno - 3
            upno = 1

        if icol == 0:
            break
    #print("gmat value is: ")
    #printmat(gmat)
    #print("\ninv value is: ")
    #printmat(inv)
    #print("\nthematrix value is: ")
    #printmat(thematrix)
    return inv
    

def coreminor(m1, ijcoordinates):
    '''
    A Composite List of Column Vectors and i,j coordinate list as Arguments
    where 
    i = column number and,
    j = row number
    Returns the Minor Matrix (Matrix formed by hiding i-th column and j-th row)
    '''
    minoris = []
    for cols in range(len(m1)):
        if cols != ijcoordinates[0]:
            newcol = m1[cols]
            newcol = newcol[:ijcoordinates[1]] + newcol[ijcoordinates[1]+1:]
            minoris.append(newcol)
    return minoris


def minor(m1):
    '''
    One Second Order Composite List of Column Vectors as Argument
    Returns the Matrix of Determinants of Minors by taking Determinants of Matrices formed by hiding the
    column and row of the current element for every element of the Matrix
    '''
    try:
        if qivar == True:
            pass
    except:
        print("Use det(m1) as determinant for calculating inverse later..")
    if len(m1) == 2 and len(m1[0]) == 2:
        return [det(m1)]
    resmat = []
    for col in range(len(m1)):
        resmatcol = []
        for row in range(len(m1[0])):
            minormat = coreminor(m1, [col, row])
            resmatcol.append(det(minormat))
        resmat.append(resmatcol)
    return resmat


def rowscalmult(m1,r1,sc):
    '''
    A Second Order Composite List of Column Vectors as Argument
    A row number argument (count from zero)
    and a Numerical argument
    Multiplies a row of a Matrix by a scalar, and returns it
    '''
    for a in m1:
        a[r1] *= sc
    return m1


def matscalmult(sc, m1):
    '''
    A Numerical Argument
    and a Second Order Composite List of Column Vectors as Argument
    Multiplies a scalar to every element of the Matrix
    FYI --> Scaling every vector doesn't effect the spectrum of Transformation
    '''
    for col in range(len(m1)):
        for row in range(len(m1[col])):
            m1[col][row] *= sc
    return m1


def addmat(m1, m2):
    '''
    Two Second Order Composite Lists of Column Vectors as Arguments
    Returns the Sum of Matrices
    '''
    if (len(m1) != len(m2)) or (len(m1[0]) != len(m2[0])):
        print("The Matrices should be of same dimensions, number of rows and columns should be equal respectively")
        return None
    else:
        cols = len(m1)
        rows = len(m1[0])
        nmat = []
        for a in range(cols):
            ncol = []
            for b in range(rows):
                ncol.append(m1[a][b] + m2[a][b])
            nmat.append(ncol)
        return nmat


def rowswap(m1,r1,r2):
    '''
    A Second Order Composite List of Column Vectors Argument
    Two Row no. Arguments (count from zero)
    Swaps two Rows of a Matrix
    .. and returns it
    '''
    for a in m1:
        k = a[r2]
        a[r2] = a[r1]
        a[r1] = k
    return m1


def rowmultadd(m1,r1,r2,sc):
    '''
    A Second Order Composite List of Column Vectors
    Two Row no. Arguments (count from zero)
    and a Numerical Argument Needed
    basically r1 ---> r1 + sc * r2
    Adds a scalar multiple of the r2-th row to the r1-th row (count from zero)
    .. and returns it
    '''
    for a in m1:
        a[r1] += sc*a[r2]
    return m1


def idmat(n):
    '''
    An Integer Argument
    Generates an Identity Matrix of n by n type for the input Integer 'n'
    '''
    m1 = []
    for a in range(n):
        col = []
        for b in range(n):
            if a == b:
                col.append(1)
            else:
                col.append(0)
        m1.append(col)
    return m1


def zeromat(n):
    '''
    An Integer Argument
    Generates a Zero Matrix of n by n type for the input Integer 'n'
    '''
    m1 = []
    for a in range(n):
        col = []
        for b in range(n):
            col.append(0)
        m1.append(col)
    return m1


def cofactors(m1):
    '''
    One Second Order Composite List of Column Vectors Argument
    Makes a Matrix of Cofactors a.k.a Chessboard Negation
    '''
    s = -1
    c = 0
    for a in range(len(m1)):
        for b in range(len(m1[0])):
            m1[a][b] = m1[a][b]*((s)**(c))
            c += 1
        if len(m1)%2 == 0:
            c += 1
    return m1


def adj(m1):
    '''
    One Second Order Composite List of Column Vectors Argument
    Returns the adjoint of a Matrix
    '''
    n = len(m1)-1
    for a in range(1,n+1):
        for b in range(int(a/2)+1):
            k = m1[b][a-b]
            m1[b][a-b] = m1[a-b][b]
            m1[a-b][b] = k
            if a != n:
                l = m1[n-b][n-(a-b)]
                m1[n-b][n-(a-b)] = m1[n-(a-b)][n-b]
                m1[n-(a-b)][n-b] = l
    return m1


def inverse(m1):
    '''
    Calculates the Inverse of a Linear Transformation
    Argument --> One Second Order Composite List of Column Vectors (Square Matrix)
    Returns --> Inverse of the Matrix in the same List form
    '''
    global qivar
    qivar = True
    
    if len(m1) != len(m1[0]):
        print("\n\n")
        print("ERROR: Inverse of a Non - Square Matrix Doesn't Exist")
        return None
    elif len(m1) == 2 and len(m1) == 2:
        res22 = matscalmult(1/det(m1), [[m1[1][1], -m1[0][1]],[-m1[1][0], m1[0][0]]])
        return res22
    
    detinv = det(m1)
    if detinv == 0:
        print("Inverse of a determinant 0 matrix is not defined, because the rank is being changed after transformation")
        return None

    print("\n")
    final = gauss_inv(m1)
    if triangcheck(coremult(final, m1), tol = 1e-6, idcheck = True):
        print("The Inverse by Gauss-Jordan Elimination is: ")
        return final
    else:
        detinv = det(m1)
        if detinv == 0:
            print("Inverse of a determinant 0 matrix is not defined, because the rank is being changed after transformation")
            return None
        final = matscalmult(1/detinv, adj(cofactors(minor(m1))))
        print("The Inverse by Martin Method / Cramer Rule is: ")
        return final


def triangcheck(m1, tol = 1e-12, idcheck = False):
    '''
    Checks if entered matrix is triangular or not
    --> values below tolerance "tol" will be considered zero
    Argument --> A second order composite list of Column Vectors
    Returns --> True or False
    '''
    lowertri = True
    uppertri = True
    diagisunit = True
    for col in range(len(m1)):
        for row in range(col):
            if m1[col][row] > tol or m1[col][row] < - tol: # --> equivalent to saying if |element| > tolerance
                lowertri = False

    for col in range(len(m1)):
        for row in range(len(m1)-1, col, -1):
            if m1[col][row] > tol or m1[col][row] < - tol:
                uppertri = False
    
    for icol in range(len(m1)):
        if m1[icol][icol] > (1 + tol) or m1[icol][icol] < 1 - tol:
            diagisunit = False

    if lowertri and uppertri and idcheck and diagisunit:
        return True
    if ((not lowertri) or (not uppertri) or (not diagisunit)) and idcheck:
        return False
    elif lowertri or uppertri:
        srsly = True
    else:
        srsly = False
    return srsly
   

def diagminor(m1):
    '''
    Argument --> A Square Matrix in the form of a Composite list of Column Vectors
    Returns --> The sum of determinants of minors along the diagonal
    '''
    s = 0
    for index in range(len(m1)):
        s += det(coreminor(m1,[index,index]))
    return s


def transpose(m1):
    '''
    Argument --> A Matrix in the form of a Composite list of Column Vectors
    Returns --> The Transpose of the Matrix in the same list form
    '''
    newmat = []
    for row in range(len(m1[0])):
        newcol = []
        for col in range(len(m1)):
            newcol.append(m1[col][row])
        newmat.append(newcol)
    return newmat


def matrixpadleft(m1):
    '''
    Pads the Matrix to the up and left so that-
    'A' becomes
    _________
    | 1 0 --|
    | 0 A --|
    | | \   |
    | |  \  |
    ---------

    Argument --> A Matrix in the form of a Composite list of Column Vectors
    Returns --> A Padded Matrix of bigger order
    '''
    for col in range(len(m1)): 
        m1[col] = [0] + m1[col]
    e2 = [1] + ([0]*(len(m1[0])-1))
    #print(e2)
    m1 = [e2] + m1
    return m1


def gsprojection(vec1, vec2):
    '''
    Arguments;
    --> vec1 : The vector which is to be projected
    --> vec2 : The vector on which vec1 is to be projected
    Returns --> The projection [<vec1, vec2>/<vec1, vec1>] * vec1
    '''
    try:
        _ = vec1[1]    #errors out if the vector is not a 1st order list
        _ = vec2[1]    #to check if the vector is in the form of a matrix
        vec1 = [vec1]
        vec2 = [vec2]
    except:
        pass
    #print(vec1, "dot", vec2)
    num = coremult(transpose(vec1), vec2)[0][0]
    #print(vec1, "dot", vec1)
    den = coremult(transpose(vec1), vec1)[0][0]
    sc = num/den
    return matscalmult(sc, vec1)


def roundmat(m1, roundto = 3):
    '''
    Rounds off the given Matrix elements to a certain decimal place, defaulted to 3
    Argument --> A Matrix in the form of a Composite list of Column Vectors
    Returns --> Rounded off Matrix in the same list Form
    '''
    m2 = []
    for icol in range(len(m1)):
        m2c = []
        for irow in range(len(m1[0])):
            m2c.append(round(m1[icol][irow], roundto))
        m2.append(m2c)
    return m2


def idvec(n):
    '''
    Returns the vector of ones [1, 1, ...] till n
    '''
    vec1 = []
    for a in range(n):
        vec1.append(0.5)
    return [vec1]


def polynomial(lofcfs, value):
    """
    Takes the sorted list of coefficients of a polynomial in the form --> 
    coefficient of x^n, coefficient x^n-1, ... up until the constant form
    Also takes the value at which the polynomial is to be evaluated
    Returns the thus evaluated value
    """
    val = 0
    n = len(lofcfs)
    for icf in range(n):
        val += lofcfs[icf] * (value ** (n-1-icf))
    return val


def eigenvalues(m1, acc = 3000):
    '''
    Argument --> A Second Order Composite List of Column Vectors
    Gives the eigenvalues as Output
    '''
    
    lofeigs = []
    tolerance = (1/acc) * 1e-6
    if len(m1) == len(m1[0]):
        #finding characteristic polynomial by
        #Faddeev-LeVerrier Method
        n = len(m1)
        coeff = 1
        auxmat = zeromat(n)
        lofcoeffs = [coeff]
        for k in range(1, n+1):
            auxmat = addmat(coremult(m1,auxmat), matscalmult(coeff,idmat(n)))
            coeff = (-1/k)*trace(coremult(m1,auxmat))
            lofcoeffs.append(coeff)
        polystring = ''
        for coeffind in range(len(lofcoeffs)):
            if lofcoeffs[coeffind] >= 0:
                polystring += '+ {}x{} '.format(lofcoeffs[coeffind],len(lofcoeffs)-coeffind-1)
            else:
                polystring += ' {}x{} '.format(lofcoeffs[coeffind],len(lofcoeffs)-coeffind-1)
        #print("The Characteristic Equation is --> ", polystring)
        if len(lofcoeffs) == 3:
            quadsolve = solveqn(lofcoeffs)
            #print("eigenvalues: ", quadsolve)
            lofeigs = quadsolve
    
    if triangcheck(m1, tol = tolerance) and (not lofeigs):
        lofeig = []
        for index in range(len(m1)):
            lofeig.append(m1[index][index])
    elif not lofeigs:
        
        # --- QR Factorization Method for eigenvalues for above 2 by 2
        
        m1a = []
        for icol in range(len(m1)):
            vec = []
            for irow in range(len(m1[0])):
                vec.append(m1[icol][irow])
            m1a.append(vec)

        c = 0
        sorryshaktimaan = False
        noqr = False
        while not triangcheck(m1a, tol = tolerance):
            ##############################################################
            qrtup = qrhousie(m1a)   # --> Factorize Matrix m1a to Q and R   ########## Prone to Divide by zero error
            ##############################################################
            if qrtup:
                q = qrtup[0]
                r = qrtup[1]
            else:
                print("Eigenvalue didn't converge, complex vectors maybe")
                noqr = True # --> QR Factorization failure Flag
                break
            m1a = coremult(r, q)   # --> Multiply RQ and then find its QR again till its triangular...
            #print("rq at the end of {}thiteration: ".format(c+1))
            #printmat(m1a)
            c += 1
            if c > acc or (not m1a):
                sorryshaktimaan = True   # --> QR Failure Flag due to complex values (almost triangular, but no)
                break

        if not sorryshaktimaan and not noqr:
            for index in range(len(m1a)):
                lofeigs.append(m1a[index][index])  # --> Extracting Diagonal Elements of Triangular Matrix
        elif sorryshaktimaan and not noqr:
            complexcase = None
            for index in range(len(m1a)):
                if index < len(m1a[0])-1 and (mod(m1a[index][index+1]) > tolerance):
                    complexcase = index
                    #print("complexcase is: ", complexcase)
                    #print("old lofeigs : ", lofeigs )
                    lofeigs += (eigenvalues([[ m1a[index][index], m1a[index][index+1] ], [m1a[index+1][index], m1a[index+1][index+1] ] ]))
                    #print("new lofeigs: ", lofeigs)
                elif index-1 != complexcase:
                    lofeigs.append(m1a[index][index])
        else:
            return None
    #try:
        #printmat(m1a)
    #except:
        #pass

    if lofeigs == []:
        print("Sorry We couldn't find the eigenvalues, try solving the characteristic equation")
        print(polystring)
    else:
        #print(polystring)
        return lofeigs

def qrhousie(m1a):
    '''
    QR Decomposition using Householder reflections
    Argument --> A Matrix in the form of a Composite list of Column Vectors
    Returns --> A Tuple of Q and R Matrices or None if things go sideways
    '''
    m2z = []
    for icol in range(len(m1a)):
        vec1b = []
        for irow in range(len(m1a[0])):
            vec1b.append(m1a[icol][irow])
        m2z.append(vec1b)
    
    m1b = []
    for icol in range(len(m1a)):
        vec2b = []
        for irow in range(len(m1a[0])):
            vec2b.append(m1a[icol][irow])
        m1b.append(vec2b)
    #print("m1b is: ", m1b)
    m1list = []
    q = idmat(len(m1b[0]))
    for xind in range(len(m1a)):
        x = m1b[xind][xind:]
        x = [x]
        #print("x vector and column starting index is: ", x, xind)
        e = [[1] + ([0]*(len(x[0])-1))]
        #print("e is: ", e)
        v = addmat(x, matscalmult((-1)**xind, matscalmult(v_norm(x), e)))
        #print("v is: ", v)
        try:
            ################################################
            c = 2/(coremult(transpose(v), v)[0][0])         ############ Prone to Divide by Zero error
            ################################################
        except:
            return qrhousie2(m2z)
        #print("c is: ", c)
        h = addmat(idmat(len(v[0])), matscalmult(-1, matscalmult(c, coremult(v, transpose(v)))))
        h1 = []
        for vec in range(len(h)):
            h1.append(h[vec])
        #print("Initial Householder: ")
        #printmat(h1)
        while len(h1) != len(m1b[0]):   # Why didn't I just use the matrixpadleft function?
            for vec in range(len(h1)):
                h1[vec] = [0] + h1[vec]
            frontvec = [1] + ([0]*(len(h1[0])-1))
            frontvec = [frontvec]
            h1 = frontvec + h1
        #print("Padded Householder")
        #printmat(h1)
        q = coremult(q, h1)
        #print("Matrix from which m1c is being carved: ")
        #printmat(m1b)
        m1c = []
        for vec in range(xind, len(m1b)):
            m1c.append(m1b[vec][xind:])
        #print("m1c is: ")
        #printmat(m1c)
        m1b = coremult(h, m1c)
        #print("Transformed half-Vector")
        #printmat(m1b)
        if len(m1b) != len(m1a):
            m1d = []
            for vec in m1list[-1]:
                m1d.append(vec)
            for i in range(-1, -len(m1b)-1, -1):
                #print("static element:", m1d[i][:(len(m1d[0]) - len(m1b[0]))])
                #print("Updating vector", m1b[i])
                m1d[i] = m1d[i][:(len(m1d[0]) - len(m1b[0]))] + m1b[i] # updating the matrix
                #print("modifying m1d: ")
                #printmat(m1d)
            m1b = m1d
        #print("m1b for next iteration is: ")
        #printmat(m1b)
        m1list.append(m1b)
    #print("Q: ")
    #printmat(q)
    #print("R: ")
    #printmat(m1b)
    #print("QR: ")
    qr = coremult(q, m1b)
    #printmat(qr)
    #print("rounded original matrix: ")
    #printmat(roundmat(m2z))
    #print("rounded qr: ")
    #printmat(roundmat(qr))
    if triangcheck(m1b, tol = 1e-5) and roundmat(qr) == roundmat(m2z):
        return (q, m1b)
    else:
        return qrhousie2(m2z)

def qrhousie2(m1a):
    '''
    Backup Function for qrhousie() with a slightly different approach
    QR Decomposition using Householder reflections
    Argument --> A Matrix in the form of a Composite list of Column Vectors
    Returns --> A Tuple of Q and R Matrices or None if things go sideways
    '''
    # Use as a backup function for qrhousie
    m2y = []
    
    for vec in m1a:
        nvec = []
        for el in vec:
            nvec.append(el)
        m2y.append(nvec)
    
    qa = idmat(len(m1a))
    qlist = []
    nohousie = False
    # t = min(len(m1a[0]), len(m1a))
    # counter = 0
    while (not triangcheck(qa)) or (qlist == []): # counter < t:
        x = m1a[0]
        #print("x is: ", x)
        e1 = [[1] + ([0]*(len(m1a[0])-1))]
        #print("e1 is: ", e1)
        x = [x]
        u = addmat(x, matscalmult(-1, matscalmult(v_norm(x), e1)))
        #print(u)
        v = matscalmult(1/v_norm(u), u)
        #print("v is: ", v)
        q = addmat(idmat(len(m1a)), matscalmult(-2, coremult(v, transpose(v))))
        #print("q is: ", q)
        if not (len(q) == len(qa)): # and len(q[0]) == len(qa[0])):
            q = matrixpadleft(q)
            #print("q has been changed to: ", q)
        qlist.append(q)
        #print("m1 to be multiplied to is: ", m1a)
        if len(qlist) == 1:
            qa = coremult(q, m1a)
        else:
            qa = coremult(q, qa)
        #print("q{}...q1a is: ".format(len(qlist)), qa)
        m1a = coreminor(qa, [0,0])
        #print("'A' for the next iteration is: ", m1a)
        # counter += 1
        if len(qlist) > 1000:
            nohousie = True
            break
    if not nohousie:
        q = transpose(qlist[-1])
        for qmat in range(len(qlist)-2, -1, -1):
            q = coremult(transpose(qlist[qmat]), q)
        #print("\nQ:\n")
        #printmat(q)
        r = qa
        #print("\nR:\n")
        #printmat(r)
        #print("\nRounded QR:\n")
        qr = coremult(q, r)
        #printmat(roundmat(qr))
        #print("\nRounded Original Matrix")
        #print(m2y)
        #printmat(roundmat(m2y))
        if triangcheck(r, tol = 1e-5) and roundmat(qr) == roundmat(m2y):
            #print("Triangular checked and rounded mats are same")
            return (q, r)
        else:
            print("\n\n--\\_(^^)_/--\n\n")
            return None
    else:
        print("\n\n--\\_(^^)_/--\n\n")
        return None
        
        

def printmat(lotups, neatness = 3):
    """
    Displays the things, takes in two positional parameters:
    ==> lotups, which is the general list of rows e.g. --> [('sn', 'name', 'job'), (1, 'john', 'assassin')] form
    """
    lotups = transpose(lotups)
    # making a list of longest words in a column for later use
    grlist = []
    for icol in range(len(lotups[0])):
        
        # searching for longest word in each column
        gr8 = lotups[0][icol]
        for irow in range(len(lotups)):
            if len(str(lotups[irow][icol])) > len(str(gr8)):
                gr8 = lotups[irow][icol]
        
        # extra whitespaces for neat-ness
        grlist.append(len(str(gr8)) + neatness)
    
    # actual printing process row by row
    
    print('\n', '-' * (sum(grlist) + neatness * (len(grlist) - 1)) , sep = '')
    
    for row in range(len(lotups)):
        for icol in range(len(lotups[row])):
            print(lotups[row][icol], end = ' ' * ( grlist[icol] - len(str(lotups[row][icol])) ) )
        print()
    
    print('-' * (sum(grlist) + neatness * (len(grlist) - 1)) , '\n', sep = '')

def run():
    '''
    User Menu function
    '''
    print("""

\    /_ | _ _ ._ _  _  _|_ _  |\/| _._|_ _ ._  _
 \/\/(/_|(_(_)| | |(/_  |_(_) |  |(_| |_(_)|_)_>
                                           |

    Report bugs at aryan_sidhwani@protonmail.com
    Solve some general problems in linear algebra here
    Leave the calculations to the computer !
""")
    default = []
    while True:
        print("""
    Things you can try out:
    1. solveqn --> Display roots a polynomial of any degree (!!!)
    2. det --> Displaying the Determinant of a Matrix
    3. inverse --> Converting the Matrix to it's Inverse (!)
    4. eigenvalues --> Displaying eigenvalues of the matrix (!!!)
    5. multiply -->  PreMultiplying the matrix to some other Matrix (!)
    6. powermat --> Exponentiating a Matrix (!)
    7. minor --> Converting the Matrix to it's Matrix of Minors (!)
    8. cofactors --> Chessboard - Negating the given Matrix (!)
    9. transpose --> Transposing the given Matrix (!)
    10. trace --> Displaying the trace of a matrix
    11. matscalmult --> Multiplying each element by a scalar (number) (!)
    12. rowmultadd -->  Adding scalar multiple of a row to another row (!)
    13. rowscalmult --> Multiplying a scalar to a row (!)
    14. rowswap --> Swapping two rows (!)
    15. diagminor --> Displays the sum of minors along a diagonal
    16. matrixpadleft --> Pads the matrix up and left (!)
    17. descartes_ros --> Displays the number of roots of a polynomial
    18. qrhousie --> Displaying the QR Factorization of the Matrix
    19. roundmat --> Rounding off the Matrix to a certain decimal place (!)
    
    --> (!) : This operation changes the default matrix you are working with
    --> (!!!) : Error prone but cool functions
        """)
        inp = input("Choose an option number [or type \'exit\' to exit]: ")
        if default == [] and inp != "exit":
            noread = False
            inpread = input("Do You want to fill in the default matrix you want to work on? [Y/n] ")
            if inpread == "n" or inpread == "no" or inpread == "No" or inpread == "N":
                noread = True
            if not noread:
                inpreadl = input("Do You want enter a python list of column vectors as Matrix? [Y/n] ")
                if inpreadl == "y" or inpreadl == "yes" or inpreadl == "Yes" or inpreadl == "Y":
                    default = eval(input("Enter a python list of lists: "))
                else:
                    default = read()
        if inp == "1":
            inpreadl = input("Do You want enter a python list of coefficients (of nth power to constant)? [Y/n] ")
            if inpreadl == "y" or inpreadl == "yes" or inpreadl == "Yes" or inpreadl == "Y":
                coeffs = eval(input("Enter a python list of decimals or integers: "))
            else:
                coeffs = []
                i = 0
                while True:
                    inp1 = input("coefficient of x ^ {}: [type \'exit\' to exit] ".format(i))
                    i += 1
                    if inp1 == "exit":
                        break
                    try:
                        inp1 = float(inp1)
                        coeffs.insert(0, inp1)
                    except:
                        print("Invalid input, repeating this one")
                        i -= 1

            sols = solveqn(coeffs)
            if sols:
                print("\nRoots: \n")
                for sol in sols:
                    if "complex" in str(type(sol)) and sol.imag >= 0:
                        print("--> {0} + {1} i".format(sol.real, sol.imag))
                    elif "complex" in str(type(sol)):
                        print("--> {0} - {1} i".format(sol.real, -sol.imag))
                    else:
                        print("--> {}".format(sol))
                print("\nPython List: ", sols)
            else:
                print ("No Solutions Found")

        elif inp == "2":
            print("\nThe Determinant is: \n", det(default))
        elif inp == "3":
            print("\nOriginal Matrix: ")
            printmat(default)
            default = inverse(default)
        elif inp == "4":
            try:
                egvs = eigenvalues(default)
                print("Eigenvalues are: ")
                for egv in egvs:
                    if "complex" in str(type(egv)) and egv.imag >= 0:
                        print("--> {0} + {1} i".format(egv.real, egv.imag))
                    elif "complex" in str(type(egv)):
                        print("--> {0} - {1} i".format(egv.real, -egv.imag))
                    else:
                        print("--> {}".format(egv))
                print("\nPython List: ", egvs)
            except:
                print("\nDivide by zero Error")
        elif inp == "5":
            print("\nEnter the Matrix to which the default matrix is to be PreMultiplied: ")
            inpreadl = input("Do You want enter a python list of Column Vectors as Matrix? [Y/n] ")
            if inpreadl == "y" or inpreadl == "yes" or inpreadl == "Yes" or inpreadl == "Y":
                tbm = eval(input("Enter a python list of decimals or integers: "))
            else:
                tbm = read()
            print("\nOriginal Matrix: ")
            printmat(default)
            default = multiply(default, tbm)
        elif inp == "6":
            while True:
                inp6 = input("\nEnter the integer to exponentiate to: [type \'exit\' to exit]")
                try:
                    inp6 = int(inp6)
                    break
                except:
                    print("\nInvalid input, try again")
            print("\nOriginal Matrix: ")
            printmat(default)
            default = powermat(default, inp6)
        elif inp == "7":
            print("\nOriginal Matrix: ")
            printmat(default)
            default = minor(default)
        elif inp == "8":
            print("\nOriginal Matrix: ")
            printmat(default)
            default = cofactors(default)
        elif inp == "9":
            print("\nOriginal Matrix: ")
            printmat(default)
            default = transpose(default)
        elif inp == "10":
            print("\nThe trace is: ", trace(default))
        elif inp == "11":
            while True:
                inp11 = input("\nEnter the scalar to multiply: ")
                try:
                    inp11 = float(inp11)
                    break
                except:
                    print("\nInvalid input, try again")
            print("\nOriginal Matrix: ")
            printmat(default)
            default = matscalmult(inp11, default)
        elif inp == "12":
            while True:
                print("r2 ---> r2 + Lr1")
                inp12a = input("\nEnter the row number of the to-be-changed (r2): ")
                inp12b = input("Enter the scalar (L): ")
                inp12c = input("Enter the row number of the to-be-scaled (temporarily) (r1): ")
                try:
                    inp12a = int(inp12a)-1   # --- minus 1 Because its a natural number, acc to end user
                    inp12c = int(inp12c)-1
                    inp12b = float(inp12b)
                    break
                except:
                    print("\nInvalid input, try again")
            print("\nOriginal Matrix: ")
            printmat(default)
            default = rowmultadd(default, inp12a, inp12c, inp12b)
        elif inp == "13":
            while True:
                print("r1 ---> Lr1")
                inp13a = input("\nEnter the row number of the to-be-scaled (r1): ")
                inp13b = input("Enter the scalar (L): ")
                try:
                    inp13a = int(inp13a)-1
                    inp13b = float(inp13b)
                    break
                except:
                    print("Invalid input, try again")
            print("\nOriginal Matrix: ")
            printmat(default)
            default = rowscalmult(default, inp13a, inp13b)
        elif inp == "14":
            while True:
                print("r2 <---> r1")
                inp14a = input("\nEnter the row number of r2: ")
                inp14b = input("Enter the row number of r1: ")
                try:
                    inp14a = int(inp14a)-1
                    inp14b = int(inp14b)-1
                    break
                except:
                    print("\nInvalid input, try again")
            print("\nOriginal Matrix: ")
            printmat(default)
            default = rowswap(default, inp14a, inp14b)
        elif inp == "15":
            print("\nSum of minors along diagonal is: ", diagminor(default))
        elif inp == "16":
            print("\nOriginal Matrix: ")
            printmat(default)
            default = matrixpadleft(default)
        elif inp == "17":
            inpreadl = input("Do You want enter a python list of coefficients (of nth power to constant)? [Y/n] ")
            if inpreadl == "y" or inpreadl == "yes" or inpreadl == "Yes" or inpreadl == "Y":
                coeffs = eval(input("Enter a python list of decimals or integers: "))
            else:
                coeffs = []
                i = 0
                while True:
                    inp17 = input("\ncoefficient of x ^ {}: [type \'exit\' to exit]".format(i))
                    i += 1
                    if inp17 == "exit":
                        break
                    try:
                        inp17 = float(inp17)
                        coeffs.append(inp17)
                    except:
                        print("\nInvalid input, repeating this one")
                        i -= 1
                print("\nNumber of solutions: ", descartes_ros(coeffs))
        elif inp == "18":
            q,r = qrhousie(default)
            print("\nQ:")
            printmat(q)
            print("As Python List of column Vectors: ", q)
            print("R:")
            printmat(r)
            print("As Python List of column Vectors: ", r)
        elif inp == "19":
            while True:
                inp19 = input("\nEnter the number of decimals you want to round to: ")
                try:
                    inp19 = int(inp19)
                    break
                except:
                    print("\nInvalid input, try again")
            default = roundmat(default)
            printmat(default)
        elif inp == "exit":
            if default != []:
                print("\nDefault Matrix: ")
                printmat(default)
                print("\nPython List as list of column vectors: \n", default)
                print("\nPython List as list of rows: \n", transpose(default))
            break
        else:
            print("\nInvalid Input, Try Again")

        
        if default != []:
            print("\nDefault Matrix right now: ")
            printmat(default)
            print("As Python List of column Vectors: ", default)


#mytup = qrhousie([[200, 1312, 1451], [721, 1329, 2745], [41455, 671, 1039]])
#for mat in mytup:
#    printmat(mat)
#print(mytup, type(mytup))
#run()
#x = input("Successful End of Program")

"""
tbd:
    eigenvectors from eigenvalues
    factorial of floats using gamma function
"""
