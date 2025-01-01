#include<iostream>

/*
    Please find the line ( y = ax + b ) that fits the following four points (1,6), (2, 5), (3, 7), (4, 10).
    -> y = 1.4x + 3.5

    Please find the least squares solution of k of the following equations.
    4k = 3, 7k = 5, 11k = 8
    -> least squares solution of k = 0.725806

    Let P(x) be the interpolating polynomial for the data (0, 0), (0.5, y), (1, 3), and (2, 2). The coefficient of x3 in P(x) is 6. Find y.
    -> y = 13/4
    P(x) = 0*L1(x) + y*L2(x) + 3*L3(x) + 2*L4(x)
    L2(x) = (x - 0)(x - 1)(x - 2) / (0.5 - 0)(0.5 - 1)(0.5 - 2)
    L3(x) = (x - 0)(x - 0.5)(x - 2) / (1 - 0)(1 - 0.5)(1 - 2)
    L4(x) = (x - 0)(x - 0.5)(x - 1) / (2 - 0)(2 - 0.5)(2 - 1)
*/

/*
    Suppose we are given a set of n data points (t1, y1), … (tn, yn).
    We choose the power law model 𝑦=𝑐1𝑡^𝑐2 to fit the data points.
    Below are the steps to compute c1, c2 and the least squares solution 𝑦̃𝑖 using linearization. 
    Fill the ‘?’ parts in steps 3 and 4.

    1. Apply the natural logarithm to the model: ln 𝑦=ln(𝑐1𝑡^𝑐2)=ln𝑐1+𝑐2ln𝑡
    2. Let 𝑘=ln𝑐1, and therefore, ln 𝑦=𝑘+𝑐2ln𝑡.
    3. Construct a matric A=[?]
    4. Construct a vector b=[?]
    5. Let 𝑥=[𝑘 𝑐2]^T .
        Solve the system of normal equations 𝐴𝑇𝐴𝑥=𝐴𝑇𝑏 . Compute 𝑐1=𝑒^𝑘.
    6. Compute the least squares solution 𝑦̃𝑖=𝑐1𝑡^𝑐2.

    -> A =  [1 ln(t1)]
            [1 ln(t2)]
            [1 ln(t3)]
            [1 ln(tn)]
    -> b =  [ln(y1)]
            [ln(y2)]
            [ln(y3)]
            [ln(yn)]
*/

/*
    We choose the exponential model 𝑦=𝑐1∙𝑒^𝑐2∙𝑡 to fit the data points.
    1. Apply the natural logarithm to the model: ln 𝑦=ln(𝑐1∙𝑒^𝑐2∙𝑡)=ln𝑐1+𝑐2∙𝑡
    2. Let 𝑘=ln𝑐1, and therefore, ln 𝑦=𝑘+𝑐2∙𝑡.
    3. Construct a matric A=[?]
    4. Construct a vector b=[?]
    5. Let 𝑥=[?]. Solve the system of normal equations 𝐴𝑇𝐴𝑥=𝐴𝑇𝑏. Get 𝑐1=𝑒^𝑘

    -> A =  [1 t1]
            [1 t2]
            [1 t3]
            [1 tn]
    -> b =  [ln(y1)]
            [ln(y2)]
            [ln(y3)]
            [ln(yn)]
    -> x =  [k]
            [c2]
*/

/*
    The three-part problem is stated as follows: 
    Divide 10 into three parts such that they shall be 
    in continued proportion to each other and the product 
    of the first and the last two shall be 6. 
    Taking x, y, and z as three parts, this problem can be represented
    as a system as follows:
        x + y + z = 10
        x / y = y / z
        x * y = 6
        y * z = 6

    -> Dr() = [1      1         1    ]
              [1/y (-x/y^2 - 1/z) y/z^2]
              [y      x         0    ]
              [0      z         y    ]
    -> r() = [x + y + z - 10]
             [x/y - y/z]
             [xy - 6]
             [yz - 6]
*/

/*
    Assume that the polynomial P9(x) interpolates the function f(x) = e^(-2x)
    at the 10 evenly spaced points x = 0, 1/9, 2/9, …, 1.
    Find an upper bound for the error |f(1/2) - P9(1/2)|.
    
    By newton's divided difference
    -> f(x) - P(x) = (x - x1)(x - x2)...(x - xn) / n! * f^(n)(c)
    -> f(x) = e^(-2x)
    then f^(n)(c) = (-2)^n * e^(-2c)
    -> |f(1/2) - P9(1/2)| <= (1/2)(1/2 - 1/9)(1/2 - 2/9)...(1/2 - 1) / 10! * (-2)^10 * e^(-2c)
*/

/*
    Given (-1,1) (0,-1) (1,-1)
    (a) Lagrange interpolation
    (b) Newton's interpolation

    -> P(x) = 1 * L1(x) + (-1) * L2(x) + (-1) * L3(x)
    -> L1(x) = (x - 0)(x - 1) / (-1 - 0)(-1 - 1) = (x - 0)(x - 1) / 2
    -> L2(x) = (x + 1)(x - 1) / (0 + 1)(0 - 1) = (x + 1)(x - 1) / -1
    -> L3(x) = (x + 1)(x - 0) / (1 + 1)(1 - 0) = (x + 1)(x - 0) / 2
    -> P(x) = x^2 - x -1

    -> P(x) = f[x1] + f[x1,x2](x - x1) + f[x1,x2,x3](x - x1)(x - x2)
            = 1 + (x + 1) * (-1 - 1) / (0 - (-1)) + (x + 1)(x - 0) * ( [-1-(-1)/(1-0)] - [(-1-1)/(0-(-1))] ) / (1 - (-1))
            = x^2 - x - 1
*/

/*
    Please find the matrix Dr needed for applying  Gauss-Newton iteration
    to the model-fitting problem with three data points (t1,y1) (t2,y2) (t3,y3)
    (a) translated exponential y = c3 + c1 * exp(c2*t)
    (b) power law y = c1 * t^c2

    -> y1 = c3 + c1 * exp(c2 * t1)
    -> y2 = c3 + c1 * exp(c2 * t2)
    -> y3 = c3 + c1 * exp(c2 * t3)
    Dr() =  [exp(c2*t1)  t1 * c1 * exp(c2*t1)  1]
            [exp(c2*t2)  t2 * c1 * exp(c2*t2)  1]
            [exp(c2*t3)  t3 * c1 * exp(c2*t3)  1]

    -> y1 = c1 * t1^c2
    -> y2 = c1 * t2^c2
    -> y3 = c1 * t3^c2
    Dr() =  [t1^c2  c1 * t1^c2 * ln(t1)]
            [t2^c2  c1 * t2^c2 * ln(t2)]
            [t3^c2  c1 * t3^c2 * ln(t3)]
*/

/*
    
*/

/*
    Cubic spline:
        A cubic spline has the property that its second derivative is zero 
        at the endpoints of the interval.
*/
