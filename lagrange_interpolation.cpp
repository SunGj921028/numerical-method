#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

/*
    1. Interpolation is the reverse of evaluation
    2. Interpolation is the process of finding a polynomial that passes through a given set of points
    3. Why do we need interpolation?
        -> we want informations as if the data were continuous(not discrete)
*/

/*
    Interpolation:
        The new value x ≠ xi is inside the range of the interpolation points x0 , x1 , …, xn
    Extrapolation:
        The new value z is outside the range
    Approximation:
        Some norm ||v - y|| of the difference of the two vectors is minimized
*/

/*
    Given n data points
    -> it have finite number of polynomials that interpolate the data points
    -> How many polynomials that interpolate the points? -> 1
    -> How many polynomials of degree n-1 or less that interpolate the points? -> 1
    -> Is P(X) unique? -> Yes
*/

/*
    Lagrange interpolation:
        P(x) = y1 * L1(x) + y2 * L2(x) + ... + yn * Ln(x)

        Lk(x) = (x - x1)(x - x2)...(x - xk-1)(x - xk+1)...(x - xn)
            --------------------------------------------------------- 
                (xk - x1)(xk - x2)...(xk - xk-1)(xk - xk+1)...(xk - xn)
    *Degree of P(x) = at most (n - 1)
        Ex:
            (0,1) (2,2) (3,4)
            P(x) = 1 * L1(x) + 2 * L2(x) + 4 * L3(x)
            L1(x) = (x - 2)(x - 3) / (0 - 2)(0 - 3) = (x - 2)(x - 3) / 6
            L2(x) = (x - 0)(x - 3) / (2 - 0)(2 - 3) = (x - 0)(x - 3) / -2
            L3(x) = (x - 0)(x - 2) / (3 - 0)(3 - 2) = (x - 0)(x - 2) / 3
*/

/*
    A degree d polynomial can have at most d roots, 
    unless it is the identically zero polynomial.
*/

/*
    Please estimate f(1.09)
    x =    {1.00, 1.05, 1.09, 1.10, 1.15}
    f(x) = {0.1924, 0.2414, ?, 0.2933, 0.3492}
*/

int main(){
    // n -> number of point we have
    int n = 11;
    // cin >> n;
    // vector<double> x{1.00,1.05,1.10,1.15};
    // vector<double> y{0.1924, 0.2414, 0.2933, 0.3492};
    // vector<double> x{1800,1850,1900,2000};
    // vector<double> y{280,283,291,370};
    // vector<double> x{8.1, 8.3, 8.6, 8.7};
    // vector<double> y{16.94410, 17.56492, 18.50515, 18.82091};
    vector<double> x{-5,-4,-3,-2,-1,0,1,2,3,4,5};
    vector<double> y{5,5,5,5,5,5,5,5,5,5,42};
    // init guess
    double ans_x = 0.0;
    // the point we want to estimate f(x)
    // double init_x = 1.09;
    double init_x = 6;
    // double init_x = 8.4;

    for(int i = 0;i < n; i++){
        double Up = 1.0;
        double Low = 1.0;
        for(int j = 0;j < n; j++){
            if(j != i){
                Up *= (init_x - x[j]);
                Low *= (x[i] - x[j]);
            }
        }
        ans_x += ( y[i] * Up / Low );
    }
    printf("ans = %lf\n",ans_x);
    return 0;
}