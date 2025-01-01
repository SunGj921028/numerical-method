#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    fixed point iteration:
        f(x) = 0
        -> g(x) = x
        -> xk+1 = g(xk) for k = 0,1,2,3...

    we should select a g(x) that is easy to compute and 
    converges to the root of f(x)
*/

/* 
    The following four fixed point iteration functions methods are proposed t ocompute 21^(1/3)
    1. g(x) = (20x + (21/x^2)) / 21 -> speed of convergence = 6/7
    2. g(x) = x - (x^3 - 21) / (3x^2) -> speed of convergence = 0
    3. g(x) = x - ((x^4 - 21x) / (x^2 - 21)) -> speed of convergence = 5.~
    4. g(x) = (21/x)^(1/2) -> speed of convergence = 1/2
    Write down their speed of convergence and rank them based on their speed of convergence, assumeing that x0 = 1.

    rank = 2 > 4 > 1 > 3
*/

/*
    Example:
        x = 1-x^3 -> x = g(x) = 1-x^3 -> not converge
        x = 1-x^(2/3) -> x = g(x) = 1-x^(2/3) -> converge
        x = (1+2x^3) / (1+3x^2) -> x = g(x) = (1+2x^3) / (1+3x^2) -> converge

        convergence -> |g'(x)| < 1
*/

/*
    Let f be a continuous function on the interval [a ,b].
    There exists a number c between a and b such that 
    f(c) = ( f(b) - f(a) ) / (b - a)

    Let xi denote the iterate at step i.
    There exists a number c between xi and r such that
    g'(x) = ( g(xi) - g(r) ) / (xi - r)
          = ( xi+1 - r ) / (xi - r)

    (xi+1 - r) = g'(x) * (xi - r)
    e(i+1) = g'(c) * e(i)
*/

/*
    Bisextion v.s. Fixed point iteration:
    Which one is faster? -> Depends on the S = |g'(c)| is smaller or bigger than 1/2
    *Faster -> Fixed point iteration (S < 1/2)
*/

// e^x + x = 7 -> x = log(7 - x)
// x^1.4 - x^(1/2) + 1/x - 100 = 0 -> x = (x^(1/2) - 1/x + 100)^(1/1.4)

int main(){
    // f(x) = x^3 = 2x + 2 -> x = g(x) = cbrt(2x+2)
    double a = 0.0;
    double x = 0.0;
    cin >> x;
    // a = cbrt((2 * x) + 2);
    // a = (exp((-1) * x));
    // a = (sqrt(3 - log(x))); // x^2 + ln(x) = 3
    a = pow((sqrt(x) - 1.0/x + 100),1.0/1.4);
    // a = log(7 - x); // e^x + x = 7
    while(fabs(a-x) > 1e-10){
        x = a;
        // a = cbrt((2 * a) + 2);
        // a = (exp(-1 * x));
        // a = (sqrt(3 - log(x)));
        a = pow((sqrt(x) - 1.0/x + 100),1/1.4);
        // a = log(7 - x);
    }
    printf("ans: %.4lf\n",a);
    // cout << "ans = " << a << endl;
}