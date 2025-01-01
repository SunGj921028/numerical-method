#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

/*
    x* = argmin f(x)
    x' = x - learning_rate * f'(x)

    learning rate -> step size
*/

/*
    Warning:
        Gradient descent is a greedy method!
        It may return a local minimum!
*/

double gradient_descent(double x) {
    // the polynomial is after first derivative to get the slope
    // return ( 4 * pow(x, 3) - ( 6 * pow(x, 2) ) );
    // return (4 * pow(x,3) + 9 * pow(x - 2, 2) - (30 * x));
    // return (4*pow(x,3) + 9*pow(x-2,2) - 30*x);
    return (3 * cos(x));
}

int main() {
    double x = -2;
    double learning_rate = 0.0001;
    // int count = 0;

    // gradient descent
    // find closest local minimum
    while(fabs(gradient_descent(x)) > 1e-6) {
        // count++;
        x -= (learning_rate * gradient_descent(x));
    }

    printf("x = %.4lf\n", x);
    // cout << count << endl;
    return 0;
}