#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    Newton’s method:
    -> xk+1 = xk - F(xk) / F'(xk)

    The method is linearly convergent if limi->∞ ei+1 / ei = S < 1
    The method is quadratically convergent if limi->∞ ei+1 / ei^2 = M (a constant)
                    F(x1)
    x2 = x1 - (---------------)
                   F'(x1)
*/

/*
    Secant method:
    -> xk+1 = xk - ( F(xk)*(xk - xk-1) ) / ( F(xk) - F(xk-1) )
    -> Need two initial guesses x0 and x1
        The secant method is superlinearly convergent ,
        meaning that it lies between linearly and quadratically
        convergent methods.
                               (x1 - x2) 
        x2 = x1 - f(x1) * ( ----------------- )
                             f(x1) - f(x2)
*/

double arr[2][100] = {{0}};
// Horner's rule
// 計算原方程和微分後方程的值
double cal(int type, double x, int n){ // Horner's rule
    double ans=0.0;
    for(int i=n;i>=0;i--){
        ans = (ans*x) + arr[type][i];
    }
    return ans;
}

double newton(double x, int n){
    double x1 = 0.0;
    while(1){
        x1 = x - (cal(0,x,n) / cal(1,x,n));
        if(fabs(x1-x)<=0.0001) break;
        x = x1;
    }
    return x1;
}

double f(double x){
    return (pow(x, 5.0) + x - 1);
}

double secant(double a, double b){
    while(abs(a - b) > 0.000001) {
        b = a - (f(a) * (a - b)) / (f(a) - f(b));
        swap(b, a);
    }
    cout << "ans: " << a << endl;
    return a;
}

int main(){
    int n = 0;
    cout << "please input the degree of eq: ";
    cin >> n;
    arr[1][n] = 0;
    for(int i = n;i >= 0 ;i--){
        cout << "Enter the cofficient of x^" << i << ": ";
        cin >> arr[0][i];
        arr[1][i-1] = arr[0][i] * i;
    }

    double x = 0.0;
    cout << "Input the value of x: ";
    cin >> x;
    double root = newton(x,n);
    printf("Roots is %.4lf\n",root);
    // cout << "Roots is " << root << endl;
    return 0;
}

// int main() {
//     double xi, xi1;

//     xi = 10, xi1 = 9;
//     while(abs(xi - xi1) > 0.00001) {
//         xi1 = xi - (func1(xi) * (xi - xi1)) / (func1(xi) - func1(xi1));
//         swap(xi1, xi);
//     }
//     cout<<"Q1: "<<xi<<endl;


//     xi = 10, xi1 = 9;
//     while(abs(xi - xi1) > 0.00001) {
//         xi1 = xi - (func2(xi) * (xi - xi1)) / (func2(xi) - func2(xi1));
//         swap(xi1, xi);
//     }
//     cout<<"Q2: "<<xi<<endl;

//     xi = 10, xi1 = 9;
//     while(abs(xi - xi1) > 0.00001) {
//         xi1 = xi - (func3(xi) * (xi - xi1)) / (func3(xi) - func3(xi1));
//         swap(xi1, xi);
//     }
//     cout<<"Q3: "<<xi<<endl;
//     return 0;
// }