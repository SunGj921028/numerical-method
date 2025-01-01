#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Multivariate Newton’s method
// F(x) is to store the value of the equation with initial guess
// Df(xk)s = -F(xk)
//   A   x =    b   -> use gassian elimination to solve x
/* Every round plus s to xk
    -> xk+1 = xk + s   for k = 0,1,2,3...
*/

/*
    Jacobian matrix:
    ->below is the example of 3 variables
    Df(x) = [偏微變數1 偏微變數2 偏微變數3]
            [偏微變數1 偏微變數2 偏微變數3]
            [偏微變數1 偏微變數2 偏微變數3]...
    
*/

// to find s (naive)
// vector<double> gass(vector<vector<double> > a, vector<double> b) {
//     int n = a.size();
//     vector<double> ans(n, 0);
//     for(int i = 0; i < n; i++) {
//         for(int j = i + 1; j < n; j++) {
//             double times = a[j][i] / a[i][i];
//             for(int k = i; k < n; k++) {
//                 a[j][k] -= a[i][k] * times;
//             }
//             b[j] -= b[i] * times;
//         }
//     }

//     for(int i = n - 1; i >= 0; i--) {
//         for(int j = i + 1; j < n; j++) {
//             b[i] -= a[i][j] * ans[j];
//         }
//         ans[i] = b[i] / a[i][i];
//     }
//     return ans;
// }

// (complete)
vector<double> gass(vector<vector<double> > a, vector<double> b) {
    int n = a.size(); // row of A
    vector<double> ans(n, 0);
    for(int i = 0; i < n; i++) {
        // sort
        double MAX = __DBL_MIN__;// __DBL_MIN__
        int index = i;
        for(int j = i; j < n; j++) {
            if(MAX < fabs(a[j][i])) { // row
                MAX = fabs(a[j][i]);
                index = j;
            }
        }
        for(int j = 0; j < n; j++) {
            swap(a[i][j], a[index][j]);
        }
        swap(b[i], b[index]);
        // operation
        for(int j = i + 1;j < n;j++){
            double t = a[j][i] / a[i][i];
            for(int m = i;m < n; m++){
                a[j][m] -= (a[i][m] * t);
            }
            b[j] -= b[i] * t;
        }
    }
    // calculate
    for(int i = n - 1;i>=0;i--){
        for(int j = i + 1;j<n;j++){
            b[i] -= a[i][j] * ans[j];
        }
        ans[i] = b[i] / a[i][i];
    }
    return ans;
}
// 3xy^2 -> 對x微分 = 3y^2
int main(){
    // number of variables
    int n = 2;
    // initial guess( size is the number of variables )
    double ans[2] = {5,10};
    // Jacobian matrix with variables
    vector<vector<double> > Df(n, vector<double>(n,0));
    // First value of the equation
    vector<double> f(n);

    for(int i = 0;i < 5000;i++){
        for(int j = 0;j < n;j++){
            if(j==0){
                Df[j][0] = 3*pow(ans[0],2) - 3*pow(ans[1],2);
                Df[j][1] = -6 * ans[0] * ans[1];
                // Df[j][2] = 1.0 / ans[2];
                // Df[j][0] = (4 * ans[0] - 4);
                // Df[j][1] = 2 * ans[1];
                // Df[j][2] = (6 * ans[2] + 6);
                f[j] = (-1) * (pow(ans[0],3) - 3*ans[0]*pow(ans[1],2) - 0.5);
                // f[j] = (-1) * (2*pow(ans[0],2) - 4*ans[0] + pow(ans[1],2) + 3*pow(ans[2],2) + 6*ans[2] + 2);
            }else if(j==1){
                // Df[j][0] = 3 *(ans[0]);
                // Df[j][1] = 2 * ans[1];
                Df[j][0] = 6*ans[0]*ans[1];
                Df[j][1] = 3*pow(ans[0],2) - 3*pow(ans[1],2);
                // Df[j][2] = (-3) * pow(ans[2],2);
                // Df[j][0] = 2*ans[0];
                // Df[j][1] = 2*ans[1] - 2;
                // Df[j][2] = 4*ans[2];
                f[j] = (-1) * (3*pow(ans[0],2)*ans[1] - pow(ans[1],3) - sqrt(3)/2.0);
                // f[j] = (-1) * (pow(ans[0],2) + pow(ans[1],2) - 2*ans[1] + 2*pow(ans[2],2) - 5);
            }
            // }else{
            //     Df[j][0] = 1;
            //     Df[j][1] = 1;
            //     Df[j][2] = 1;
            //     // Df[j][0] = 6*ans[0] - 12;
            //     // Df[j][1] = 2*ans[1];
            //     // Df[j][2] = 6*ans[2];
            //     f[j] = (-1) * (ans[0] + ans[1] + ans[2] - 5);
            //     // f[j] = (-1) * (3*pow(ans[0],2) - 12*ans[0] + pow(ans[1],2) + 3*pow(ans[2],2) + 8);
            // }
        }
        // s is the answer of Df(xk)s = -F(xk)
        vector<double> s;
        s = gass(Df,f);
        // update ans
        for(int i = 0;i<n;i++){
            ans[i] += s[i];
        }
    }
    // print the answer
    for(int i = 0; i < n; i++) {
        printf("%.4lf\n",ans[i]);
    }
    return 0;
}

/*
    Example:
        f1(x1,x2,x3) = sin(x) + y^2 + ln(z) - 7 = 0
        f2(x1,x2,x3) = 3x + 2y - z^3 + 1 = 0
        f3(x1,x2,x3) = x + y + z - 5 = 0
        initial guess: (0,2,2)T
        answer: (0.6331,2.3934,1.9735)
        Df(x) = [cos(x) 2y 1/z]
                [3     2  -3z^2]
                [1     1   1  ]
*/