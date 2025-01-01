#include<iostream>
#include<vector>

using namespace std;

/*
    newton's divided-difference interpolation:
        P(x) = c0 + c1(x - x1) + c2(x - x1)(x - x2) + ... + cn(x - x1)(x - x2)...(x - xn-1)
        P(x) = f[x1] + f[x1,x2](x - x1) + f[x1,x2,x3](x - x1)(x - x2) + ... + f[x1,x2,...,xn](x - x1)(x - x2)...(x - xn-1)
            denote that f[x1 x2 x3 ... xn] = cn-1
            0-th divided difference: f[xi] = P(xi) = yi
            1-st divided difference: f[xi,xi+1] = (f[xi+1] - f[xi]) / (xi+1 - xi)  -> (yi+1 - yi) / (xi+1 - xi)
            k-th divided difference: f[xi,xi+1,...,xi+k] = (f[xi+1,xi+2,...,xi+k] - f[xi,xi+1,...,xi+k-1]) / (xi+k - xi)

            Ex: f[x1 x2 x3 x4] = (f[x2 x3 x4] - f[x1 x2 x3]) / (x4 - x1)
*/

/*
    Triangular array:
        x1 | *f[x1]
                    *f[x1 x2] = f[x2] - f[x1] / x2 - x1
        x2 | f[x2]                                      *f[x1 x2 x3]
                    f[x2 x3] = f[x3] - f[x2] / x3 - x2              *f[x1 x2 x3 x4]
        x3 | f[x3]                                      f[x2 x3 x4]
                    f[x3 x4] = f[x4] - f[x3] / x4 - x3
        x4 | f[x4]
*/

/*
    interpolation error (importance exam question):
        interpolation error at x(the target point) is 
            f(x) - P(x) =  (x - x1)(x - x2)...(x - xn)
                          ---------------------------  * f^(n)(c)
                                      (n)!
            *where c is between the smallest and largest of xi
            *n -> number of point we have

        Ex:
            for point xi: -1 -0.5 0 0.5 1
            f(x) - P(x) = (x + 1)(x + 0.5)(x)(x - 0.5)(x - 1) / 5! * f^(5)(c)
                -> range of c: -1 <= c <= 1
                -> maximal value of |f^(5)(c)| on [-1,1] is e^1 = e
                -> |f(x) - P(x)| <= (x + 1)(x + 0.5)(x)(x - 0.5)(x - 1) / 5! * e
                *at x = 0.75, |f(x) - P(x)| <= 0.002323 = error
                *at x = 0.25, |f(x) - P(x)| <= 0.000995 = error
        The interpolation error will be smaller close
        to the center of the interpolation interval.
*/

/*
    Runge phenomenon:
        polynomial wiggle near the ends of the interpolation interval.
*/

double newton_divided_diffenence(vector<vector<double> > data, double x){
    // n -> number of point we have
    int n = data.size();
    vector<double> divided(n,0);
    // for ci(divided diffenence)
    vector<double> divided_val(n,0);

    // c0 = y1 = f(x1)
    divided_val[0] = data[0][1];
    for(int i = 0;i < n; i++){
        // divided = f(x1), f(x2), f(x3), f(x4), f(x5)
        divided[i] = data[i][1];
    }

    for(int i = 1; i < n; i++) {
        for(int j = 0; j + i < n; j++) {
            // triangle divided diffenence (array)
            divided[j] = (divided[j + 1] - divided[j]) / (data[j + i][0] - data[j][0]);
        }
        // only take the value of divided[0] (c1, c2, c3, c4)
        divided_val[i] = divided[0];
    }

    // print the ci
    for(int i = 0;i < n; i++){
        cout << "c" << i << ": " << divided_val[i] << endl;
    }

    // Similar horner's rule
    double ans1 = divided_val[n - 1];
    for(int i = n - 2; i >= 0; i--) {
        ans1 = (ans1 * (x - data[i][0])) + divided_val[i];
    }
    return ans1;
}

int main(){
    // vector<vector<double> > data{
    //     {1800.0, 280.0},
    //     {1850.0, 283.0},
    //     {1900.0, 291.0},
    //     {2000.0, 370.0},
    // };
    // vector<vector<double> > data{
    //     {8.1, 16.94410},
    //     {8.3, 17.56492},
    //     {8.6, 18.50515},
    //     {8.7, 18.82091},
    // };
    // vector<vector<double> > data{
    //     {-1, 1},
    //     {0, 1},
    //     {1, 2},
    //     {2, 0},
    // };
    // vector<vector<double> > data{
    //     {1, 0.7652},
    //     {1.3, 0.6201},
    //     {1.6, 0.4554},
    //     {1.9, 0.2818},
    //     {2.2, 0.1104},
    // };
    vector<vector<double> > data{
        {-1, 1},
        {0, 3},
        {1, 5},
        {2, 2},
        {3, 0.5},
    };
    // vector<vector<double> > data{
    //     {-5, 5},
    //     {-4, 5},
    //     {-3, 5},
    //     {-2, 5},
    //     {-1, 5},
    //     {0, 5},
    //     {1, 5},
    //     {2, 5},
    //     {3, 5},
    //     {4, 5},
    //     {5, 42},
    // };
    // the point we want to estimate f(x)
    // double x = 1950.0;
    // double x = 8.4;
    double x = 6;
    printf("f(%lf) = %lf\n",x,newton_divided_diffenence(data, x));
    return 0;
}