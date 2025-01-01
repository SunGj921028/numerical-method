#include<iostream>
#include<vector>
#include<cmath>

using namespace std;
// Multiples in Gaussian elimination should be kept as small as possible to avoid swamping.

/*
    partial pivoting:
    The index p of the pivot in the k th step of GE is determined by
    |apk^(k-1)| = max(|aik^(k-1)|) for i = k,k+1,...,n

    -> swap the rows to put the largest element in the pivot position.
    Forces the absolute value of multiples to be no larger than 1 (|Iik| <= 1).
*/

/*
    PALU decomposition (partial pivoting):
    -> Ax = b
    -> PAx = Pb
    -> LUx = Pb
    -> Ly = Pb -> (Ux = y)
    solve:
    -> Ly = Pb for y
    -> Ux = y for x
*/

/*
    matrix norm:
    -> ||A||inf = max(sum(|aij|)) for i = 1,2,3...n
    Example:
    -> A  = [1 2 3]
            [4 5 6]
            [7 8 9]
        -> ||A||inf = 24
    
    infinity norm:
    -> ||x||inf = max(|xi|) for i = 1,2,3...n

    backward error (small) -> use approximation:
    -> ||b - Ax||inf = ||r||inf

    forward error (large) -> use correct answer:
    -> ||x - xk||inf

    What does this mean (backward error (small) and forward error (large))?
    -> Even though the “approximate solution” is
    relatively far from the exact solution, it nearly
    lies on both lines (slope nearly same)!
*/

/*
    Find the relation between forward error (FE) and backward error (BE) for finding the root of the linear function f(x) = ax − b.
    BE = |f(xa)| = |axa − b| = |axa − ax| = |a(xa − x)| = |a|*|xa − x| = |a|*FE
    -> BE = |a|*FE
    FE = |(b/a) - xa|

    xi+1 = xi - f(xi) / f'(xi) = xi - (axi+b) / a
    Relation = (xi - (axi+b)/a - r) / (xi - r)
*/

/*
    error magnification factor:
        (relative forward error / relative backward error)
        -> (||x-xk||inf / ||x||inf) 
          ---------------------------- = ||A||inf * ||A^-1||inf = condition number of A
           (||b-Ax||inf / ||b||inf)
    Example:
        x = [1 1]T
        xa = [-1 3.0001]T
        b = [2 2.0001]T
        relative forward error = ||x-xa||inf / ||x||inf = 2.0001 / 1 = 2.0001
        relative backward error = ||b-Axa||inf / ||b||inf = 0.0001 / 2.0001 = 0.00005
        error magnification factor = 2.0001 / 0.00005 = 40004.0001

        error magnification factor = condition number of A = ||A||inf * ||A^-1||inf
*/

int main(){
    // augmented matrix
    // A is a n*n matrix
    int n = 4; // can be input
    // double a[4][5] = {
    //     {1.19,2.11,-100,1,1.12},
    //     {14.2,-0.122,12.2,-1,3.44},
    //     {0,100,-99.9,1,2.15},
    //     {15.3,0.11,-13.1,-1,4.16}
    //     }; //equation
    // double a[4][5] = {
    //     {1.19,2.11,-100,1,1.12},
    //     {14.2,-0.122,12.2,-1,3.44},
    //     {0,100,-99.9,1,2.15},
    //     {15.3,0.110,-13.1,-1,4.16}
    //     }; //equation
    // double a[4][5] = {
    //     {1,0,0,1,0},
    //     {-1,1,0,1,0},
    //     {-1,-1,1,1,0},
    //     {-1,-1,-1,1,0}
    //     }; //equation
    // double a[4][5] = {
    //     {4,2,-1,3,16.9},
    //     {3,-4,2,5,-14},
    //     {-2,6,-5,-2,25},
    //     {5,1,6,-3,9.4}
    //     }; //equation
    
    // double a[4][5] = {
    //     {M_PI,sqrt(2),-1,1,0},
    //     {M_E,-1,1,2,1},
    //     {1,1,-1*sqrt(3),1,2},
    //     {-1,-1,1,-1*sqrt(5),3}
    //     }; //equation

    double a[4][5] = {
        {1,-2,3,0,0},
        {1,-2,3,1,0},
        {1,-2,2,-2,0},
        {2,1,3,-1,0}
        }; //equation


    double b[10][10] = {{0}};

    vector<vector<double> > P(n, vector<double>(n + 1, 0));
    // permutation matrix
    for(int i = 0; i < n; i++) {
        P[i][i] = 1;
    }
    // swap
    for(int i = 0; i < n; i++) {
        printf("\n");
        for(int l = 0;l<4;l++){
            for(int k = 0;k<4;k++){
                printf("%lf ",a[l][k]);
            }
            printf("\n");
        }
        printf("\n");
        double MAX = __DBL_MIN__;
        int index = i;
        for(int j = i; j < n; j++) {
            if(MAX < abs(a[j][i])) {
                MAX = abs(a[j][i]);
                index = j;
            }
        }
        for(int j = 0; j < n + 1; j++) {
            swap(a[i][j], a[index][j]);
            swap(P[i][j], P[index][j]);
            swap(b[i][j], b[index][j]);
        }
        // swap(b[i], b[index]);
        for(int j = i;j<n;j++){
            if(i==j){
                b[i][j] = 1;
                continue;
            }
            double t = a[j][i] / a[i][i];
            b[j][i] = t;
            // printf("%lf %lf %.20lf\n",a[j][i],a[i][i],t);
            printf("%lf\n",t);
            // when m == n, it is the b vector
            for(int m = i;m <= n; m++){
                a[j][m] -= a[i][m] * t;
            }
        }
    }

    cout<<"P:"<<endl;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout<<P[i][j]<<" ";
        }
        cout<<endl;
    }
    cout << "-------------------\n";
    cout << "A\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j=0;j<n + 1;j++){
            cout << a[i][j] << " ";
        }
        cout << "]";
        cout << endl;
    }
    printf("--------------------\n");
    vector<double > ans(n,0);
    for(int i = n - 1;i>=0;i--){
        for(int j = i + 1;j<n;j++){
            // a[i][n] -> b[i]
            a[i][n] -= a[i][j] * ans[j];
        }
        // no solution
        if(a[i][i] == 0) {
            cout << "No solution" << endl;
            return 0;
        }
        ans[i] = a[i][n] / a[i][i];
    }
    cout << "Upper:\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j=0;j<n;j++){
            // cout << a[i][j] << " ";
            printf("%lf ",a[i][j]);
        }
        cout << "]";
        cout << endl;
    }
    cout << "Lower:\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j = 0;j<n;j++){
            // cout << b[i][j] << " ";
            printf("%lf ",b[i][j]);
        }
        cout << "]";
        cout << endl;
    }
    cout << "ans:\n";
    cout << "x1 = " << ans[0] << endl;
    cout << "x2 = " << ans[1] << endl;
    cout << "x3 = " << ans[2] << endl;
    // cout << "x4 = " << ans[3] << endl;
    return 0;
}