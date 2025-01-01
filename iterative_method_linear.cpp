#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    Jacobi method:
        Ax = b;
        (L+D+U)x = b;
        Dx = b - (L+U)x;
        xk+1 = D^-1(b - (L+U)xk); for k = 0,1,2...
    -> xk+1 = D^-1 * (b - (L + U)xk)
    -> D: diagonal matrix
    -> L: lower triangular matrix
    -> U: upper triangular matrix
    -> D^-1: inverse matrix of D
    -> b: constant vector
    -> xk: initial guess

    Gauss-Seidel method:
    -> xk+1 = (D - L)^-1 * (b - Uxk)
    -> D: diagonal matrix
    -> L: lower triangular matrix
    -> U: upper triangular matrix
    -> D^-1: inverse matrix of D
    -> b: constant vector
    Uses the most *recently updated values of the
    knowns, even if the updating occurs in the
    current step!
*/

/*
    strictly diagonally dominant matrix:
    -> |aii| > sum(|aij|)  for (j != i)
    -> Jacobi method will converge if A is strictly diagonally dominant.
*/

/*
    Starting from initial guess x0 , generate iterates x1 ,
    x2, …, xk , hopefully converging to solution x.

    benefits:
        Can be faster if the input matrix is large.
        A good approximation to the solution is already known.
        The input matrix is sparse.
*/

bool check(vector<vector<double>> A, int n){
    for(int i = 0; i < n;i++){
        double dia_num = A[i][i];
        double sum_other = 0.0;
        for(int j = 0;j < n;j++){
            if(j!=i){
                sum_other += fabs(A[i][j]);
            }
        }
        if(fabs(dia_num) < sum_other){
            return false;
        }
    }
    return true;
}

vector<vector<double> > matrix_add(vector<vector<double> > &a, vector<vector<double> > &b) {
    vector<vector<double> > ans(a.size(), vector<double>(a.front().size(), 0));
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < a.front().size(); j++) {
            ans[i][j] = a[i][j] + b[i][j];
        }
    }
    return ans;
}

vector<vector<double> > matrix_min(vector<vector<double> > &a, vector<vector<double> > &b) {
    vector<vector<double> > ans(a.size(), vector<double>(a.front().size(), 0));
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < a.front().size(); j++) {
            ans[i][j] = a[i][j] - b[i][j];
        }
    }
    return ans;
}

vector<vector<double> > matrix_mul(vector<vector<double> > &a, vector<vector<double> > &b) {
    vector<vector<double> > ans(a.size(), vector<double>(a.front().size(), 0));
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < b.size(); j++) {
            for(int k = 0; k < a.front().size(); k++) {
                ans[i][j] += (a[i][k] * b[k][j]);
            }
        }
    }
    return ans;
}

int main(){
    int n = 0;
    cout << "input n: ";
    cin >> n;
    // n = 10;
    
    vector<vector<double > > U(n, vector<double>(n,0));
    vector<vector<double > > L(n, vector<double>(n,0));
    vector<vector<double > > b(n, vector<double>(1,0));
    vector<vector<double > > D(n, vector<double>(n,0));

    // double array[6][6] = {
    //     {3, -1, 0, 0, 0, 0.5},
    //     {-1, 3, -1, 0, 0.5, 0},
    //     {0 , -1, 3, -1, 0, 0},
    //     {0 , 0, -1, 3, -1, 0},
    //     {0 , 0.5, 0, -1, 3, -1},
    //     {0.5, 0, 0, 0, -1, 3}
    // };

    // input matrix
    vector<vector<double> > array(n, vector<double>(n, 0));
    cout << "input matrix:\n";
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            cin >> array[i][j];
        }
    }
    // for(int i = 0;i<n;i++){
    //     for(int j = 0;j<n;j++){
    //         if(i == j){
    //             array[i][j] = 3;
    //         }else if(fabs(i - j) == 1.0){
    //             array[i][j] = -1.0;
    //         }else if(fabs(i - j) > 1.0){
    //             array[i][j] = 0;
    //         }
    //     }
    // }
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            printf("%.1lf ",array[i][j]);
        }
        printf("\n");
    }
    printf("--------------------\n");
    if(check(array,n) == false){
        cout << "not strictly diagonally dominant matrix" << endl;
        return 0;
    }

    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            if(i == j){
                D[i][j] = (1.0/array[i][j]);
            }else if(i < j){
                U[i][j] = array[i][j];
            }else if(i > j){
                L[i][j] = array[i][j];
            }
        }
    }

    // for(int i = 0;i<n;i++){
    //     D[i][i] = 1.0/3;
    // }
    // for(int i = 0; i < n-1; i++) {
    //     U[i][i + 1] = -1.0;
    // }
    // U[0][n-1] = 0.5;
    // U[1][n-2] = 0.5;
    // for(int i = 1; i < n; i++) {
    //     L[i][i - 1] = -1.0;
    // }
    // L[n-1][0] = 0.5;
    // L[n-2][1] = 0.5;
    cout << "D^-1" << endl;
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            printf("%.2lf ",D[i][j]);
        }
        printf("\n");
    }
    cout << "U" << endl;
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            printf("%.2lf ",U[i][j]);
        }
        printf("\n");
    }
    cout << "L" << endl;
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            printf("%.2lf ",L[i][j]);
        }
        printf("\n");
    }
    cout << "b" << endl;
    for(int i = 0;i<n;i++){
        cin >> b[i][0];
    }
    vector<vector<double> > X(n, vector<double>(1, 0));
    for(int i = 0;i<1500;i++){
        vector<vector<double> > UL_plus = matrix_add(U, L);
        vector<vector<double> > UL_mult_X = matrix_mul(UL_plus,X);
        vector<vector<double> > b_sub_ULX = matrix_min(b,UL_mult_X);
        X = matrix_mul(D,b_sub_ULX);
        // find specify turn to print
        if(i <= 6){
            cout << i + 1 << " turn" << endl;
            for(int j = 0;j<n;j++){
                printf("x%d = %.4lf ",j+1,X[j][0]);
            }
            printf("\n");
        }
    }
    for(int i = 0;i<n;i++){
        printf("u%d = %.4lf\n",i+1,X[i][0]);
    }
}