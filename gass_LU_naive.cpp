#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    gassian elimination (naive):
    -> A: coefficient matrix
    -> b: constant vector
    -> x: solution vector
    -> n: size of matrix
    -> A * x = b

    LU Factorization:
    -> A = L * U
    -> L: lower triangular matrix
    -> U: upper triangular matrix
    -> (L * U) * x = b
        Define y = U * x
            -> (L * y = b) -> (y = L^-1 * b) for y
            -> (U * x = y) -> (x = U^-1 * y) for x
*/

/*
    Operation count of Gaussian elimination:
    -> O(2/3 * n^3) + O(n^2) -> O(n^3)
    Complexity of Gaussian elimination:
    Elimination: O(n^3)
    Back substitution: O(n^2)
*/

/*
    1.
    Swamping:
        ex: [0 1 3] -> 用第一列消除第二列時，應將a2j減去(a21/a11)*a1j，但此時a11為0，會出錯
            [2 5 1]
            [3 0 2]
        -> Solution: partial pivoting
            Focus on the absolute value of multiples t obe no longer larger than 1.
            -> |Iik| <= 1

    2.
    x = [u v w] , initial guess = [0 0 0]
    Df(x) = [2 -1 0]    F(x) = [2u-v]
            [0 2 -1]           [-u+2v-w-2]
            [-1 0 2]           [-v+2w]
    u = (v)/2 , v = (u+2+w)/2 , w = (v)/2
    Gass seidel method:
        [0 0 0] -> [0 1 1/2] -> [1/2 3/2 3/4] -> 拿上一個最近的結果
    Jacobi method:
        [0 0 0] -> [0 1 0] -> [1/2 1 1/2] -> 拿上一組的結果
    
    3.
    If newton's method fail -> use bisection method or select a new g(x) which |g'(x)| < 1, then again using FPI
*/

// Does a matrix have a LU factorization? -> Yes, if A is invertible.

// naive gassian elimination
vector<double> gass(vector<vector<double> > a, vector<double> b) {
    int n = a.size();
    vector<vector<double> > L(n, vector<double>(n, 0));
    vector<double> ans(n, 0);
    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            // diagonal (if no need L matrix, above for j = i + 1)
            if(i == j){
                L[i][j] = 1;
                continue;
            }
            double times = a[j][i] / a[i][i];
            L[j][i] = times;
            for(int k = i; k < n; k++) {
                a[j][k] -= a[i][k] * times;
            }
            b[j] -= b[i] * times;
        }
    }

    for(int i = n - 1; i >= 0; i--) {
        for(int j = i + 1; j < n; j++) {
            b[i] -= a[i][j] * ans[j];
        }
        ans[i] = b[i] / a[i][i];
    }
    cout << "Upper:\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j=0;j<n;j++){
            cout << a[i][j] << " ";
        }
        cout << "]";
        cout << endl;
    }
    printf("--------------------\n");
    cout << "Lower:\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j = 0;j<n;j++){
            cout << L[i][j] << " ";
        }
        cout << "]";
        cout << endl;
    }
    return ans;
}

int main(){
    int n = 3; // can be input
    // double a[3][3] = {
    //         {3,-1,2},
    //         {6,-1,5},
    //         {-9,7,3},
    //     };
    vector<vector<double > > a(n, vector<double>(n,0));
    cout << "input A matrix:\n";
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            cin >> a[i][j];
        }
    }
    // double b[3] = {10, 22, -7};
    vector<double > b(n,0);
    cout << "input B vector:\n";
    for(int i = 0;i<n;i++){
        cin >> b[i];
    }
    vector<vector<double > > L(n, vector<double>(n,0));

    cout << "A\n";
    for(int i = 0;i<n;i++){
        cout << "[";
        for(int j=0;j<n;j++){
            cout << a[i][j] << " ";
        }
        cout << "]";
        cout << endl;
    }
    printf("--------------------\n");
    vector<double > ans(n,0);
    ans = gass(a, b);
    cout << "ans:\n";
    cout << "x1 = " << ans[0] << endl;
    cout << "x2 = " << ans[1] << endl;
    cout << "x3 = " << ans[2] << endl;
    return 0;
}

/*
    Example:
        3x1 - x2 + 2x3 = 10
        6x1 - x2 + 5x3 = 22
        -9x1 + 7x2 + 3x3 = -7
        ans:
        x1 = 1
        x2 = -1
        x3 = 3
*/