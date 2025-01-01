#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>

using namespace std;

/*
    Orthogonal matrix:
        The column vectors form an orthogonal set.
        --> Q^T = Q^-1
    
    QR factorization (avoid computing A^T * A):
        A = QR
        Q: orthogonal matrix -> Q^T * Q = I (Q^T = Q^-1)
        R: upper triangular matrix
        A: m * n
        Q: m * n
        R: n * n

        -> Ax = b
        -> QRx = b
        -> Q^-1 * QRx = Q^-1 * b = Q^T * b
        -> Rx = Q^T * b

    Gram-Schmidt process (Orthogonalization a set of vectors):
        input: v1, v2, ..., vn
        q1 = v1
        u2 = v2 - (v2 * q1) * q1
        q3 = v3 - (v3 * q1) * q1 - (v3 * q2) * q2
        ...
        qn = vn - (vn * q1) * q1 - (vn * q2) * q2 - ... - (vn * qn-1) * qn-1
        q1, q2, ..., qn -> orthogonal basis
        q1', q2', ..., qn' -> orthonormal basis

        Q = [q1', q2', ..., qn'] -> after normalization
        For R:
            rjj = ||qj||2
            rij = qj^T * vi

    From reduced to full QR:
        Add a arbitrary vector A3

    Summary:
        QR factorization -> recording the orthogonalization of a matrix
*/

/*
    Solve inconsistent system of linear equations:
        QR factorization:
            Ax = b
            QRx = b
            Q^-1 * QRx = Q^-1 * b = Q^T * b
            Rx = Q^T * b
*/

/*
    Conditioning of least squares:
        Recall that we compute cond(A) for error estimation on sloving (A^T)x = b:
            cond(A) = ||A|| * ||A^-1||
    The cond(A) is too large for least square
*/

void show_matrix(vector<vector<double> > A) {
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A.front().size(); j++) {
            // cout << setw(10) << to_string(A[i][j]) << " ";
            printf("%10.4f ", A[i][j]);
        }
        cout << endl;
    }
    return;
}

vector<vector<double> > transform(vector<vector<double> > A) {
    // change matrix to A^T
    // A.size() -> row
    // A.front().size() -> col == A[0].size()
    // change (m * n) to (n * m)
    vector<vector<double> > result(A.front().size(), vector<double>(A.size(), 0.0));
    for(int i = 0; i < A.front().size(); i++) {
        for(int j = 0; j < A.size(); j++) {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

vector<double> subtract_vector(vector<double> A, vector<double> B) {
    for(int i = 0; i < A.size(); i++) {
        A[i] -= B[i];
    }
    return A;
}

// 內積
double inner_mul_vector(vector<double> A, vector<double> B) {
    double sum = 0;
    for(int i = 0; i < A.size(); i++) {
        sum += (A[i] * B[i]);
    }
    return sum;
}

// 常數乘上向量
vector<double> mul_scalar_vector(vector<double> A, double scalar) {
    for(int i = 0; i < A.size(); i++) {
        A[i] *= scalar;
    }
    return A;
} 

// calculate the length of column vector
double normalize(vector<double> &v) {
    double sum = 0;
    for(int i = 0; i < v.size(); i++) {
        sum += (v[i] * v[i]);
    }
    sum = sqrt(sum);
    for(int i = 0; i < v.size(); i++) {
        double tmp = v[i] / sum;
        v[i] = tmp;
    }
    return sum;
}

void gram_schmidt(vector<vector<double> > A) {
    // A -> m * n
    int n = A.front().size();
    vector<vector<double> > Q;
    // Q -> n * m -> A^T
    // because A's v1 is the u1 of Q
    // so, we can use the column of A as the row of Q
    // then, Q[0] is the u1 of Q
    // Q[1] is the u2 of Q
    Q = transform(A);
    vector<vector<double> > R(n, vector<double>(n, 0.0));

    for(int i = 0;i < Q.size(); i++){
        for(int j = 0;j < i;j++){
            // the Q[j] before changing is the v2, v3, v4, ...
            // we calculate the inner product of u1 . v2
            double inner = inner_mul_vector(Q[j], Q[i]);
            Q[i] = subtract_vector(Q[i], mul_scalar_vector(Q[j], inner)); // inner is Q1^T * v2 and Q[j] is u1, Q[i] is u2
            R[j][i] = inner;// place the inner product in the R[j][i] because R is the upper triangular matrix
        }
        R[i][i] = normalize(Q[i]); // place the length
    }
    Q = transform(Q); // we need to transform the Q to the original form as (m * n) matrix

    cout<<endl<<"Q: "<<endl;
    show_matrix(Q);
    cout<<endl<<"R: "<<endl;
    show_matrix(R);
}

int main(){
    int m = 0;
    int n = 0;
    // vector<vector<double> > A;
    // for(int i = 0;i < m;i++){
    //     for(int j = 0 ; j< n;j++){
    //         cin>>A[i][j];
    //     }
    // }
    // vector<vector<double> > A{
    //     {4, 2, 3, 0},
    //     {-2, 3, -1, 1},
    //     {1, 3, -4, 2},
    //     {1, 0, 1, -1},
    //     {3, 1, 3, -2}
    // };
    // vector<vector<double> > A{
    //     {2, -1, 0, -1},
    //     {2, 2, 1, 0},
    //     {0, 0, 0, 0},
    //     {2, -1, 1, 1},
    //     {1, 1, 0, 1}
    // };
    // vector<vector<double> > A{
    //     {1, -2, -1},
    //     {2, 0, 1},
    //     {2, -4, 2},
    //     {4, 0, 0}
    // };
    vector<vector<double> > A{
        {2, 2, 2, 3},
        {3, -5, -2, 3},
        {2, -2, 5, -3},
        {-1, -5, -5, 0},
        {2, -4, -1, -1},
        {-3, 3, -1, 1}
    };
    // vector<vector<double> > A{
    //     {1, -4},
    //     {2, 3},
    //     {2, 2}
    // };
    gram_schmidt(A);

    return 0;
}
