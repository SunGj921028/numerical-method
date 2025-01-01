#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    The normal equation:
    -> (Ax)^T (b - A~x) = 0 (because b - A~x is orthogonal to the column space of A)
    -> x^T A^T (b - A~x) = 0
    -> A^T (b - A~x) = 0
    -> (A^T) (A~x) = (A^T)b
    -> (A^T)A(~x) = (A^T)b   (~x) -> least square solution

    Given an inconsistent system Ax = b -> A(m*n) x(n*1) = b(m*1)
        solve (A^T)A(n*n) (~x)(n*1) = (A^T)b(n*1) for the least square solution (~x)
        that minimizes the Euclidean lenth(norm) ||b - A~x|| of the residual vector r = b - A(~x). where A(~x) = (~b)
*/

/*
    Measure the residual size:
        1. Euclidean length(2-norm) ||r||2 =  ||b - A~x||2 = sqrt(r1^2 + r2^2 + ... + rm^2)
        2. squared error SE = r1^2 + r2^2 + ... + rm^2
        3. Root mean square error RMSE = sqrt(SE / m) -> m is the number of data points
*/

/*
    Periodic model:
        y = c1 + c2 * cos2PI(x) + c3 * sin2PI(x)

    **Expoential model:
        y = c1 * exp(c2 * t) -> c2 may not be linearly
        -> ln(y) = ln(c1) + c2 * t ln = loge
        -> let k = ln(c1) -> ln(y) = k + c2 * t
*/

/*
    *Model linearization changes the least square problem !
        The original problem:
            min (c1 * e^(c2 * t1) - ym)^2 + ... + (c1 * e^(c2 * tn) - ym)^2
        The linearized problem:
            min (ln(c1) + c2 * t1 - ln(ym))^2 + ... + (ln(c1) + c2 * tn - ln(ym))^2
                -> errors in "log space"
*/

/*
    Practice:
        - A list of 50 numbers of atmospheric carbon dioxide yi, record at year 1961 to 2010 ti. -> means 50 data points
        - Subtract the background level of 279 ppm
        - suntract the background year 1960

        Fit the data to an Exponential model (y - 279 = c1 * exp(c2 * (t - 1960)))
        Report c1,c2, and RMSE
*/

/*
    1 2 x 1 3
    3 4   2 4
*/

vector<double> gass(vector<vector<double> > a, vector<double> b){
    int n = a.size();
    vector<double> ans(n, 0);
    double B[10][10] = {{0}};

    vector<vector<double> > P(n, vector<double>(n , 0));
    // permutation matrix
    for(int i = 0; i < n; i++) {
        P[i][i] = 1;
    }

    // swap
    for(int i = 0; i < n; i++) {
        // printf("\n");
        // for(int l = 0;l < n;l++){
        //     for(int k = 0;k < n;k++){
        //         printf("%lf ",a[l][k]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");
        double MAX = __DBL_MIN__;
        int index = i;
        for(int j = i; j < n; j++) {
            if(MAX > abs(a[j][i])) {
                MAX = abs(a[j][i]);
                index = j;
            }
        }
        for(int j = 0; j < n; j++) {
            swap(a[i][j], a[index][j]);
            swap(P[i][j], P[index][j]);
            swap(B[i][j], B[index][j]);
        }
        swap(b[i], b[index]);
        for(int j = i;j < n;j++){
            if(i==j){
                B[i][j] = 1;
                continue;
            }
            double t = a[j][i] / a[i][i];
            B[j][i] = t;
            for(int m = i;m < n; m++){
                a[j][m] -= a[i][m] * t;
            }
            b[j] -= b[i] * t;
        }
    }

    for(int i = n - 1;i>=0;i--){
        for(int j = i + 1;j<n;j++){
            // a[i][n] -> b[i]
            b[i] -= a[i][j] * ans[j];
        }
        // no solution
        if(a[i][i] == 0) {
            cout << "No solution" << endl;
            exit(0);
        }
        ans[i] = b[i] / a[i][i];
    }
    return ans;
}

/*
    please find the least squares solution of x1 and x2 of the following equations:
        x1 - 0.4 = 0;
        x2 - 0.8 = 0;
        x1^2 + x2^2 - 1 = 0;

    -> A =  [1 0]
            [0 1]
            [2x1 2x2]
    -> b =  [0.4]
            [0.8]
            [1]
*/

int main(){
    // num -> the number of data points
    int num = 0;
    // cin >> n;
    // num = 50;
    num = 3;
    // cout << pow(0.4243,2) + pow(0.8482, 2) << endl;
    // cout << cbrt(16.0/25) << endl;
    // return 0;
    // for the first column of A
    // 1 - 50
    // vector<double > x;
    // for(int i = 0;i < num;i++){
    //     x.push_back(i + 1);
    // }

    // vector<double > x{0,1,2,3};

    // vector<double > b{10, 5, 2, 1};
    vector<double > b{0.4, 0.8, 1};
    // vector<double > b{10, 10, -5, 15, 0};
    // vector<double > b{3,5,8};
    // vector<double > b{
    //     320.58 ,
    //     321.01 ,
    //     322.25 ,
    //     322.24 ,
    //     322.16 ,
    //     324.01 ,
    //     325.00 ,
    //     325.57 ,
    //     327.34 ,
    //     328.07 ,
    //     328.92 ,
    //     330.07 ,
    //     332.48 ,
    //     333.09 ,
    //     333.97 ,
    //     334.87 ,
    //     336.75 ,
    //     338.01 ,
    //     339.47 ,
    //     341.46 ,
    //     342.91 ,
    //     344.14 ,
    //     345.75 ,
    //     347.43 ,
    //     348.93 ,
    //     350.21 ,
    //     351.84 ,
    //     354.22 ,
    //     355.67 ,
    //     357.16 ,
    //     359.34 ,
    //     359.66 ,
    //     360.28 ,
    //     361.68 ,
    //     363.79 ,
    //     365.41 ,
    //     366.80 ,
    //     369.30 ,
    //     371.00 ,
    //     371.82 ,
    //     374.02 ,
    //     375.55 ,
    //     378.35 ,
    //     380.61 ,
    //     382.24 ,
    //     384.94 ,
    //     386.43 ,
    //     388.49 ,
    //     390.18 ,
    //     393.22
    // };

    // for(int i = 0;i < n;i++){
    //     cin >> b[i];
    // }

    vector<double > b_origin(num);
    b_origin.assign(b.begin(), b.end());

    // the equation is ln(y) = ln(c1) + c2 * (t - 1960)
    // y = b - 279
    // x[i] = ti - 1960

    // for(int i = 0;i < num;i++){
    //     // b[i] -= 279;
    //     b[i] = log(b[i]);
    // }

    // 未知數數量 (c1, c2)
    int n = 0;
    n = 2;

    // A -> 50 * 2
    vector<vector<double> > A(b.size(), vector<double> (n));
    // vector<vector<double> > A{{4}, {7}, {11}};

    A = {
        {1, 0},
        {0, 1},
        {0.4, 0.8}
    };

    // A = {
    //     {3,-1,2},
    //     {4,1,0},
    //     {-3,2,1},
    //     {1,1,5},
    //     {-2,0,3}
    // };

    // if j has three column -> then the third column is x^2
    // for(int i = 0;i < num;i++){
    //     for(int j = 0;j < n;j++){
    //         if(j == 0){
    //             A[i][j] = 1;
    //         }else{
    //             A[i][j] = x[i];
    //         }
    //     }
    // }

    cout << "ATA:\n";
    vector<vector< double > > ATA(n, vector<double>(n, 0));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            double sum = 0;
            // for(int k = 0; k < b.size(); k++) {
            //     sum += ( A[k][i] * A[k][j] );
            // }
            for(int k = 0; k < b.size(); k++) {
                sum += ( A[k][i] * A[k][j] );
            }
            ATA[i][j] = sum;
            cout << ATA[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout  << "ATb: " << endl;
    vector<double > ATb(n, 0);
    for(int i = 0; i < n; i++) {
        double sum = 0;
        // for(int j = 0; j < b.size(); j++) {
        //     sum += ( A[j][i] * b[j] );
        // }
        for(int j = 0; j < b.size(); j++) {
            sum += ( A[j][i] * b[j] );
        }
        ATb[i] = sum;
        cout << ATb[i] << " ";
    }
    cout << endl;
    vector<double > ans = gass(ATA, ATb);
    // cout << ans[0] << endl;
    // cout << endl;
    // the ans[0] is ln(c1)
    // so c1 = exp(ans[0])
    // ans[0] = exp(ans[0]);
    for(int i = 0; i < n; i++) {
        cout<< "c" << i+1 << ": " << ans[i] << endl;
    }
    // for(int i = 0;i<ans.size();i++){
    //     cout<<ans[i]<<" ";
    // }
    cout<<endl;

    // bonus
    // vector<double > pre(num);
    // for(int i = 0;i < num;i++){
    //     pre[i] = 279 + (ans[0] * exp(ans[1] * x[i]));
    // }
    // vector<double > minus(50);
    // double sum = 0;
    // for(int i = 0;i<50;i++){
    //     minus[i] = pow(b_origin[i] - pre[i], 2);
    //     sum += minus[i];
    // }
    // cout << "RMSE: \n";
    // cout << sqrt(sum / 50) << endl;
    return 0;
}

// int main() {
//     vector<vector<double>> node = {{-1, 0, 1}, {1, 0.5, 0.5}, {1, -0.5, 0.5}, {0, 1, 0.5}};
//     int n = node.size(), m = node.front().size();
//     vector<vector<double>> Dr(n, vector<double>(m, 0));
//     vector<double> ans = {0, 0, 0};

//     for(int t = 0; t < 1000; t++) {
//         vector<double> rx(n, 0);
//         for(int i = 0; i < n; i++) {
//             Dr[i][0] = (ans[0] - node[i][0]) / sqrt( pow(ans[0] - node[i][0], 2) + pow(ans[1] - node[i][1], 2));
//             Dr[i][1] = (ans[1] - node[i][1]) / sqrt( pow(ans[0] - node[i][0], 2) + pow(ans[1] - node[i][1], 2));
//             Dr[i][2] = -1;
//             rx[i] = (-1) * (sqrt( pow(ans[0] - node[i][0], 2) + pow(ans[1] - node[i][1], 2) ) - (node[i][2] + ans[2]));
//         }
//         //display(Dr);
//         vector<vector<double>> DrT = transpose(Dr);
//         vector<vector<double>> DrTDr = mul_matrix(DrT, Dr);
//         vector<double> DrTrx = mul_vector(DrT, rx);
//         vector<double> v = guss(DrTDr, DrTrx);
//         for(int i = 0; i < m; i++) {
//             ans[i] += v[i];
//         }
//     }

//     for(int i = 0; i < m; i++) {
//         cout<<ans[i]<<" ";
//     }
//     cout<<endl;

//     return 0;
// }