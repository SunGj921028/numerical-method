#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

/*
    P(x) = 2x^4 + 3x^3 - 3x^2 + 5x - 1
    P'(x) = 8x^3 + 9x^2 - 6x + 5
    -> -1 + x(5 + x(-3 + x(3 + x(2))))

    P(x) = 4x^5 + 7x^8 - 3x^11 + 2x^14
    P'(x) = 20x^4 + 56x^7 - 33x^10 + 28x^13
    -> x^5(4 + x^3(7 + x^3(-3 + x^3(2))))
    -> addition time = 3
    -> multiplication time = 4
    -> total time = 7
*/

/*
    Please provide a solution to fix the divergence problem.
    -> 1. change the iteration function (find a better convergence g(x)), then use FPI method again.)
    -> 2. change the iteration method (bisection)
*/

int arr[2][100] = {{0}};
// Horner's rule
// 計算原方程和微分後方程的值
double cal(int type, double x, int n){ // Horner's rule
    double ans=0.0;
    for(int i=n;i>=0;i--){
        ans = (ans*x)+arr[type][i];
    }
    return ans;
}

int main(){
    int n=0;
    cout << "input max degree: \n";
    cin>>n;
    int **arr = new int*[2];
    int num = 0;
    cout << "how many equation: \n";
    cin >> num;
    cout << "input coefficient of equation from smallest power: \n";
    for(int i=0;i<num;i++){
        arr[i] = new int[n+1];
    }
    for(int i=0;i<num;i++){
        for(int j=0;j<n+1;j++){
            cin>>arr[i][j];
        }
    }
    double x=0.0;
    cout<<"input x: \n";
    cin>>x;
    double f=cal(0,x,n);
    cout<<f<<endl;
    if(num == 2){
        double df=cal(1,x,n);
        cout<<df<<endl;
    }

    // double ans = 0.0;
    // ans = 8 * pow(x, 5) - pow(x,4) - 3 * pow(x, 3) + pow(x, 2) - 3 * x + 1;
    // cout << "ans = " <<ans << endl;
    return 0;
}