#include<iostream>
#include<vector>
#include<cmath>
#include<limits>

using namespace std;

#define EPSILON 1e-6

static double fval(double x){
    // return (x - cbrt(x) - 2);
    // return (pow(x,3) - 2 * pow(x,2) - 5);
    // return (pow(x,2) - 4*x + 4 - log(x));
    // return (exp(x) + pow(2,-1*x) + 2 * cos(x) - 6);
    // return (exp(-1*x) - x);
    // return (pow(x,3) + 4*pow(x,2) -10);
    // return (exp(x-2) + pow(x,3) - x);
    // return (x-pow(2,-1*x));

    // return (pow(x,1.4) - sqrt(x) + 1.0/x - 100);
    return (cos((pow(x,2)+5)/(pow(x,4)+1)));
    // return (1.0/ (100 + sqrt(x) - pow(x,1.4)));
}

/** get the sign(+.1 or -.-1 =.0) of any given number */
template <typename T>
int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

/*
    f(0)f(1) < 0 -> root within (0,1)
*/

// double bisection(double min , double max){
//     double mid = 0.0;
//     double t_min = min;
//     double t_max = max;
//     while(t_max-t_min>1e-6){
//         mid = (t_max + t_min) / 2;
//         double val = mid - cbrt(mid) - 2;
//         double lval = t_min-cbrt(t_min)-2;
//         // printf("%lf\n",val);
//         if(val*lval < 0){
//             t_max = mid;
//         }else{
//             t_min = mid;
//         }
//     }
//     return t_min;
// }

double bisection(double lpos,double rpos){
    double mid;
    double mval;
    while(rpos - lpos > EPSILON && lpos < rpos){
        mid = (lpos + rpos) / 2.0;
        mval = fval(mid);
        std::cout << "\nz: " << mval << "\t[" << lpos << " , " << rpos
                  << " | Bisect: " << mid << "]\n";
        if(mval < 0){
            lpos = mid;
        }else{
            rpos = mid;
        }
    }
    return mid;
}

int main(){
    double a = 0, b = 1, x, z;
    // int i;

    // for(int i = 0;i<50000;i++){
    //     z = fval(a);
    //     x = fval(b);
    //     if(sgn(z) == sgn(x)){ // same signs, increase interval
    //         b++;
    //         a--;
    //     }else{
    //         break;
    //     }
    // }

    cout << "First initial: " << a << endl;
    cout << "Second initial: " << b << endl;
    cout << endl;
    printf("Roots: %.4lf\n",bisection(a,b));
    // printf("%lf\n",0.6412 - pow(2,-0.6412));
    // printf("%lf\n",pow(27.8235,1.4)-sqrt(27.8235)+1.0/27.8235-100);
    printf("%lf\n",cos((pow(-0.6123,2)+5)/(pow(-0.6123,4)+1)));
    return 0;
}