///Thuật toán chính chạy code newton tiếp tuyến

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double f(double x) {
    
}

double df(double x) {
    double h = 1e-7;
    return ((f(x + h) - f(x - h)) / (2 * h));
}

double ddf(double x) {
    float h = 1e-3;
    return ((df(x + h) - df(x - h)) / (2 * h));
}
///min ham so
double min(double f(double x), double a, double b){
    double x0 = a, alpha;
    double min = a;
    alpha = (b - a) / 10000;
    for (int i = 0; i < 10000; i++){
        if (f(min) > f(x0))
            min = x0;
        x0 = x0 + alpha;
    }
    return f(min);
}
///max ham so
double max(double f(double x), double a, double b){
    double x0 = a, alpha;
    double max = a;
    alpha = (b - a) / 10000;
    for (int i = 0; i < 10000; i++){
        if (f(max) < f(x0))
            max = x0;
        x0 = x0 + alpha;
    }
    return f(max);
}
//Sign of a number
double sign(double x)
{
    if (x>=0) return 1;
    else return -1;
}

double f1(double x){
    return fabs(df(x));
}
double f2(double x){
    return fabs(ddf(x));
}
///Su dung cong thuc sai so muc tieu
//Newton Raphson method
double sol(double a, double b, double e) {
    fstream output;
    output.open("output.txt", ios::out);
    double s, x0, m1, m2, c;
    int i = 1;
    if (sign(f(a)) == sign(ddf(a))) x0 = a;
    else x0 = b;

    m1 = min(f1, a, b);
    m2 = max(f2, a, b);
    c = sqrt(2.0 * m1 * e / m2);

    output << "\nDiem Fourier la: " << x0;
    do {
        s = x0;
        x0 = x0 - f(x0) / df(x0);
        output <<"\nLan thu: " << i;
        output <<" x = " << setprecision(20) << x0;
        output <<" sai so = " << setprecision(20) << fabs(x0 -s);
        output <<" fx = " << setprecision(20) << f(x0);

        i++;
    }
    while (fabs(x0 - s) >= c);

    if (fabs(x0 - s) < c) {
        output << "\nVay so lan lap: " << i - 1;
        output << " nghiem x = " << x0;
    }
    output.close();
    return (x0);
}

int main() {
    double a, b, e;
    //Verify that this interval is "good"
    do {
        printf("\n Nhap a: ");
        scanf("%lf", &a);
        printf("\n Nhap b: ");
        scanf("%lf", &b);

        if (sign(f(a)) == sign(f(b))) {
            printf("Day khong phai khoang cach ly nghiem! Hay nhap lai!");
        }

        if (a == b) {
            printf("a phai khac b! Hay nhap lai!");
        }
    }
    while ((sign(f(a)) == sign(f(b))) || (a == b));

    if (a > b) {
        a = a + b;
        b = a - b;
        a = a - b;
    }

    printf("\n Nhap gia tri sai so: ");
    scanf("%lf", &e);

    //Check that f' and f'' preserve sign on this interval
    if (sign(max(ddf, a, b)) != sign(min(ddf, a, b))) {
        printf("Dao ham cap hai doi dau tren doan nay. Thuat toan khong thuc hien duoc");
    }
    else {
        if (sign(max(df, a, b)) != sign(min(df, a, b))) {
            printf("Dao ham cap mot doi dau tren doan nay. Thuat toan khong thuc hien duoc");
        }
        else {
            //Find root in the interval
            sol(a, b, e);
        }
    }
    return 0;
}