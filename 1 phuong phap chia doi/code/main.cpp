#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

double f(double x){
    return (7.0/8)*3.1416+(1.0/2)*sin(x)-x/2;
}

bool check_DK(double a, double b);
void halving_method(double a, double b, double epsi);

int main() {
    double a, b, epsi;
    cout << "\nNhap vao khoang cach li nghiem: ";
    cout << "\na = ";
    cin >> a;
    cout << "\nb = ";
    cin >> b;
    cout << "\nNhap vao sai so epsilon: ";
    cin >> epsi;
    halving_method(a, b, epsi);
    return 0;
}

bool check_DK(double a, double b){
    if (f(a) * f(b) >= 0)
        return false;
    else
        return true;
}

void halving_method(double a, double b, double epsi) {
    fstream out;
    out.open("output.txt", ios::out);
    out << "Phương pháp chia đôi giải gần đúng nghiệm bài toán f(x) = 0";
    out << "\nInput: ";
    out << "\nPhương trình: x^3 - 4x^2 + x + 1 = 0";
    out << "\nKhoảng cách li nghiệm: (a, b) = (" << a << ", " << b << ")";
    out << "\nSai số epsilon: " << epsi << "\n";
    out << "\nChạy code theo công thức sai số hậu nghiệm \n";
    if (check_DK(a, b)) {
        int lan_lap = 0;
        double a1 = a;
        double b1 = b;
        double c, delta;
        do {
            lan_lap ++;
            out << "\nLần lặp thứ " << lan_lap << ": ";
            c = (a1 + b1) / 2;
            if (f(c) == 0)
                out << "Nghiệm x = " << fixed << setprecision(20) << c;
            else {
                if (f(a1) * f(c) < 0)
                    b1 = c;
                else
                    a1 = c;
            }
            delta = fabs(a1 - b1);
            out << "Nghiệm x =  " << fixed << setprecision(20) << c << ",  Sai số = " << delta;
        }
        while (delta > epsi);
        out << "\n\nKết luận: \nSố lần lặp: " << lan_lap << "\nNghiệm gần đúng: " << c;
    }
    else
        out << "\n Không sử dụng được phương pháp !";
}