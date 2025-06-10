#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;
#define N 100
#define max_loop 1000000
#define input "input.txt"

void inputDaThuc(double coeff_poly[N], int &n, double &epsi);
void output(double coeff_poly[N], int n);

double derivative(double x, double coeff_poly[N], int n);
double function_value(double x, double coeff_poly[N], int n);
void domain_solution(double &lower_bround, double &upper_bround, double coeff_poly[N], int n);
double halving_method(double x0, double x1, double coeff_poly[N], int n, double epsi);
void solve_polynomial(double coeff_poly[N], int n, double epsi);

int main(){
    double coeff_poly[N], epsi;
    int n;
    inputDaThuc(coeff_poly, n, epsi);
    output(coeff_poly, n);
    cout << "\nSai so la: " << epsi;
    cout << "\nGiai phuong trinh da thuc: ";
    cout << "\nDung phuong phap Gradient Descent de tim cuc tri. \nDung phuong phap chia doi de tim nghiem\n";
    solve_polynomial(coeff_poly, n, epsi);
    cout << "\n";
    system("pause");
    return 0;
}

void inputDaThuc(double coeff_poly[N], int &n, double &epsi){
    fstream in;
    in.open(input, ios::in);
    int c1;
    double c2;
    in >> c1;
    n = c1;
    for (int i = 0; i <= n; i++){
        in >> c2;
        coeff_poly[i] = c2;
    }
    in >> c2;
    epsi = c2;
    in.close();
}
void output(double coeff_poly[N], int n){
    cout << "\n\nPhuong trinh da thuc vua nhap la: ";
    cout << coeff_poly[0] << "x^" << n;
    for (int i = 1; i<= n; i++)
        cout << " + " << coeff_poly[i] << "x^" << n - i;
}

double derivative(double x, double coeff_poly[N], int n) {
    double temp = 0;
    for (int i = n; i >= 1; i--) {
        temp += coeff_poly[n - i] * i * pow(x, i - 1);
    }
    return temp;
}
double function_value(double x, double coeff_poly[N], int n) {
    double temp = 0;
    for (int i = n; i >= 0; i--) {
        temp += coeff_poly[n - i] * pow(x, i);
    }
    return temp;
}
void domain_solution(double &lower_bround, double &upper_bround, double coeff_poly[N], int n) {
    double temp[N], max, k;
    upper_bround = -1;
    lower_bround = 1;
    for (int i = 0; i <= n; i++) {
        temp[i] = coeff_poly[i];
    }
    do {
        max = 0;
        k = 0;
        if (temp[0] < 0) {
            for (int i = 0; i <= n; i++) {
                temp[i] = -temp[i];
            }
        }
        for (int i = 1; i <= n; i++) {
            if (temp[i] < 0) {
                k = i;
                break;
            }
        }
        for (int i = 1; i <= n; i++) {
            if (temp[i] < 0) {
                if (fabs(temp[i]) > max) max = fabs(temp[i]);
            }

        }
        if (max == 0) {
            if (upper_bround == -1) upper_bround = 0;
            else lower_bround = 0;
        }
        else {
            if (upper_bround == -1) upper_bround = 1 + pow((max / temp[0]), 1 / k);
            else lower_bround = -(1 + pow((max / temp[0]), 1 / k));
        }
        for (int i = 0; i <= n; i++) {
            if (i % 2 == 1)
                temp[n - i] = -temp[n - i];
        }
    }
    while (lower_bround > 0);
}
double halving_method(double x0, double x1, double coeff_poly[N], int n, double epsi) {
    double a1 = x0;
    double b1 = x1;
    double c, value, temp;
    int count = 0;
    do {
        count ++;
        c = (a1 + b1) / 2;
        cout << "\nlan lap thu " << count << " x = " << c << endl;
        value = function_value(c, coeff_poly, n);
        temp = function_value(a1, coeff_poly, n);
        if (value == 0)
            break;
        else {
            if (temp * value < 0) b1 = c;
            else a1 = c;
        }
    }
    while (fabs(a1 - b1) >= epsi);
    cout << "\n\nVoi khoang cach li nghiem: (" << setprecision(-log10(epsi) + 1) << x0 << " , " << x1 << ")";
    cout << "\nNghiem cua phuong trinh la: x = " << c;
    cout << "\nSo vong lap: " << count << "\nSai so: " << fabs(a1 - b1);
    cout << "\n\n*****************************************************************************";
    return c;
}
void solve_polynomial(double coeff_poly[N], int n, double epsi) {
    double survey_value[N];
    double x0, x1, sign, temp0, temp1, value1, value2, lower_bround, upper_bround;
    double eta = 1e-11;
    int k = 1;
    domain_solution(lower_bround, upper_bround, coeff_poly, n);
    if (n == 0) {
        if (coeff_poly[0] == 0)
            cout <<"\nPhuong trinh vo so nghiem !!";
        else
            cout <<"\nPhuong trinh vo nghiem !!";
    }
    else if (n == 1)
        cout << "\nPhuong trinh co nghiem duy nhat: " << - coeff_poly[1] / coeff_poly[0];
    else if (lower_bround == upper_bround) {
        if (function_value(lower_bround, coeff_poly, n) == 0) {
            cout << "\nNghiem cua phuong trinh la: " << lower_bround;
        }
        else
            cout <<"\nPhuong trinh vo nghiem !!";
    }
    else {
        x1 = lower_bround;
        while (x1 < upper_bround) {
            x0 = x1;
            temp0 = derivative(x0, coeff_poly, n);
            if (temp0 < 0)
                sign = -1;
            else
                sign = 1;
            x1 = x0 + sign * eta * temp0;
            temp1 = derivative(x1, coeff_poly, n);
            for (int i = 0; i < max_loop; i++) {
                if (temp0 * temp1 > 0) {
                    while (eta < 0.008) {
                        eta = eta * 2;
                        x1 = x0 + sign * eta * temp0;
                        if (derivative(x1, coeff_poly, n) * temp0 < 0) {
                            eta = eta / 2;
                            break;
                        }
                    }
                }
                else {
                    while (eta > 0) {
                        eta = eta / 2;
                        x1 = x0 + sign * eta * temp0;
                        if (derivative(x1, coeff_poly, n) * temp0 > 0)
                            break;
                    }
                }
                x1 = x0 + sign * eta * temp0;
                x0 = x1;
                if (fabs(derivative(x1, coeff_poly, n)) < 1e-4) {
                    survey_value[k] = x1;
                    k++;
                    break;
                }
                eta = 1e-11;
                temp0 = derivative(x0, coeff_poly, n);
                x1 = x0 + sign * eta * temp0;
                temp1 = derivative(x1, coeff_poly, n);
                if (x1 > upper_bround)
                    break;
            }
            x1 = x1 + 0.001;
        }
        survey_value[0] = lower_bround;
        survey_value[k] = upper_bround;
        cout << "\nMien chua nghiem cua phuong trinh da thuc la: ";
        cout << "\nCan duoi : " << lower_bround << "\nCan tren : " << upper_bround;
        cout << "\n*****************************************************************************";
        for (int i = 0; i < k; i++) {
            value1 = function_value(survey_value[i], coeff_poly, n);
            value2 = function_value(survey_value[i + 1], coeff_poly, n);
            if (fabs(value1) <= epsi) {
                cout << "\nNghiem cua phuong trinh gan voi cuc tri, hoac chinh la cuc tri: Nghiem do la: ";
                cout << survey_value[i];
            }
            else if (fabs(value2) <= epsi) {
                cout << "\nNghiem cua phuong trinh gan voi cuc tri, hoac chinh la cuc tri: Nghiem do la: ";
                cout << survey_value[i];
                i++;
            }
            else if (value1 * value2 < 0) {
                halving_method(survey_value[i], survey_value[i + 1], coeff_poly, n, epsi);
            }
        }
    }
}
