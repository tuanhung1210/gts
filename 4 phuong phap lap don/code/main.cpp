//#include <bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#define  delta  1.0e-6
#define  eps  1.0e-5
using namespace std;
double a,b,xo,q;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
double g(double x) // nhap ham g(x)
{
	double mu=(double) 1/3;
	double t=4*pow(x,2) - x - 1;
	if (t<0) return -exp((log(-t)*mu));
	else
		return exp((log(t)*mu));
   // return exp(log(4*pow(x,2)-x-1) * coe);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
double abs_g_phay(double x) // khoi tao ham |g'(x)|
{
    double dy = g(x + delta) - g(x - delta), dx = 2 * delta;
    return fabs(dy / dx);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
double max_abs_g_phay(double a, double b) // khoi tao ham tim q = max|g'(x)| tren [a,b]
{
    double max = abs_g_phay(a), i = a + eps;
    do {
        if (abs_g_phay(i) > max)
            max = abs_g_phay(i);
        i += eps;
    }
    while (i <= b);
    return max;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
// int kiemtra(double a,double b,double xo) // check input
// {
// 	if (a>=b) printf("Dau vao khong hop le! Yeu cau can duoi < can tren (a<b)\n");
// 	else if ((xo<a)||(xo>b)) printf ("Xo nam ngoai [a,b]. Yeu cau nhap lai\n");
// 	else return 0;
// 	printf("Nhap lai can duoi a=");scanf("%lf",&a);
// 	printf("Nhap lai can tren b=");scanf("%lf",&b);
// 	printf("Nhap lai xap xi dau xo=");scanf("%lf",&xo);
// 	kiemtra(a,b,xo);
//     return 1;
// }
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
// double lapdontiennghiem(double xo, double q,double epsilon)
// {
// 	double x=g(xo);
// 	if (x==xo) printf("Tinh theo cong thuc sai so tien nghiem, nghiem cua phuong trinh chinh la xo= %.15lf\n",xo);
// 	else 
// 	{
// 		int i,n;
// 		epsilon*=(1-q);
// 		double step=ceil(log(fabs(epsilon/(x-xo)))/log(q));
// 		n=(int) step;
// 		printf("Theo cong thuc sai so tien nghiem-----------------------------------------------------------------\n");
// 		printf("X1= %.15lf\n",x);
// 		for (i=2;i<=n;i++)
// 		{
// 			x=g(x);
// 			printf("X%d= %.15lf\n",i,x);
// 		}
// 		printf("So lan lap : %d\n",n);
// 		printf("Nghiem gan dung X= %.15lf\n",x);
// 	}	
// }
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
void lapdonhaunghiem(double xo,double q,double epsilon) {
    fstream out;
    out.open("output.txt", ios::out);
    out << "\nPhương pháp lặp đơn giải gần đúng nghiệm bài toán f(x) = 0";
    out << "\nInput: ";
    out << "\nPhương trình: x^3 - 4x^2 + x + 1 = 0";
    out << "\nKhoảng cách li nghiệm: (a, b) = (" << a << ", " << b << ")";
    out << "\nSai số epsilon: " << epsilon << "\n";
    out << "\nChạy code theo công thức sai số hậu nghiệm \n";
    out << "\nChọn điểm Fourier xấp xỉ đầu bất kì : " << xo;
    out << "\nGiá trị co q = " << setprecision(20) << q;
    double x = g(xo);
    int step = 1;
    double del = (epsilon * q) / (1 - q);
    out << "\n\nLần lặp thứ " << step << ": nghiệm x = " << fixed << setprecision(15) << x << ",  sai số = " << del;
    //printf("\n\nLan lap thu %d: nghiem x = %.15lf,  sai so = %.15lf ", step, x, del);
    while (fabs(x - xo) >= epsilon) {
        xo = x;
        x = g(xo);
        step += 1;
        out << "\nLần lặp thứ " << step;
        out << ": nghiệm x = " << setprecision(20) << x;
        out << ",  sai số = " << setprecision(20) << fabs(x - xo);
        //printf("\nLan lap thu %d: nghiem x = %.15lf,  sai so = %.15lf ", step, x, del);
    }
    out << "\n\nKết luận: ";
    out << "\nSố lần lặp: " << step << "\nNghiệm gần đúng: " << x;
   // printf("\n\nKet luan: \nSo lan lap : %d", step);
    //printf("\nNghiem gan dung X= %.15lf", x);
    out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
void dieu_kien_co(double q, double epsilon) // kiem tra dieu kien co
{
    if ((q <= 0) || (q >= 1)) {
        printf("Ham g(x) khong thoa man dieu kien anh xa co\n");
    }
    if((g(a) - a) * (g(b) - b) > 0) {
        cout << "khoang li nghiem khong hop le";
    }
    printf("Ham g(x) thoa man dieu kien co\n");
    //lapdontiennghiem(xo, q, epsilon);
    lapdonhaunghiem(xo, q, epsilon);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
int main()
{
    double epsilon;
	printf("Nhap can duoi a=");scanf("%lf",&a);
	printf("Nhap can tren b=");scanf("%lf",&b);
    printf("Nhap vao sai so epsilon =");scanf("%lf", &epsilon);
    xo = (a + b)/2;
	q=max_abs_g_phay(a,b);
	printf("\nGia tri co: %.15lf\n", q);
	dieu_kien_co(q, epsilon);
    return 0;
}


