#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define LD long double

using namespace std;

const LD h = 1e-7;
LD eta = 0.001;
const int MAX_LOOP = 1000000;
LD a, b, x, x_old, fa, fb, fx, eps, m1, M1;

LD f(LD x) 											//f(x)=0
{
    return pow(x,3) - 0.2*pow(x,2) + 1;
}

LD f1(LD x) 										//f`(x)
{
    return (f(x + h) - f(x - h)) / (2 * h);
}

LD f2(LD x) 										//f``(x)
{
    return (f1(x + h) - f1(x - h)) / (2 * h);
}

LD random_num(LD a, LD b)
{
    return a + (LD)rand() / RAND_MAX * (b - a);
}

int re_sign(LD (*g)(LD x)) 							//Gradient-Decent
{
    int i;
    for (i = 1; i < 100; i++)						//Kiểm tra hàm đang xét có là hàm hằng không
        if (g(random_num(a, b)))
            break;
    if (i == 100)
        return 0;
    LD gx, sign;
    x = a;
    gx = g(x);
    if (gx == 0)
    {
        x += eps;
        gx = g(x);
    }
    sign = (gx > 0) ? 1 : -1;
    int cnt = 0;
    while (x < b && ++cnt <= MAX_LOOP) 				//Gradient-Decent rất hiếm khi có trường hợp g(x)=0
    {								   				//Tuy nhiên để đề phòng lặp vô hạn tôi vẫn giới hạn số lấn lặp
        x += eta * fabs(gx);
        if (x > b)
            x = b;
        gx = g(x);
        if (gx==0) {x+=eps;continue;}				//Vạn nhất g(x) = 0 thật thì tôi vẫn muốn hàm cho KQ đúng
        if (sign * gx < 0)
            return 1;
    }
    return 0;
}

LD sai_so_1()
{
    return fabs(fx) / m1;
}

LD sai_so_2()
{
    return (M1-m1) / m1 * fabs(x-x_old);
}

LD sai_so;

void The_Secant_Method() 							//PP_dây_cung
{
    fstream out;
    out.open("output.txt", ios::out);
    x = a - fa * (a - b) / (fa - fb);
    fx = f(x);
    if (f(x) * fa < 0)								//Xác định điểm Fourier
    {												//Nếu a là điểm Fourier gán b = a
        b = a;										//Như vậy điểm b luôn là điểm Fourier
        fb = fa;
    }
    out << "\nPhương pháp dây cung giải gần đúng nghiệm bài toán f(x) = 0";
    out << "\nInput: ";
    out << "\nPhương trình: x^3 - 4x^2 + x + 1 = 0";
    out << "\nKhoảng cách li nghiệm: (a, b) = (" << a << ", " << b << ")";
    out << "\nSai số epsilon: " << eps << "\n";
    out << "\nChạy code theo công thức sai số hai lần lặp liên tiếp: \n";
    out << "\nĐiểm Fourier : " << b <<"\n";
    int cnt = 0;
    sai_so = sai_so_2();							//muốn dùng sai số nào chỉnh ở đây
    //cout << cnt <<". f(" << x << ") = " << fx <<  " | sai so : "<<sai_so<<'\n';
    while (sai_so > eps && cnt <= MAX_LOOP)
    {
        cnt ++;
        x_old = x;
        x = x - fx * (x - b) / (fx - fb);
        fx = f(x);
        sai_so = sai_so_2();//và ở đây
        out << "\nLần lặp thứ " << cnt << ": nghiệm x = " << fixed << setprecision(15) << x << ", Sai số: " << sai_so;
        //cout << cnt <<". f(" << x << ") = " << fx <<  " | sai so : "<<sai_so<<'\n';
    }
    out << "\n\nKết luận: \nSố lần lặp: " << cnt << "\nNghiệm gần đúng: " << x;
}

void Input()
{
    // freopen("INP.txt", "r", stdin);freopen("OUT.txt", "w", stdout);ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0); //đặt dòng này vào comment nếu muốn nhập từ console
    cout << "\nnhap a: ";
    cin >> a;
    cout << "\nnhap b: ";
    cin >> b;
    cout << "\nnhap eps: ";
    cin >> eps;
    if (a > b)
        swap(a, b);
    //cout << fixed << setprecision(-log10(eps) + 3);
}
int Input_Invalid() //kiểm tra dữ liệu đầu vào
{
    fa = f(a);
    fb = f(b);
    if (fa * fb > 0)
    {
        cout << "khoang li nghiem ko hop le";
        return 1;
    }
    if (re_sign(f1))
    {
        cout << "ham khong don dieu";
        return 1;
    }
    if (re_sign(f2))
    {
        cout << "dao ham khong don dieu";
        return 1;
    }

    m1 = fabs(f1(a));
    M1 = fabs(f1(b));
    if (m1 > M1)
        swap(m1, M1);
    return 0;
}
int main()
{
    Input();

    srand(time(NULL));							//câu lệnh này chỉ dùng cho hàm random không liên quan thuật toán chính

    if (Input_Invalid())
        return 0;

    The_Secant_Method();

    return 0;
}
