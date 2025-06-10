using namespace std;
double eps;
int iterations;
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "matrix.h" 
#define SumNorm 3
#define Vao "input.txt"
#define Ra "output.txt"

vector<int> arr;    // luu dia chi cua gia tri lon nhat trong hang hoac cot

int soptu()     //dem so phan tu trong file input
{
    fstream ab;
    double a;
    int dem=0;
    ab.open(Vao, ios::in);
    while(!ab.eof())
    {
        ab>>a;
        dem++;
    }
    ab.close();
    return dem;
}

template <typename T>
T Pow(const T &A ,int n) // ham binh phuong nhanh
{
    if (n==1) return A;
    T tmp = Pow(A, n>>1);
    return  (n&1) ? tmp * tmp * A : tmp * tmp ;
}

/*template<typename T>
int MinNorm(const matrix<T> &A, double &q)
{
    int k = 1;
    double crr[5];
    for(int i = 1; i <= SumNorm; ++i) crr[i] = getNorm(A,i);
    q = crr[1];
    for(int i = 2; i <= SumNorm; i++)
    {
        if(crr[i] < q);
        {
            q = crr[i];
            k = i;
        }
    }
    if(q < 1 ) return k;
    return 0;
}*/
template<typename T>
void Swap(T &a, T &b)       // doi vi tri cá»§a a va b
{
    T c = a;
    a = b;
    b = c;
}

template<typename T>
bool checkRow(const matrix<T> &A)       // kiem tra xem co dua duoc ve chuan hang khong. Luu dia chi troi
{
    arr.resize(A.row + 1);
    vector<int> use;
    use.resize(A.row + 1);
    for(int i = 0; i < use.size(); i++) 
        use[i] = 0;
    double norm1;
    double Max;
    int k = -1;
    for (int i = 0; i < A.row; ++i)
    {
        Max = 0;
        norm1 = 0;
        for (int j = 0; j < A.col; ++j)
        {
            {
                norm1 = norm1 + abs(A[i][j]);
                if(abs(A[i][j]) > Max)
                {
                    Max = abs(A[i][j]);
                    k = j;
                }
            }
        }
        if(abs(A[i][k]) <= norm1 / 2 || use[k] == 1)
        {
            return false;
        }
        else
        {
            arr[i] = k;
            use[k] = 1;
        }

    }
    return true;
}

template<typename T>
bool checkCol(const matrix<T> &A)   //ham kiem tra xem dua duoc ve chuan cot hay khong. Luu dia chi troi
{
    return checkRow(!A);
}

template<typename T>
int dominant(const matrix<T> &A)        //ham kiem tra tinh troi cua ma tran
{
    double a = 0;
    bool check = true;
    for (int i = 0; i < A.row; ++i)
    {
        a = abs(A[i][i]);
        for (int j = 0; j < A.col; ++j)
        {
            if (i == j)
                continue;
            a = a - abs(A[i][j]);
            if (a <= 0)
            {
                check = false;
                break;
            }
        }
    }
    if (check) return 1;        // chuan vo cung
    check = true;
    for (int i = 0; i < A.col; ++i)
    {
        a = abs(A[i][i]);
        for (int j = 0; j < A.row; ++j)
        {
            if (i == j)
                continue;
            a = a - abs(A[j][i]);
            if (a <= 0)
            {
                check = false;
                break;
            }
        }
    }
    if (check) return 2;        // chuan 1
    if(checkRow(A)) return 3;       // dua ve duoc chuan vo cung
    if(checkCol(A)) return 4;       // dua ve duoc chuan 1
    return 0;
}

template<typename T>
void Check(matrix<T> &A, matrix<T> &B, matrix<T> &X0, double &eps)      //ham kiem tra dau vao
{
    cout<<endl;
    if (A.col != A.row)
    {
        cout << "Ma tran khong vuong!!";
        exit(0);
    }
    if(A.row != B.row)
    {
        cout<<"So hang cua A hoac B khong hop le";
        exit(0);
    }
    if(X0.col != B.col)
    {
        cout << "So cot cua B hoac X khong hop le!!";
        exit(0);
    }
    if (A.row != B.row || A.row != X0.row)
    {
        cout << "So hang cua B hoac X khong hop le!!";
        exit(0);
    }
    if (eps < 0)
    {
        cout << "Sai so khong hop le";
        exit(0);
    }

    int a;
    a = A.row * A.col + B.row * B.col + 7;
    int b = soptu();
    if(a != b)
    {
        cout<<"Input khong hop le";
        exit(0);
    }

}
template<typename T>
matrix<T> loop(matrix<T> &A, matrix<T> &B, matrix<T> &X0, double &eps, int &type, double &q, int &normType, double w = 1)   //ham lap tim nghiem
{
    fstream wi;
    matrix<T> X(X0.row, X0.col);
    switch (type)
    {
        case 1:
        {
            wi.open(Ra, ios::out);
            matrix<T> X1(X0.row, X0.col);
            double w1 = q / (1 - q);
            X = X0;
            int loopNumber = 0;
            do
            {
                ++loopNumber;
                X1 = X;
                X = (A * X1) + B;
                wi << "Ket qua cua vong lap thu " << loopNumber  << endl;
                wi << fixed << setprecision(-log10(eps) + 3) << X;
            }
            while (w * w1 * getNorm((X-X1),normType) > eps);
            cout << "So vong lap la : " << loopNumber << "\n";
            cout << "he so co q = " << q << endl;
            wi << "Sai so : " << (w * w1 * getNorm((X-X1),normType));
            wi << endl;
            wi.close();
            return X;
            break;
        }
        case 2:
        {
            wi.open(Ra, ios::out);
            matrix<T> X1(X0.row, X0.col);
            X = X0;
            X1 = (A * X) + B;
            int loopNumber = ceil((log((eps * (1 - q)) / (getNorm((X-X1),normType) * w))) / (log(q)));
            for (int i = 1; i <= loopNumber; ++i)
            {
                X = (A * X) + B;
                wi << "Ket qua vong lap thu " << i << endl;
                wi << fixed << setprecision(-log10(eps) + 3) << X;
            }
            cout << "So vong lap la : " << loopNumber << "\n";
            cout << q << endl;
            wi << "Sai so la : " << (w*(pow(q, loopNumber)/(1-q))*getNorm((X1 - X0), normType));
            wi << endl;
            wi.close();
            return X;
            break;
        }
        case 3:
        {
            wi.open(Ra, ios::out);
            matrix<T> X1(X0.row, X0.col);
            X = X0;
            X1 = (A * X) + B;
            int loopNumber = ceil((log((eps * (1 - q)) / (getNorm((X-X1),normType) * w))) / (log(q)));
            matrix<matrix<T> > P(2,2);
            matrix<T> O(A.row, A.col);
            matrix<T> I(A.row, A.col);
            for (int i = 0; i < A.row; i++) I[i][i] = 1;
            P[0][0] = A; P[0][1] = I;
            P[1][0] = O;     P[1][1] = I;
            P = Pow(P,loopNumber);
            X = P[0][0] * X + P[0][1] * B;
            cout << "So vong lap la : " << loopNumber << "\n";
            cout << q << endl;
            wi << fixed << setprecision(-log10(eps) + 3) << X;
            wi << "Sai so la : " << (w*(pow(q, loopNumber)/(1-q))*getNorm((X1 - X0), normType));
            wi << endl;
            wi.close();
            return X;
            break;
        }
    }
    wi.close();
    return X;
}
template<typename T>
matrix<T> singleloop(matrix<T> &A, matrix<T> B, matrix<T> &X0, double &eps) //lap don
{
    Check(A, B, X0, eps);
    int type;
    double q;
    int normType = 0;
    matrix<T> C(A.row, A.col);
    for(int i = 0; i < A.row; ++i)
        for(int j = 0; j < A.col; ++j)
        {
            if(i == j)
                C[i][j] = 1 - A[i][j];
            else
                C[i][j] = -A[i][j];
        }
    //normType = MinNorm(C, q);
    for(int i = 1; i <= SumNorm; i++)
    {
        if(getNorm(C, i) < 1)
        {
            normType = i;
            q = getNorm(C, i);
            break;
        }
    }
    if (normType != 0)
    {
        cout << endl << "1.Hau nghiem" << endl;
        cout << "2.Tien nghiem" << endl;
        cout << "3.Tien nghiem 2" << endl;
        cout << "Chon cong thuc sai so : ";
        cin >> type;
        return loop(C, B, X0, eps, type, q, normType);
    }
    for(int i = 0; i < A.row; ++i)
    {
        for(int j = 0; j < A.col; ++j)
        {
            if(i == j)
                C[i][j] = 1 + A[i][j];
            else
                C[i][j] = A[i][j];
        }
    }
    //normType = MinNorm(C, q);
    for(int i = 1; i <= SumNorm; i++)
    {
        if(getNorm(C, i) < 1)
        {
            normType = i;
            q = getNorm(C, i);
            break;
        }
    }
    if(normType == 0)
    {
        cout << "Khong the giai bang lap don!!";
        exit(0);
    }
    else
    {
        for(int i = 0; i < B.row; ++i)
            for(int j = 0; j < B.col; ++j)
            {
                B[i][j] = -B[i][j];
            }
    }
    cout << endl << "1.Hau nghiem" << endl;
    cout << "2.Tien nghiem" << endl;
    cout << "3.Tien ngiem 2" << endl;
    cout << "Chon cong thuc sai so : ";
    cin >> type;
    return loop(C, B, X0, eps, type, q, normType);
}
template<typename T>
matrix<T> jacobiloop1(matrix<T> A, matrix<T> B, matrix<T> &X0, double &eps) // lap jacobi
{
    Check(A, B, X0, eps);
    double q;
    int type;
    int normType = dominant(A);
    if (normType == 0)
    {
        cout << "Ma tran khong troi";
        exit(0);
    }
    cout << endl << "1.Hau nghiem" << endl;
    cout << "2.Tien nghiem" << endl;
    cout << "3.Tien nghiem 2" << endl;
    cout << "Chon cong thuc sai so : ";
    cin >> type;
    matrix<T> C(A.row, A.col), D(B.row, B.col);
    if(normType == 3 || normType ==4)
    {
        for(int i = 0; i < A.row; ++i)
        {
            int a = 0;
            if(a != arr[a])
            {
                A[a].swap(A[arr[a]]);
                B[a].swap(B[arr[a]]);
                Swap(arr[a], arr[arr[a]]);
            }
            else a++;
        }
        normType = normType - 2;
    }
    int Domi = normType;
    for (int i = 0; i < A.row; ++i)
    {
        for (int j = 0; j < A.col; ++j)
        {
            if(i!=j)
            {
                C[i][j] = -(A[i][j]) / (A[i][i]);
            }
        }
    }
    for (int i = 0; i < B.row; ++i)
    {
        for (int j = 0; j < B.col; ++j)
        {
            D[i][j] = (B[i][j]) / (A[i][i]);
        }
    }
    double w = 1;
    q = getNorm(C, normType);
    //normType = MinNorm(C, q);
    if(Domi == 2)
    {
        matrix<T> F(A.row, A.col);
        for (int i = 0; i < A.row; ++i)
            for (int j = 0; j < A.col; ++j)
            {
                if(i!=j){
                    F[i][j] = -(A[i][j]) / (A[j][j]);
                }
            }
        q = getNorm(F,normType);
        //normType = MinNorm(F, q);
        double t1 = abs(A[0][0]);
        double t2 = abs(A[0][0]);
        for (int i = 1; i < A.col; i++)
        {
            t1 = max(t1, abs(A[i][i]));
            t2 = min(t2, abs(A[i][i]));
        }
        w = t1 / t2;
    }
    return loop(C, D, X0, eps, type, q, normType, w);
}

template<typename T>
matrix<T> loopIterations(matrix<T> &A, matrix<T> &B, matrix<T> &X0, int iterations, int normType = 1)
{
    fstream wi;
    wi.open(Ra, ios::out); // Ra: tên file output (d?m b?o b?n dã khai báo tru?c)

    matrix<T> X(X0.row, X0.col);
    matrix<T> X1(X0.row, X0.col);

    X = X0;

    for (int i = 1; i <= iterations; ++i)
    {
        X1 = X;
        X = (A * X1) + B;
        double error = getNorm(X - X1, normType); // Tính sai s? gi?a 2 l?n l?p

        wi << "Ket qua cua vong lap thu " << i << endl;
        wi << fixed << setprecision(6) << X;
        wi << "Sai so: " << error << "\n\n";
    }

    wi.close();
    return X;
}



template<typename T>
matrix<T> jacobiloop2(matrix<T> A, matrix<T> B, matrix<T> &X0, int iterations) // l?p Jacobi theo s? bu?c
{
    Check(A, B, X0, eps); // eps không c?n thi?t n?a
    int type;
    int normType = dominant(A);
    if (normType == 0)
    {
        cout << "Ma tran khong troi";
        exit(0);
    }

    matrix<T> C(A.row, A.col), D(B.row, B.col);

    if (normType == 3 || normType == 4)
    {
        for (int i = 0; i < A.row; ++i)
        {
            int a = 0;
            if (a != arr[a])
            {
                A[a].swap(A[arr[a]]);
                B[a].swap(B[arr[a]]);
                Swap(arr[a], arr[arr[a]]);
            }
            else a++;
        }
        normType = normType - 2;
    }

    int Domi = normType;

    for (int i = 0; i < A.row; ++i)
        for (int j = 0; j < A.col; ++j)
            if (i != j)
                C[i][j] = -(A[i][j]) / (A[i][i]);

    for (int i = 0; i < B.row; ++i)
        for (int j = 0; j < B.col; ++j)
            D[i][j] = B[i][j] / A[i][i];

    double q = getNorm(C, normType);
    double w = 1;

    if (Domi == 2)
    {
        matrix<T> F(A.row, A.col);
        for (int i = 0; i < A.row; ++i)
            for (int j = 0; j < A.col; ++j)
                if (i != j)
                    F[i][j] = -(A[i][j]) / (A[j][j]);

        q = getNorm(F, normType);

        double t1 = abs(A[0][0]), t2 = abs(A[0][0]);
        for (int i = 1; i < A.col; i++)
        {
            t1 = max(t1, abs(A[i][i]));
            t2 = min(t2, abs(A[i][i]));
        }
        w = t1 / t2;
    }

    return loopIterations(C, D, X0, iterations); // G?i hàm l?p v?i s? l?n c? d?nh
}

int main(){
    int select;
    int n, m, k, t, a, b;
    cout << "1.singleloop" << endl;
    cout << "2.jacobiloop sai so" << endl;
    cout << "3.jacobiloop so lan lap" << endl;
    cout << "Nhap lua chon : ";
    cin >> select;
    fstream re,wi;
    re.open(Vao,ios::in);
    wi.open(Ra,ios::app);
    if (select == 1)
    {
        cout << "Giai pt AX+B " << endl;
        cout << "Nhap vao so hang va so cot cua ma tran A : " << endl;
        re >> n >>m;
        cout <<"Nhap vao so hang va so cot cua X0: " << endl;
        re >> k >> t;
        cout <<"Nhap vao so hang va so cot cua B: " << endl;
        re >> a >>b;
        matrix<double> A(n, m), B(a, b), X0(k, t), X(k, t);
        cout << "Nhap ma tran A : " << endl;
        re >> A;
        cout << "Nhap ma tran B : " << endl;
        re >> B;
        cout << "Nhap vao sai so : ";
        re >> eps;
        X = singleloop(A, B, X0, eps);
        wi << "Ket qua cua X :" << endl;
        wi << fixed << setprecision(-log10(eps) + 3) << X;
        wi << "Ket qua cua A * X -B la : " << endl;
        wi << setprecision(-log10(eps) + 3) << (A * X - B);
    }
    else if (select == 2)
    {
        cout << "Giai pt AX=B " << endl;
        cout << "Nhap vao so hang va so cot cua ma tran A : " << endl;
        re >> n >>m;
        cout <<"Nhap vao so hang va so cot cua X0: " << endl;
        re >> k >> t;
        cout <<"Nhap vao so hang va so cot cua B: " << endl;
        re >> a >>b;
        matrix<double> A(n, m), B(a, b), X0(k, t), X(k, t);
        cout << "Nhap ma tran A : " << endl;
        re >> A;
        cout << "Nhap ma tran B : " << endl;
        re >> B;
        cout << "Nhap vao sai so : ";
        re >> eps;
        X = jacobiloop1(A, B, X0, eps);
        wi << "Ket qua cua X :" << endl;
        wi << fixed << setprecision(-log10(eps) + 3) << X;
        wi << "ket qua cua A * X - B la : " << endl;
        wi << setprecision(-log10(eps) + 3) << (A * X - B);
    }
    else {
    	cout << "Giai pt AX=B " << endl;
        cout << "Nhap vao so hang va so cot cua ma tran A : " << endl;
        re >> n >>m;
        cout <<"Nhap vao so hang va so cot cua X0: " << endl;
        re >> k >> t;
        cout <<"Nhap vao so hang va so cot cua B: " << endl;
        re >> a >>b;
        matrix<double> A(n, m), B(a, b), X0(k, t), X(k, t);
        cout << "Nhap ma tran A : " << endl;
        re >> A;
        cout << "Nhap ma tran B : " << endl;
        re >> B;
        cout << "Nhap vao so lan : ";
        re >> iterations;
        eps = 0;
        X = jacobiloop2(A, B, X0, iterations);
        wi << "Ket qua cua X :" << endl;
        wi << fixed << setprecision(7) << X;			//<=== chinh so sau dau phay 
        wi << "ket qua cua A * X - B la : " << endl;	// Neu de bao so chu so co nghia thi chinh so sau dau phay to ra
        wi << setprecision(7) << (A * X - B);			// roi tu viet so chu so co nghia vao bai thi
		}	
    wi.close();
    re.close();
    system("notepad output.txt");
    return 0;
}
