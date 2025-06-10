#define alpha 1e-5

template <typename T>
class matrix
{
public:
    typedef vector<T> Row;
    vector<Row> data;
    int row, col;

    matrix()                     : row(0), col(0), data(0)               {} //khoi tao ma tran
    matrix(T d)                  : row(1), col(1), data(1, Row(1,d))     {} // ma tran 1 phan tu la d
    matrix(int r, int c)         : row(r), col(c), data(r, Row(c))       {} //ma tran toan 0
    matrix(int r, int c, T d)    : row(r), col(c), data(r, Row(c,d))     {} // ma tran toan d
    matrix(const matrix<T> &A)   : row(A.row), col(A.col) , data(A.data) {} //khoi gan cho ma tran

    matrix<T> operator=(const matrix<T> &A)  //nap chong toan tu gan
    {
        row=A.row;
        col=A.col;
        data=A.data;
        return A;
    }
    const Row &operator[](int i) const { return data[i]; } // nap chong toan tu truy cap
    Row &operator[](int i)             { return data[i]; }
};


//operator for vector
template <typename T>
vector<T> operator-(const vector<T> &A)
{
    vector<T> C(A.size());
    for (int i = 0; i < A.size(); i++) 
        C[i] = -A[i];
    return C;
}

template <typename T>
vector<T> operator+(const vector<T> &A, const vector<T> &B)
{
    int m = min(A.size(), B.size());
    vector<T> C(max(A.size(), B.size()));
    for (int i = 0; i < m; i++) C[i] = A[i] + B[i];
    for (int i = m; i < A.size(); i++) C[i] = A[i];
    for (int i = m; i < B.size(); i++) C[i] = B[i];
    return C;
}

template <typename T>
vector<T> operator-(const vector<T> &A, const vector<T> &B)
{
    vector<T> C(max(A.size(), B.size()));
    C = A + (-B);
    return C;
}

template <typename T, typename X>
vector<T> operator*(const vector<T> &A, X B)
{
    vector<T> C(A.size());
    for (int i = 0; i < A.size(); i++) C[i] = A[i] * B;
    return C;
}

//operator for matrix

template <typename T>
matrix<T> operator+(const matrix<T> &A, const matrix<T> &B)
{
    matrix<T> C;
    C.data = A.data + B.data;
    C.row = max(A.row, B.row);
    C.col = max(A.col, B.col);
    return C;
}

template <typename T>
matrix<T> operator-(const matrix<T> &A, const matrix<T> &B)
{
    matrix<T> C;
    C.data = A.data - B.data;
    C.row = max(A.row, B.row);
    C.col = max(A.col, B.col);
    return C;
}

template <typename T>
matrix<T> operator*(const matrix<T> &A, const matrix<T> &B)
{
    matrix<T> C(A.row, B.col);
    C.row = A.row;
    C.col = B.col;
    for (int i = 0; i < C.row; i++)
        for (int j = 0; j < C.col; j++)
            for (int k = 0; k < A.col; k++)
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
    return C;
}

template <typename T,typename X>
matrix<T> operator*(const matrix<T> &A, X B)
{
    matrix<T> C(A.row,A.col);
    C.data = A.data * B;
    return C;
}


template <typename T>
matrix<T> operator!(const matrix<T> &A)
{
    matrix<T> C(A.col, A.row);
    for (int i = 0; i < C.row; i++)
        for (int j = 0; j < C.col; j++)
            C[i][j] = A[j][i];
    return C;
}

template <typename T>
istream& operator>>(istream& is,matrix<T>& A) // nap chong nhap ma tran
{   
    for (int i = 0; i < A.row; i++)
        for (int j = 0; j < A.col; j++)
            is>>A[i][j];
    return is;
}

template <typename T>
ostream& operator<<(ostream& os,const matrix<T>& A) //nap chong in ma tran
{
    for (int i = 0; i < A.row; i++)
    {
        os<<A[i][0];
        for (int j = 1; j < A.col; j++)
            os<<" , "<<A[i][j];
        os<<'\n';
    }
    return os;
}
template<typename T>
double getNorm(const matrix<T> &A, const int &normType) //ham tinh chuan cua ma tran loai chuan 1 la hang, 2 la chuan cot, 3 la chuan euclid, 4 la chuan tri rieng
{
    double norm1;
    double Max = 0;
    switch (normType)
    {
        case 1:
            for (int i = 0; i < A.row; ++i)
            {
                norm1 = 0;
                for (int j = 0; j < A.col; ++j)
                {
                    norm1 = norm1 + abs(A[i][j]);
                }
                Max = max(Max, norm1);
            }
            return Max;
        case 2:
            for (int i = 0; i < A.col; ++i)
            {
                norm1 = 0;
                for (int j = 0; j < A.row; ++j)
                {
                    norm1 = norm1 + abs(A[j][i]);
                }
                Max = max(Max, norm1);
            }
            return Max;
        case 3:
            for (int i = 0; i < A.row; ++i)
            {
                norm1 = 0;
                for (int j = 0; j < A.col; ++j)
                {
                    norm1 = norm1 + A[i][j] * A[i][j];
                }
            }
            Max = sqrt(norm1);
            return Max;
        case 4:
        {
            matrix<T> B = !A;
            B = B * A;
            matrix<T> X(A.row, 1, 1);
            double t0 = 0;
            double t1 = 0;
            do
            {
                t0 = t1;
                X = B * X;
                double s = X[0][0];
                int j;
                for(int i = 1; i < A.row; ++i)
                {
                    if(abs(X[i][0]) > s)
                    {
                        j = i;
                        s = abs(X[i][0]);
                    }
                }
                t1 = X[j][0];
                for(int i = 0; i < A.row; ++i)
                {
                    X[i][0] = X[i][0] / t1;
                }
            } while (abs(t1 - t0) > alpha);
            return sqrt(t1);
        }
        default:
            return -1;
            break;
    }
}

