#define Id(I) for (int i = 0; i < I.row; i++) I[i][i] = 1;      // khai bao ma tran don vi

template <typename T>   // dinh nghia nap chong du lieu
class matrix {      // dinh nghia ma tran
public:
    typedef vector<T> Row;
    vector<Row> data;
    int row, col;
    
    matrix()                     : row(0), col(0), data(0)               {}     //khoi tao ma tran mac dinh
    matrix(T d)                  : row(1), col(1), data(1, Row(1,d))     {}     //khoi tao ma tran mot phan tu la d
    matrix(int r, int c)         : row(r), col(c), data(r, Row(c))       {}     //khoi tao ma tran r hang c cot
    matrix(int r, int c, T d)    : row(r), col(c), data(r, Row(c,d))     {}     //khoi tao ma tran toan d
    matrix(const matrix<T> &A)   : row(A.row), col(A.col) , data(A.data) {}     //khoi tao ma tran bang ma tran
    
    matrix<T> operator=(const matrix<T> &A)     // ma tran a bang ma tran b
    {
        row=A.row; 
        col=A.col;
        data=A.data;
        return A;
    }

    const Row &operator[](int i) const { return data[i]; }
    Row &operator[](int i)             { return data[i]; }

};

// tinh toan cho vecto
template <typename T>
vector<T> operator-(const vector<T> &A)     // gan so am cho mot vecto
{
    vector<T> C(A.size());
    for (int i = 0; i < A.size(); i++) C[i] = -A[i];
    return C;
}

template <typename T>
vector<T> operator+(const vector<T> &A, const vector<T> &B)     // cong hai vecto
{
    int m = min(A.size(), B.size());
    vector<T> C(max(A.size(), B.size()));
    for (int i = 0; i < m; i++) C[i] = A[i] + B[i];
    for (int i = m; i < A.size(); i++) C[i] = A[i];
    for (int i = m; i < B.size(); i++) C[i] = B[i]; 
    return C;
}

template <typename T>
vector<T> operator-(const vector<T> &A, const vector<T> &B)     // tru hai vec to
{
    vector<T> C(max(A.size(), B.size()));
    C = A + (-B);
    return C;
}

template <typename T, typename X>
vector<T> operator*(const vector<T> &A, X B)        // nhan hai vecto
{
    vector<T> C(A.size());
    for (int i = 0; i < A.size(); i++) C[i] = A[i] * B;
    return C;
}


// tinh toan cho ma tran
template <typename T>
matrix<T> operator+(const matrix<T> &A, const matrix<T> &B) // cong hai ma tran
{
    matrix<T> C;
    C.data = A.data + B.data;
    C.row = max(A.row, B.row);
    C.col = max(A.col, B.col);
    return C;
}

template <typename T>
matrix<T> operator-(const matrix<T> &A, const matrix<T> &B)     // tru hai ma tran
{
    matrix<T> C;
    C.data = A.data - B.data;
    C.row = max(A.row, B.row);
    C.col = max(A.col, B.col);
    return C;
}

template <typename T>
matrix<T> operator*(const matrix<T> &A, const matrix<T> &B)     // nhan hai ma tran
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
matrix<T> operator*(const matrix<T> &A, X B)        // nhan ma tran A cho B
{
    matrix<T> C(A.row,A.col);
    C.data = A.data * B;
    return C;
}

template <typename T>
matrix<T> operator|(matrix<T> A, const matrix<T> &B)        // bo sung ma tran B vao ma tran A
{
    A.col+=B.col;
    for (int i = 0; i < A.row; i++)
        A[i].insert(A[i].end(), B[i].begin(), B[i].end());
    return A;
}

template <typename T>
matrix<T> operator!(const matrix<T> &A)     // chuyen vi ma tran A
{
    matrix<T> C(A.col, A.row);
    for (int i = 0; i < C.row; i++)
        for (int j = 0; j < C.col; j++)
        C[i][j] = A[j][i];
    return C;
}

template <typename T>
istream& operator>>(istream& is,matrix<T>& A)
{
    for (int i = 0; i < A.row; i++)
        for (int j = 0; j < A.col; j++)
            is>>A[i][j];
    return is;
}

template <typename T>
ostream& operator<<(ostream& os,const matrix<T>& A)     //In ra định dạng .csv
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

