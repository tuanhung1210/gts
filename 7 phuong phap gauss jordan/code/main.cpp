#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#define ra "output.txt"
#define N 100

using namespace std;

struct matrix{
    double a[N][N];
    int row;
    int col;
};

typedef struct matrix MT;

void input_matrix(string file_name, MT &A);

/// index[j] = i biểu thị vị trí phần tử giải đang ở hàng i cột j
/// index[N] chứa vị trí hàng của phần tử giải

bool index_element_solution(MT A, int index[N], int &row_solution, int &col_solution);

void gauss_Jordan_method(MT A, MT B);
int min(int a, int b);
int main(){
    MT A, B;
    fstream inp;
    double m, n;
    inp.open("input.txt", ios::in);
    if (!inp.is_open()) {
        cout << "\nERROR !! Khong ton tai file !!\n";
        exit(0);
    }
    inp >> m;
    A.row = m;
    inp >> n;
    A.col = n;
    for (int i = 0; i < A.row; i++) { 
        for (int j = 0; j < A.col; j++) {
            inp >> A.a[i][j];
        }
    }
    for (int i = 0; i < A.row; i++) { 
        inp >> B.a[i][0]; 
    }
    inp.close();
    gauss_Jordan_method(A, B);
    return 0;
}

int min(int a, int b){
    if (a >= b) return b;
    else return a;
}

/// index[j] = i biểu thị vị trí phần tử giải đang ở hàng i cột j
/// index[N] chứa vị trí hàng của phần tử giải

bool index_element_solution(MT A, int index[N], int &row_solution, int &col_solution) {
    double max = 0;
    int check;
    for (int i = 0; i < A.row; i++) {
        check = 0;
        for (int m = 0; m < A.col; m++) {
            if (i == index[m]) {
                check = 1;
                break;
            }
        }
        if (check == 1) continue;
        else {
            for (int j = 0; j < A.col; j++) {
                if (index[j] != -1)
                    continue;
                if (fabs(A.a[i][j]) == 1) {
                    row_solution = i;
                    col_solution = j;
                    index[col_solution] = row_solution;
                    return true;
                }
                else {
                    if (fabs(A.a[i][j]) > max) {
                        max = fabs(A.a[i][j]);
                        row_solution = i;
                        col_solution = j;
                    }
                }
            }
        }
    }
    if (max != 0) {
        index[col_solution] = row_solution;
        return true;
    }
    else return false;
}

void gauss_Jordan_method(MT A, MT B) {
    fstream out;
    int index[N], index_row[N], row_solution, col_solution;
    int min_row_col = min(A.row, A.col);
    int rankA = 0, rankAB = 0;
    int pos_col_nghiem_tham_so[N], count = 0;
    double coeff, sum, X[N], temp;
    MT A1 = A;
    MT B1 = B;
    out.open(ra, ios::out);
    out << "\nSử dụng phương pháp Gauss - Jordan giải bài toán Ax = B\n";
    for (int i = 0; i < A.col; i++) {
        index[i] = -1;
    }
    for (int i = 0; i < A.row; i++){
        index_row[i] = -1;
    }

    out << "\nMa trận đầu vào A cỡ " << A.row << "x" << A.col <<": ";
    for (int i = 0; i < A.row; i++){
        out << "\n";
        for (int j = 0; j < A.col; j++){
            out << setw(8) << A.a[i][j] << " ";
        }
    }
    out << "\n\nMa trận đầu vào B cỡ " << A.row << "x" << 0 << ": ";
    for (int i = 0; i < A.row; i++){
        out << "\n" << setw(8) << B.a[i][0];
    }
    out << "\n";
    out << "\nBắt đầu quá trình biến đổi ma trận bổ sung A|B \n";
    for (int m = 0; m < min_row_col; m++) {
        out << "\nLần lặp: " << m << "\n";
        /// In ra ma trận bổ sung
        for (int i = 0; i < A.row; i++){
            for (int j = 0; j < A.col; j++){
                if (j < A.col - 1)
                    out << fixed << setw(10) << setprecision(6) << A.a[i][j] << " ";
                else if (j == A.col - 1)
                    out << fixed << setw(10) << setprecision(6) << A.a[i][j] << " | ";
            }
            out << fixed << setw(10) << setprecision(6) << B.a[i][0] << "\n";
        }
        if (index_element_solution(A, index, row_solution, col_solution)) {
            index_row[row_solution] = col_solution;
            for (int i = 0; i < A.row; i++) {
                if (i != row_solution) {
                    coeff = A.a[i][col_solution] / A.a[row_solution][col_solution];
                    for (int j = 0; j < A.col; j++) {
                        A.a[i][j] = A.a[i][j] - coeff * A.a[row_solution][j];
                    }
                    B.a[i][0] = B.a[i][0] - coeff * B.a[row_solution][0];
                }
            }
        }
        else break;
    }
    ///In ra ma trận bổ sung
    out << "Kết quả sau khi biến đổi\n";
    for (int i = 0; i < A.row; i++){
        for (int j = 0; j < A.col; j++){
            if (j < A.col - 1)
                out << fixed << setw(10) << setprecision(6) << A.a[i][j] << " ";
            else if (j == A.col - 1)
                out << fixed << setw(10) << setprecision(6) << A.a[i][j] << " | ";
        }
        out << fixed << setw(10) << setprecision(6) << B.a[i][0] << "\n";
    }
    /// Xét đến rank của ma trận A (rankA) và ma trận bổ sung A|B (rankAB)
    /// Mảng index_row[N] chứa vị trí cột giải
    for (int i = 0; i < A.row; i++){
        if (index_row[i] != -1){
            rankA ++;
            rankAB ++;
        }
        else {
            if (B.a[i][0] != 0)
                rankAB ++;
        }
    }
    out << "\nHạng của ma trận A là: " << rankA;
    out << "\nHạng của ma trận mở rộng A|B là: " << rankAB;
    if (rankA < rankAB)
        out << "\n\nKẾT LUẬN. VẬY HỆ PHƯƠNG TRÌNH VÔ NGHIỆM !";
    /// Trường hợp có nghiệm (có duy nhất nghiệm hoặc vô số nghiệm)
    else if (rankA == rankAB){
        /// Số cột <= số hàng
        if (rankA == A.col){ /// Hạng của ma trận bằng số ẩn
            out << "\n\nKẾT LUẬN. VẬY HỆ PHƯƠNG TRÌNH CÓ NGHIỆM DUY NHẤT: \n";
            for (int i = 0; i < A.col; i++){
                X[i] =  B.a[index[i]][0] / A.a[index[i]][i];
                out << "X[" << i+1 <<"] = " << fixed << setprecision(10) << B.a[index[i]][0] / A.a[index[i]][i] << "\n";
            }
        }
        else if (rankA <= min_row_col){ /// bao gồm cả trường hợp rankA = số hàng (số hàng <= số cột)
            out << "\n\nTa thấy rank A = rank A|B < số ẩn";
            out << "\n\n=> KẾT LUẬN. VẬY HỆ PHƯƠNG TRÌNH CÓ VÔ SỐ NGHIỆM.\nPhụ thuộc vào số tham số là: " << A.col - rankA;
            out << "\nCác tham số đó là: ";
            /// Cột nào không được chọn là cột giải thì x[i + 1] sẽ là tham số
            for (int i = 0; i < A.col; i++){
                if (index[i] == -1){
                    out << "x[" << i + 1 <<"]; ";
                    pos_col_nghiem_tham_so[count] = i;
                    count ++;
                }
            }
            out << "\nNghiệm của hệ phương trình là: \n";
            for (int i = 0; i < A.col; i++){
                if (index[i] == -1)
                    out << "x[" << i + 1 <<"] = x[" << i + 1 <<"]\n";
                else {
                    out << "x[" << i + 1 <<"] = " << fixed << setprecision(10) << B.a[index[i]][0] / A.a[index[i]][i];
                    for (int j = 0; j < count; j++) {
                        out << " + " << fixed << setprecision(10) << (-A.a[index[i]][pos_col_nghiem_tham_so[j]]) / A.a[index[i]][i] << " * x[" << pos_col_nghiem_tham_so[j] + 1 << "] ";
                    }
                    out <<"\n";
                }
            }
            out << "\nMột nghiệm của hệ phương trình là: \n";
            for (int i = 0; i < A.col; i++) {
                sum = 0;
                if (index[i] == -1) {
                    X[i] = 1;
                    out << "x[" << i + 1 << "] = 1\n";
                }
                else {
                    for (int j = 0; j < count; j++) {
                        sum = sum + (-A.a[index[i]][pos_col_nghiem_tham_so[j]]);
                    }
                    X[i] =  (B.a[index[i]][0] + sum) / A.a[index[i]][i];
                    out << "x[" << i + 1 << "] = " << fixed << setprecision(10) << (B.a[index[i]][0] + sum) / A.a[index[i]][i] << "\n";
                }
            }
        }
        out << "\nKiểm tra lại nghiệm của hệ phương trình: \nlấy AX - B = ";
        for (int i = 0; i < A.row; i++){
            temp = 0;
            for (int j = 0; j < A.col; j++){
                temp += A1.a[i][j] * X[j];
            }
            out << "\n" << temp - B1.a[i][0];
        }
    }
    /// Chú ý: trong code ghi x[i + 1] để in ra màn hình đúng thứ tự x[1] -> n
    /// Vì trong ma trận bắt đầu từ hệ số 0 -> n-1 nên khi muốn in ra đúng thì chỗ in ra là: i + 1
    out.close();
}