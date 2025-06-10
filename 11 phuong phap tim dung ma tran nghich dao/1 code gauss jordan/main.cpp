/// CODE TÌM MA TRẬN NGHỊCH ĐẢO BẰNG PHƯƠNG PHÁP GAUSS - JORDAN

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

#define N 100
#define input "input.txt"

struct matrix{
    double a[N][N];
    int row;
    int col;
};

typedef struct matrix MT;

void input_matrix(string file_name, MT &A);     // nhap input

void output_screen(MT A);
void print_MT_boSung(MT A, MT B);
MT swap_row(MT &A, int x, int y);

MT create_unit_matrix(MT &A);

bool index_element_solution(MT A, int index[N], int &row_solution, int &col_solution);
void inverse_matrix_GJ(MT A);


int main(){
    MT A;
    input_matrix(input, A);
    inverse_matrix_GJ(A);
    return 0;
}

void input_matrix(string file_name, MT &A) {        // nhap input
    fstream inp;
    double c;
    inp.open(file_name, ios::in);
    if (!inp.is_open()) {
        cout << "\nERROR !! Khong ton tai file !!\n";
        exit(0);
    }
    // nhap so hang
    inp >> c;
    A.row = c;
    // nhap so cot
    inp >> c;
    A.col = c;
    // nhap ma tran A
    for (int i = 0; i < A.row; i++)
        for (int j = 0; j < A.col; j++) {
            inp >> c;
            A.a[i][j] = c;
        }
    inp.close();
}

void output_screen(MT A) {      // in ra ma trận
    cout << "******************************************************************************\n";
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.col; j++) {
            cout << fixed << setw(10) << setprecision(4) << A.a[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "******************************************************************************\n";
}

MT swap_row(MT &A, int x, int y) {      // doi hang x voi hang y
    double temp;
    for (int j = 0; j < A.col; j++) {
        temp = A.a[x][j];
        A.a[x][j] = A.a[y][j];
        A.a[y][j] = temp;
    }
    return A;
}

MT create_unit_matrix(MT &A) {      // tạo ma tran don vi
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.col; j++) {
            if (i == j) A.a[i][j] = 1;
            else A.a[i][j] = 0;
        }
    }
    return A;
}

void print_MT_boSung(MT A, MT B){       // in ra ma tran bo sung
    int pos;
    cout << "*******************************************************************************\n";
    for (int i = 0; i < A.row; i++){
        pos = 0;
        for (int j = 0; j < A.row * 2; j++){
            if (j < A.row)
                cout << fixed << setw(10) << setprecision(4) << A.a[i][j];
            if (j == A.row)
                cout << "\t";
            if (j >= A.row) {
                cout << fixed << setw(10) << setprecision(4) << B.a[i][pos];
                pos++;
            }
        }
        cout << "\n";
    }
    cout  << "******************************************************************************\n";
}

bool index_element_solution(MT A, int index[N], int &row_solution, int &col_solution) {     // kiem tra xem co tim duoc phan tu giai nua khong va luu dia chi cua phan tu giai neu tim duoc
    double max = 0;
    int check;
    for (int i = 0; i < A.row; i++) {       
        check = 0;
        for (int m = 0; m < A.col; m++) {       // kiem tra hang i da la phan tu giai hay chua
            if (i == index[m]) {
                check = 1;
                break;
            }
        }
        if (check == 1) continue;
        else {
            for (int j = 0; j < A.col; j++) {
                if (index[j] != -1)     // neu index[j] da co thi tang index
                    continue;
                if (fabs(A.a[i][j]) == 1) {     // luu dia chi phan tu giai
                    row_solution = i;
                    col_solution = j;
                    index[col_solution] = row_solution;
                    return true;        // tra ve la van tim duoc phan tu giai
                }
                else {
                    if (fabs(A.a[i][j]) > max) {        // luu dia chi cua phan tu giai
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
        return true;        // tra ve la van tim duoc phan tu giai
    }
    else return false;      // tra ve la khong tim duoc phan tu giai
}

void inverse_matrix_GJ(MT A) {      // tim nghich dao cua ma tran bang gauss - jordan
    int index[N], row_solution, col_solution, m;
    double coeff;
    MT unit_matrix;
    unit_matrix.row = A.row;
    unit_matrix.col = A.col;
    create_unit_matrix(unit_matrix);        // tao ma tran don vi
    for (int i = 0; i < A.col; i++) {       // khoi tao mang luu dia chi cua phan tu giai 
        index[i] = -1;
    }
    cout << "\nma tran can tim nghich dao la: \n";
    output_screen(A);
    cout << "\nbat dau tim ma tran nghich dao\n";
    for (m = 0; m < A.row; m++) {       
        
        if (index_element_solution(A, index, row_solution, col_solution)) {     // neu tim duoc phan tu giai
            cout << "\nChon phan tu giai la: " << A.a[row_solution][col_solution] << "\n";
            for (int i = 0; i < A.row; i++) {       
                if (i != row_solution) {
                    coeff = A.a[i][col_solution] / A.a[row_solution][col_solution];
                    for (int j = 0; j < A.col; j++) {
                        A.a[i][j] = A.a[i][j] - coeff * A.a[row_solution][j];
                        unit_matrix.a[i][j] = unit_matrix.a[i][j] - coeff * unit_matrix.a[row_solution][j];
                    }
                }
            }
        print_MT_boSung(A, unit_matrix);
        }
        else break;        // neu khong tim thay tuc la da hoan thanh chuan hoa
    }
    
    cout << "\nket thuc bien doi ma tran bang phuong phap gauss - jordan\n";
    cout << endl;
    if (m < A.row)
        cout << "\nDinh thuc cua ma tran = 0 => khong ton tai ma tran nghich dao !!";
    else {
        for (int i = 0; i < A.row; i++) {       // doi vi tri cac hang sao cho a_ii != 0
            if (fabs(A.a[i][i]) < 1e-4) {       // neu a[i][i] == 0, dùng < 1e-4 la vi a[i][i] co sai so
                for (int h = i + 1; h < A.row; h++) {
                    if (fabs(A.a[h][i]) > 1e-4) {       // neu a[i][i] != 0, dùng > 1e-4 la vi a[i][i] co sai so
                        swap_row(A, h, i);  
                        swap_row(unit_matrix, h, i);
                        break;
                    }
                }
            }
        }
        cout << "\nma tran sau khi doi hang va cot: \n";
        print_MT_boSung(A, unit_matrix);
        cout << endl;
        for (int i = 0; i < unit_matrix.row; i++) {     // tinh nghiem cua he phuong trinh sau bien doi gauss - jordan
            for (int j = 0; j < unit_matrix.col; j++) {
                unit_matrix.a[i][j] = unit_matrix.a[i][j] / A.a[i][i];
            }
        }
        cout << "\nMa tran nghich dao la: \n";
        output_screen(unit_matrix);
    }
}