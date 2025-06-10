#include<bits/stdc++.h>
#define N 105

using namespace std;

long int n, k;
double a[N][N], b[N], c[N][N], x[N], y[N];
double epsi, lamda, s, zeta;

bool row_diagonal () { // kiem tra cheo troi hang
    double temp;
    for(int i = 0; i < n; i++) {
        temp = 0;
        for(int j = 0; j < n; j++)
            if  (j != i) temp = temp + fabs(a[i][j]);
        if (abs(a[i][i]) <= temp) return false;
    }
    return true;
}

bool col_diagonal () { // kiem tra cheo troi cot
    double temp;
    for(int i = 0; i < n; i++) {
        temp = 0;
        for(int j = 0; j < n; j++)
            if (j != i) temp = temp + fabs(a[j][i]);
        if(abs(a[i][i]) <= temp) return false;
    }
    return true;
}


double lamda_cot_chuan_cot () { // tinh he so lamda cua dieu kien troi cot
    double lamdaa = 0;
    double gamma, beta, temp;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (i == j) c[i][j] = 0;
            else c[i][j] = -a[i][j]/a[i][i];
        }
    }
    for(int j = 0; j < n; j++) {
        gamma = 0;
        beta = 0;
        for (int i = 1; i <= j - 1; i++)    // tinh beta
            beta = beta + fabs(c[j][i]);
        for (int i = j; i <= n; i++)    // tinh gamma
            gamma = gamma + fabs(c[j][i]);
        temp = gamma / (1 - beta);
        if (temp > lamdaa) lamdaa = temp;
    }
    return lamdaa;
}

void C_row_matrix () {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (i == j) c[i][j] = 0;
            else c[i][j] = -a[i][j]/a[i][i];
        }
    }
}

double lamda_cal () {   // tinh he so lamda cua dieu kien dung troi hang
    double lamdaa = 0;
    double gamma, beta, temp;
    C_row_matrix ();
    for(int i = 0; i < n; i++) {
        gamma = 0;
        beta = 0;
        for (int j = 1; j <= i - 1; j++)    // tinh beta
            beta = beta + fabs(c[i][j]);
        for (int j = i; j <= n; j++)    // tinh gamma
            gamma = gamma + fabs(c[i][j]);
        temp = gamma / (1 - beta);
        if (temp > lamdaa) lamdaa = temp;
    }
    return lamdaa;
}

double norm_row() {     // tinh chuan hang
    double norm = 0;
    double temp;
    C_row_matrix ();
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 0; j < n; j++) {
            temp += fabs(c[i][j]);
        }
        if(norm < temp)
            norm = temp;
    }
    return norm;
}

void C_col_matrix () {      // dua ve ma tran tuong duong troi cot
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(j != i) c[i][j] = -a[i][j]/a[j][j];
            else c[i][j] = 0;
        }
    }
}

double s_cal() {    // tinh he so s cua dieu kien dung troi cot
    double s, temp;
    C_col_matrix();
    for(int j = 0; j < n; j++) {
        temp = 0;
        for(int i = j+1; i < n; i++)
            temp = temp + fabs(c[i][j]);
        if(temp > s) s = temp;
    }
    return s;
}

double zeta_cal () {    // tinh he so zeta cua dieu kien dung troi cot
    double zeta;
    double temp1 = 0, temp2 = 0, temp3 = 0;
    C_col_matrix();
    for(int j = 0; j < n; j++) {
        temp1 = 0, temp2 = 0;
        for(int i = 0; i <= j; i++) {
            temp1 = temp1 + fabs(c[i][j]);    
        }
        for(int i = j + 1; i < n; i++){
            temp2 = temp2 + fabs(c[i][j]);
        }
        temp3 = temp1/(1 - temp2);
        // cout << endl;
        if(temp3 > zeta) zeta = temp3;
    }
    return zeta;
}

double norm_col() {     // tinh chuan cot
    double norm = 0;
    double temp;
    C_col_matrix();
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 0; j < n; j++) {
            temp += fabs(c[j][i]);
        }
        if(norm < temp)
            norm = temp;
    }
    return norm;
}

void row_solve () {     // giai bai toan troi hang bang lap hau nghiem
    double z[N];
    double temp, norm;
    lamda = lamda_cal();    // tinh he so lamda cua dieu kien dung troi hang
    k = 0;      // khoi tao de tinh so lan lap
    cout << "ma tran da cho cheo troi hang\n";
    cout << "ta dua ve ma tran dang sau x = Cx + D: \n";
    C_row_matrix ();
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            cout << fixed << setw(10) << setprecision(4) << c[i][j];
        }
        cout << "|" << fixed << setw(10) << setprecision(4) << b[i]/a[i][i] << endl;
    }
    while (true) {  // lap tim nghiem
        k = k + 1; 
        for(int i = 0; i < n; i++) z[i] = x[i];     // khoi tao gia tri ban dau cua z
        for(int i = 0; i < n; i++) {
            temp = 0;
            for(int j = 0; j < n; j++)
                if (j != i) temp = temp + (a[i][j]/a[i][i]) * x[j];
            x[i] = b[i]/a[i][i] - temp;
        }   
        norm = 0;
        for(int i = 0; i < n; i++) {       // chuan hang cua hieu (z - x)
            if (fabs(z[i] - x[i]) > norm) norm = fabs(z[i] - x[i]);
        }
        cout << "\nLan lap thu " << k << ": \n";
        for(int i = 0; i < n; i++) {
            cout << "x_" << i+1 << ": " <<  setprecision(-log10(epsi) + 3) << x[i] << endl;
        }
        if(lamda*norm < epsi*(1- lamda)) break;     // kiem tra dieu kien dung
    }
    cout << fixed;
    cout << "Nghiem cua phuong trinh la: \n";
    for(int i = 0; i < n; i++)
        cout << "x_" << i+1 << ": " <<  setprecision(-log10(epsi) + 3) << x[i] << endl;
    cout << "So lan lap: " << k;
    cout << "\nChuan cua ma tran duoc su dung la chuan hang co he so co: " << lamda;
    cout << "\nGia tri cua chuan hang la: " << norm_row();
    cout << "\nSai so thuc te la: " << lamda * norm / (1- lamda);
}

void col_solve() {      // giai bai toan troi cot bang lap hau nghiem
    double z[N];
    double norm, temp;
    k = 0;
    s = s_cal();
    zeta = zeta_cal();
    for(int i = 0; i < n; i++)
        y[i] = x[i] * a[i][i];      // dat y_i = a_ii * x_ii
    while(true) {       // lap tim nghiem
        k = k + 1;      // dem so lan lap
        for (int i = 0; i < n; i++) z[i] = y[i];
        for (int i = 0; i < n; i++) {
            temp = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) temp += (a[i][j] / a[j][j]) * y[j];
            }
            y[i] = b[i] - temp;
            x[i] = y[i]/a[i][i];    // dua tu y ve lai x
        }
        cout << "\nLan lap thu " << k << ": \n";
        for(int i = 0; i < n; i++) {
            cout << "x_" << i+1 << ": " <<  setprecision(-log10(epsi) + 3) << x[i] << endl;
        }
        norm = 0;
        for (int i = 0; i < n; i++) {       // tinh chuan cua cot cua (y - x)
            norm += fabs(y[i] - z[i]);
        }
        if (zeta * norm < epsi * (1 - s) * (1 - zeta)) break;
    }
    for(int i = 0; i < n; i++) x[i] = y[i]/a[i][i];
    cout << fixed;
    cout << "nghiem cua phuong trinh la: \n";
    for(int i = 0; i < n; i++)
        cout << "x_" << i+1 << ": " <<  setprecision(-log10(epsi) + 3) << x[i] << endl;
    cout << "so lan lap: " << k;
    cout << "\nChuan cua ma tran duoc su dung la chuan cot: " << zeta;
    cout << "\nGia tri cua chuan cot la: " << norm_col();
    cout << "\nSai so thuc te la: " << zeta * norm / ((1 - s) * (1 - zeta));
}

int main () {
    fstream input;
    char s;
    cout << "Hay nhap input theo thu tu sau vao file input\n";
    cout << "nhap n\n";
    cout << "nhap ma tran A\n";
    cout << "nhap ma tran B\n";
    cout << "nhap epsilon\n";
    cout << "go [Y] neu ban muon tiep tuc: ";
    cin >> s;
    if (s == 'Y' || s == 'y') {
        input.open("input.txt", ios:: in);
        // nhap n
        input >> n;
        // nhap ma tran A
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                input >> a[i][j];
        // nhap vecto B
        for(int i = 0; i < n; i++)
            input >> b[i];
        for(int i = 0; i < n; i++)
            x[i] = 0;
        // nhap epsilon
        input >> epsi;
        if(row_diagonal() == 1)     // kiem tra troi hang
            row_solve();    // lap tim nghiem troi hang
        else if(col_diagonal() == 1)    // kiem tra troi cot
            col_solve();        // lap tim nghiem troi cot
        else cout << "ma tran khong cheo troi";
    }
    return 0;
}
