#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>

using namespace std;

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"


class Exception : public exception {
public:
    Exception(const string& msg) : msg_(msg) {}

    string getMessage() {return(msg_);}

private:
    string msg_;
};

double module(double x) {
    return (x > 0) ? x : -x;
}

class Matrix {
protected:
    int rows;
    int columns;
    double* matrix;

public:
    Matrix(int n = 1, int m = 1) {
        rows = n;
        columns = m;
        matrix = new double[m * n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                *(this->matrix + i * n + j) = 0;
            }
        }
    }

    void copy(Matrix M) {
        rows = M.getRows();
        columns = M.getColumns();
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                *(start() + j * columns + i) = *(M.start() + j * columns + i);
            }
        }
    }


    Matrix& operator = (Matrix M2) {
        rows = M2.getRows();
        columns = M2.getColumns();
        matrix = M2.start();

        return *this;
    }

    Matrix operator - (Matrix M2) {
        if (rows != M2.getRows() || columns != M2.getColumns()) {
            throw Exception("Error: the dimensional problem occured");
        }
        Matrix res(rows, columns);
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                *(res.start() + i * rows + j) = *(matrix + i * rows + j) - *(M2.start() + i * rows + j);
            }
        }
        return res;
    }

    Matrix operator + (Matrix M2) {
        if (rows != M2.getRows() || columns != M2.getColumns()) {
            throw Exception("Error: the dimensional problem occured");
        }
        Matrix res(rows, columns);
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                *(res.start() + i * rows + j) = *(matrix + i * rows + j) + *(M2.start() + i * rows + j);
            }
        }
        return res;
    }

    Matrix operator * (Matrix M2) {
        if (columns != M2.getRows()) {
            throw Exception("Error: the dimensional problem occured");
        }
        int res_c = M2.getColumns();
        Matrix res(rows, res_c);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < res_c; j++) {
                for (int k = 0; k < columns; k++) {
                    *(res.start() + i * res_c + j) += *(matrix + i * columns + k) * *(M2.start() + k * res_c + j);
                }
            }
        }
        return res;
    }

    friend istream& operator >>(istream &in, Matrix m) {
        for (int i = 0; i < m.columns; i++) {
            for (int j = 0; j < m.rows; j++) {
                in >> *(m.start() + i * m.rows + j);
            }
        }
        return in;
    }

    friend ostream& operator <<(ostream &out, Matrix m) {
        out.setf(ios::fixed);
        out.precision(4);
        for (int i = 0; i < m.getRows(); i++) {
            for (int j = 0; j < m.getColumns(); j++) {
                if (module(*(m.start() + i * m.getColumns() + j)) < pow(0.1, 10)) {
                    out << "0.0000" << " ";
                } else out << *(m.start() + i * m.getColumns() + j) << " ";
            }
            out << endl;
        }
        return out;
    }

    Matrix transpose() {
        int new_columns = rows;
        int new_rows = columns;
        Matrix m(new_rows, new_columns);
        for (int i = 0; i < new_rows; i++ ) {
            for (int j = 0; j < new_columns; j++) {
                *(m.start() + i * new_columns + j) = *(matrix + j * new_rows + i);
            }
        }
        return m;
    }
    int getRows() const {
        return rows;
    }
    int getColumns() const {
        return columns;
    }
    double* start() {
        return matrix;
    }
};


vector<double> split(vector<string>& A) {
    vector<double> v;
    for (int i = 0; i < A.size(); i++) {
        char buff[A[i].size() + 1];
        strcpy(buff, A[i].c_str());
        char *tmp_char = strtok(buff, " ");
        while (tmp_char != NULL) {
            string s = tmp_char;
            v.push_back(atof(s.c_str()));
            tmp_char = strtok(NULL, " ");
        }
    }
    return v;
}


class SquareMatrix: public Matrix {
public:
    SquareMatrix(int n = 1): Matrix(n, n) {
    }
};


class IdentityMatrix: public Matrix {
public:
    IdentityMatrix(int n = 1): Matrix(n, n) {
        for (int j = 0; j < rows; j++) {
            *(this->matrix + j * rows + j) = 1;
        }
    }
};


class EliminationMatrix: public IdentityMatrix {
public:
    EliminationMatrix(int r, int col, int pivot_r, Matrix A): IdentityMatrix(A.getColumns()) {
        int columns = A.getColumns();
        double ratio = *(A.start() + columns * r + col) / *(A.start() + columns * pivot_r + col);
        *(this->start() + columns * r + col) = -ratio;
    }
};


class PermutationMatrix: public IdentityMatrix {
public:
    PermutationMatrix(int r1, int r2, int n): IdentityMatrix(n) {
        for (int i = 0; i < n; i++) {
            double temp = *(this->start() + r2 * columns + i);
            *(this->start() + r2 * columns + i) = *(this->start() + r1 * columns + i);
            *(this->start() + r1 * columns + i) = temp;
        }
    }
};


class ColumnVector: public Matrix {
public:
    int size = 0;

    ColumnVector(Matrix M) {
        rows = M.getRows();
        columns = M.getColumns();
        if (M.getColumns() != 1) {
            throw Exception("No viable conversion from Matrix (m x {}) to ColumnVector (m x 1)");
        }
        matrix = M.start();
    }

    ColumnVector(int r): Matrix(r, 1) {
    }

    Matrix norm() {
        double len = 0;
        Matrix norm(rows, 1);
        for (int j = 0; j < rows; j++) {
            len += pow(*(this->start() + j * columns), 2);
        }
        len = sqrt(len);
        for (int j = 0; j < rows; j++) {
            *(norm.start() + j * columns) = *(start() + j * columns) / len;
        }
        return norm;
    }

    friend ostream& operator << (ostream& out, ColumnVector c) {
        out.setf(ios::fixed);
        out.precision(4);
        for (int j = 0; j < c.rows; j++) {
            if (module(*(c.start() + j)) < pow(0.1, 10))
                out << "0.0000" << endl;
            else out << *(c.start() + j) << endl;
        }
        return out;
    }

    void add(double elem) {
        if (size == rows) {
            throw Exception("Cannot add elements to a filled ColumnVector");
        }
        *(start() + size) = elem;
        size++;
    }

};


/// Least Squares Matrix
class LSMatrix: public Matrix {
public:
    LSMatrix(int m, int n, vector<double> d): Matrix(m, n + 1) {
        for (int k = 0; k < m; k++) {
            *(start() + k * (n + 1)) = 1;
        }
        for (int i = 1; i < n + 1; i++) {
            int ind = 0;
            for (int j = 0; j < m; j++) {
                *(start() + j * (n + 1) + i) = pow(d[ind], i);
                ind += 2;
            }
        }
    }
};

vector<int> findi(Matrix M) {
    vector<int> v;
    int cols = M.getColumns();
    for (int i = 0; i < cols; i++) {
        for (int j = i + 1; j < cols; j++) {
            if (*(M.start() + cols * j + i) != 0) {
                v.push_back(j);
                v.push_back(i);
                return v;
            }
        }
    }
}


void printAugmented(Matrix M, Matrix A) {
    cout.setf(ios::fixed);
    cout.precision(2);
    int n = M.getRows();
    int cols = M.getColumns();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < cols; j++) {
            cout << *(M.start() + i * cols + j) << " ";
        }
        for (int j = 0; j < n; j++) {
            cout << *(A.start() + i * cols + j) << " ";
        }
        cout << endl;
    }
}


void permute(Matrix& M, int cur_col, int& step, Matrix& A) {
    int cols = M.getColumns();
    Matrix P = SquareMatrix(cols);
    int change_col = cur_col;
    double max_pivot = module(*(M.start() + cur_col * cols + cur_col));
    for (int i = cur_col + 1; i < cols; i++) {
        if (module(*(M.start() + i * cols + cur_col)) > max_pivot) {
            max_pivot = module(*(M.start() + i * cols + cur_col));
            change_col = i;
        }
    }
    if (change_col != cur_col) {
        P = PermutationMatrix(change_col, cur_col, cols);
        step++;
        M = P * M;
        A = P * A;
        /*printf("step #%d: permutation", step);
        cout << endl;
        cout << M;
        cout << A;*/
    }
}


void eliminate(int r, int col, Matrix& M, int& step, Matrix& A) {
    Matrix E = EliminationMatrix(r, col, col, M);
    step++;
    M = E * M;
    A = E * A;
    /*printf("step #%d: elimination", step);
    cout << endl;
    cout << M;
    cout << A;*/
}


void normalize(Matrix& M, Matrix& A){
    int cols = M.getColumns();
    int rows = M.getRows();
    for (int i = 0; i < rows; i++) {
        double pivot = *(M.start() + i * cols + i);
        for (int j = 0; j < cols; j++) {
            *(A.start() + i * cols + j) /= pivot;
        }
        *(M.start() + i * cols + i) /= *(M.start() + i * cols + i);
    }
}



Matrix gaussianElimination(Matrix M) {
    Matrix A = IdentityMatrix(M.getColumns());
    int cols = M.getColumns();
    int step = 0;
    for (int i = 0; i < cols - 1; i++) {
        permute(M, i, step, A);
        for (int j = i + 1; j < cols; j++) {
            if (*(M.start() + j * cols + i) != 0 && *(M.start() + i * cols + i) != 0)
                eliminate(j, i, M, step, A);
        }
        int flag = 0;
        for (int k = 0; k < cols; k++) {
            if (module(*(M.start() + (cols - 1) * cols + k)) < pow(0.1, 10)) {
                flag++;
            }
        }
        if (flag == cols) {
            if (module(*(A.start() + cols - 1)) < pow(0.1, 10)) {
                cout << "INF" << endl;
            } else cout << "NO" << endl;
            exit(0);
        }
    }

    for (int i = cols - 1; i >= 1; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (*(M.start() + j * cols + i) != 0) {
                eliminate(j, i, M, step, A);
            }
        }
    }



    for (int j = 0; j < cols; j++) {
        if (*(M.start() + j * cols + j) == 0) {
            if (*(A.start() + j) == 0) {
                cout << "INF" << endl;
                exit(0);
            } else {
                cout << "NO" << endl;
                exit(0);
            }
        }
    }
    normalize(M, A);
    return A;

}


ColumnVector leastSquares(Matrix A, ColumnVector b) {
    Matrix reduced = A.transpose() * A;
    cout << "A_T*A:" << endl;
    cout << reduced;
    ColumnVector c(A.transpose() * b);
    Matrix inverse = gaussianElimination(reduced);
    cout << "(A_T*A)^-1:" << endl;
    cout << inverse;
    cout << "A_T*b:" << endl;
    cout << c;
    return inverse * c;
}


int main() {
    FILE *pipe = _popen(GNUPLOT_NAME, "w");

    int m, n;
    string s;
    vector<string> A;
    vector<double> d;
    cin >> m;
    getline(cin, s);

    try {
        if (pipe != NULL) {

            for (int i = 0; i < m; i++) {
                getline(cin, s);

                A.push_back(s);
            }
            cin >> n;

            d = split(A);

            fprintf(pipe, "%s\n", "set xrange [-20:20]");
            fprintf(pipe, "%s\n", "set yrange [-20:20]");


            ColumnVector c(m);
            for (int i = 1; i < d.size(); i += 2) {
                c.add(d[i]);
            }

            Matrix LS = LSMatrix(m, n, d);
            cout << "A:" << endl;
            cout << LS;

            ColumnVector x = leastSquares(LS, c);
            cout << "x~:" << endl;
            cout << x;

            string f;
            for (int i = 0; i < x.getRows(); i++) {
                double koef = *(x.start() + i);
                string s = to_string(koef) + "*x**" + to_string(x.getRows() - i - 1);
                if (i != x.getRows() - 1) {
                    s += " + ";
                }
                f += s;
            }

            fprintf(pipe, "plot %s with lines\n", f.c_str());

            fprintf(pipe, "replot '-' with points palette cb -5 pointtype 7 pointsize 1\n");

            for (int i = 0; i < m; i++) {
                double x = d.at(i);
                double y = d.at(i+1);

                fprintf(pipe, "%f\t%f\n", x, y);
            }

            fprintf(pipe, "%s\n", "e");

            fflush(pipe);

            _pclose(pipe);
        }


    } catch (Exception &e) {
        cout << e.getMessage() << endl;
    }

    return 0;


}