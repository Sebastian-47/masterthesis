#pragma once

class matrix{
  public:
    //constructor:
    matrix();

    matrix(const matrix& m);

    matrix(int n, int m);

    matrix(int n, int m, double* data_);

    matrix(int size_n_, int size_m_, const char* name);

    //destructor:
    ~matrix();

    //return options;

    int n() const;

    int m() const;

    int length() const;

    void print() const;

    double* return_data() const;

    void to_file(const char* name) const;

    void to_csv(const char* name) const;

    void to_png(const char* name) const;

    double operator()(int a, int b) const;

    double& operator()(int a, int b);

    //operators:
    matrix& operator=(const matrix& m);

    bool operator==(const matrix& m);

    matrix operator%(matrix m);

    matrix operator*(const matrix& m);

    matrix operator*(double d);

    matrix operator/(double d);

    matrix operator+(const matrix& m);

    matrix operator-(const matrix& m);

  private:
    int size_n; 	//number of rows
    int size_m; 	//number of columns
    double* data;	//contains the matrix itself
};

matrix image_to_matrix(const char* name);
