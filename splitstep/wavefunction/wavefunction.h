#pragma once
#include <complex>

typedef std::complex<double> complex_d;

class wavefunction{
  public:
    //constructor:
    wavefunction();
    wavefunction(const wavefunction& wf);
    wavefunction(int size_n_, int size_m_, double* potential_);

    //destructor:
    ~wavefunction();

    //saving options:
    void saving_options(int frame_rate_, std::string saving_path_, std::string saving_name_);
    void set_frame_rate(int frame_rate_);
    void set_saving_path(std::string saving_path_);
    void set_saving_name(std::string saving_name_);
    void reset_frames();
    std::string generate_path();

    //functions:
    void reset();
    void reset(int i, int j);
    void reset(int i, int j, complex_d* vec_0 ,complex_d* vec_1, int n);
    void reset(int i, int j, complex_d* vec_0, complex_d* vec_1, int n, int m);
   
    void lie_trotter_timestep_schroedinger(double dt, int steps, double direction, bool create_csv = false, int fftwthreads = 1);
    void lie_trotter_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv = false, int fftwthreads = 1);
    void strang_timestep_schroedinger(double dt, int steps, double direction, bool create_csv = false, int fftwthreads = 1);
    void strang_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv = false, int fftwthreads = 1);
    void discrete_timestep_schroedinger(double dt_in, int steps, double direction, bool create_csv);
    void discrete_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv);

    void timestep(std::string potential, std::string method, double dt_in, int steps, double direction, double uu, bool create_csv, int fftwthreads = 1);

    void reduce_noise(double tol);

    //normen:
    double l2_norm() const;
    double l2_norm(const wavefunction& wf) const;
    double l2_max() const;
    void renormalize();

    //return functions:
    int return_n() const;
    int return_m() const;
    complex_d return_psi_0(int i) const;
    complex_d return_psi_1(int i) const;
    void to_csv(const char* name);
    void print() const;
    void print_norm() const;

    //operators:
    wavefunction& operator=(const wavefunction& wf);
    complex_d operator()(int psi, int i, int j) const;
    complex_d& operator()(int psi, int i,int j);

  private:
    int size_n;
    int size_m;
    double* potential;
    complex_d* psi_0;
    complex_d* psi_1;

    //options for printing the wavefunction
    int frame_rate;
    int frame_number; //used to number the frames
    std::string saving_path;
    std::string saving_name;
};
