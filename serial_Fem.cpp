#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include <chrono>
#include <omp.h>

struct mesh
{
    std::vector<double> p;
    std::vector<std::vector<int>> t;
};

double Norm(std::vector<double> vec)
{
    double results = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        results += vec[i] * vec[i];
    }

    return sqrt(results);
}

double Dot(std::vector<double> v1, std::vector<double> v2)
{
    double results = 0;
    for (int i = 0; i < v1.size(); i++)
    {
        results += v1[i] * v2[i];
    }

    return results;
}

int lin_comb(std::vector<double> &out, std::vector<double> &v1, double b, std::vector<double> &v2)
{

    for (int i = 0; i < v1.size(); i++)
    {
        out[i] = v1[i] + b * v2[i];
    }

    return 0;
}

template <typename T>
void printVec(T out)
{
    int i = 0;
    for (i = 0; i < out.size() - 1; i++)
    {
        printf("%8.4f", (double)out[i]);
    }

    printf("%8.4f\n", (double)out[i]);
}

int Quadrature_1D_nodes(std::vector<double> elem, std::vector<double> &nodes)
{

    nodes[0] = -0.93246951420315202781;
    nodes[1] = -0.66120938646626451366;
    nodes[2] = -0.23861918608319690863;
    nodes[3] = 0.23861918608319690863;
    nodes[4] = 0.66120938646626451366;
    nodes[5] = 0.93246951420315202781;

    for (int i = 0; i < nodes.size(); i++)
    {
        nodes[i] = 0.5 * (elem[1] - elem[0]) * nodes[i] + 0.5 * (elem[1] + elem[0]);
    }

    return 0;
}

int Quadrature_1D_weights(std::vector<double> elem, std::vector<double> &weights)
{

    weights[0] = 0.17132449237917034504;
    weights[1] = 0.36076157304813860757;
    weights[2] = 0.46791393457269104739;
    weights[3] = 0.46791393457269104739;
    weights[4] = 0.36076157304813860757;
    weights[5] = 0.17132449237917034504;

    for (int i = 0; i < weights.size(); i++)
    {
        weights[i] = ((elem[1] - elem[0]) / 2.0) * weights[i];
    }

    return 0;
}

double shape_fun_1D_Lagrange(double x, std::vector<double> elem, int degree, int shape_index, int dx)
{
    if (dx == 0)
    {
        if (shape_index == 1)
            return ((x - elem[1]) / (elem[0] - elem[1]));
        else
            return ((x - elem[0]) / (elem[1] - elem[0]));
    }
    if (dx == 1)
    {
        if (shape_index == 1)
            return (1.0 / (elem[0] - elem[1]));
        else
            return (1.0 / (elem[1] - elem[0]));
    }

    return 0.0;
}

mesh mesh_generator_1D(std::vector<double> interval, int n)
{
    std::vector<double> p(n + 1, 0.0);
    std::vector<std::vector<int>> t(n);

    mesh m;
    m.p = p;
    m.t = t;

    for (int i = 0; i < n + 1; i++)
    {
        std::vector<int> elem = {i, i + 1};
        m.p[i] = interval[0] + (i / static_cast<double>(n)) * (interval[1] - interval[0]);

        if (i < n)
            m.t[i] = elem;
    }

    return m;
}

int FE_matrix_1D_Lagrange_tri_loc(std::vector<std::vector<double>> &sol,
                                  std::function<double(double)> func, std::vector<double> elem, int degree, int dx1, int dx2)
{
    int Quad_size = 6;
    std::vector<double> nodes(Quad_size, 0);
    std::vector<double> weights(Quad_size, 0);
    Quadrature_1D_nodes(elem, nodes);
    Quadrature_1D_weights(elem, weights);

    int idof = (degree + 1);

    std::vector<std::vector<double>> M(idof, std::vector<double>(idof, 0));
    std::vector<double> ibas_val(Quad_size, 0);
    std::vector<double> jbas_val(Quad_size, 0);

    std::vector<double> fun_val(Quad_size, 0);
    for (int k = 0; k < Quad_size; k++)
    {
        fun_val[k] = func(nodes[k]);
    }

    for (int i = 1; i <= idof; i++)
    {

        for (int k = 0; k < Quad_size; k++)
        {
            ibas_val[k] = shape_fun_1D_Lagrange(nodes[k], elem, degree, i, dx2);
        }

        for (int j = 1; j <= idof; j++)
        {

            for (int k = 0; k < Quad_size; k++)
            {
                jbas_val[k] = shape_fun_1D_Lagrange(nodes[k], elem, degree, j, dx1);
            }

            for (int k = 0; k < Quad_size; k++)
            {
                M[i - 1][j - 1] += weights[k] * fun_val[k] * ibas_val[k] * jbas_val[k];
            }
        }
    }

    sol = M;

    return 0;
}

int FE_matrix_1D(std::vector<std::vector<double>> &sol,
                 std::function<double(double)> func, mesh ms, int degree, int dx1, int dx2)
{
    int size = (degree * ms.t.size() + 1);
    int idof = (degree + 1);

    std::vector<std::vector<double>> M(size, std::vector<double>(size, 0));

    std::vector<std::vector<double>> local_stiffness;

    int nextrow = 0;

    for (int i = 0; i < ms.t.size(); i++)
    {

        std::vector<int> index = ms.t[i];
        std::vector<double> elem = {ms.p[index[0]], ms.p[index[1]]};
        FE_matrix_1D_Lagrange_tri_loc(local_stiffness, func, elem, degree, dx1, dx2);

        for (int k = 0; k < idof; k++)
        {
            for (int j = 0; j < idof; j++)
            {

                M[nextrow + k][nextrow + j] += local_stiffness[k][j];
            }
        }
        nextrow += degree;
    }

    sol = M;

    return 0;
}

int FE_vec_1D_Lagrange_tri_loc(std::vector<double> &sol,
                               std::function<double(double)> func, std::vector<double> elem, int degree, int dx)
{
    int Quad_size = 6;
    std::vector<double> nodes(Quad_size, 0.0);
    std::vector<double> weights(Quad_size, 0.0);
    Quadrature_1D_nodes(elem, nodes);
    Quadrature_1D_weights(elem, weights);

    int idof = (degree + 1);

    std::vector<double> b(idof, 0.0);
    std::vector<double> ibas_val(Quad_size, 0.0);

    std::vector<double> fun_val(Quad_size, 0.0);
    for (int k = 0; k < Quad_size; k++)
    {
        fun_val[k] = func(nodes[k]);
    }
    // std::cout << "idof: " << idof << std::endl;
    for (int i = 1; i <= idof; i++)
    {
        for (int k = 0; k < Quad_size; k++)
        {
            ibas_val[k] = shape_fun_1D_Lagrange(nodes[k], elem, degree, i, dx);
        }

        for (int k = 0; k < Quad_size; k++)
        {
            b[i - 1] += weights[k] * fun_val[k] * ibas_val[k];
        }
    }

    sol = b;

    return 0;
}

int FE_vec_1D(std::vector<double> &sol,
              std::function<double(double)> func, mesh ms, int degree, int dx1)
{
    int size = (degree * ms.t.size() + 1);
    int idof = (degree + 1);

    std::vector<double> b(size, 0.0);
    std::vector<double> local_vec;

    std::vector<int> index;
    int nextrow = 0;

    for (int i = 0; i < ms.t.size(); i++)
    {

        index = ms.t[i];
        std::vector<double> elem = {ms.p[index[0]], ms.p[index[1]]};
        FE_vec_1D_Lagrange_tri_loc(local_vec, func, elem, degree, dx1);

        for (int j = 0; j < idof; j++)
        {

            b[nextrow + j] += local_vec[j];
        }
        nextrow += degree;
    }

    sol = b;

    return 0;
}

double FE_evaluate_1D(std::vector<double> &b, double x, mesh ms, int degree, int dx)
{

    int idof = (degree + 1);

    std::vector<double> vu_loc(idof, 0.0);
    std::vector<int> index;
    std::vector<double> elem;

    int i = 0;
    for (i = 0; i < ms.t.size(); i++)
    {

        index = ms.t[i];
        std::vector<double> element = {ms.p[index[0]], ms.p[index[1]]};
        if (x >= element[0] && x <= element[1])
        {
            elem = element;
            break;
        }
    }

    int offset = i * degree;
    for (int j = 0; j < idof; j++)
    {
        vu_loc[j] = b[offset + j];
    }

    double total = 0.0;


    for (i = 1; i <= idof; i++)
    {
        total += vu_loc[i - 1] * shape_fun_1D_Lagrange(x, elem, degree, i, dx);
    }

    return total;
}


std::vector<double> conjugate_gradient(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& b,
                                         int max_iter = 1000,
                                         double tol = 1e-10)
{
    int n = static_cast<int>(b.size());
    std::vector<double> x(n, 0.0);      
    std::vector<double> r = b;          
    std::vector<double> p = b;          
    std::vector<double> Ap(n, 0.0);
    double rsold = 0.0;

   
    
    for (int i = 0; i < n; ++i)
    {
        rsold += r[i] * r[i];
    }

    for (int k = 0; k < max_iter; ++k)
    {
        
        
        for (int i = 0; i < n; ++i)
        {
            Ap[i] = 0.0;
            for (int j = 0; j < n; ++j)
            {
                Ap[i] += A[i][j] * p[j];
            }
        }

     
        double dot_pAp = 0.0;
       
        for (int i = 0; i < n; ++i)
        {
            dot_pAp += p[i] * Ap[i];
        }

    
        double alpha = rsold / dot_pAp;

   
        
        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = 0.0;
        
        for (int i = 0; i < n; ++i)
        {
            rsnew += r[i] * r[i];
        }

      
        if (std::sqrt(rsnew) < tol)
            break;

        
        
        for (int i = 0; i < n; ++i)
        {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
    }

    return x;
}

int getSubmatrix(std::vector<std::vector<double>> &M)
{
    std::vector<std::vector<double>> d(M.size() - 2, std::vector<double>(M.size() - 2, 0));
    for (int i = 0; i < M.size() - 2; i++)
    {
        for (int j = 0; j < M.size() - 2; j++)
        {
            d[i][j] = M[i + 1][j + 1];
        }
    }
    M = d;
    return 0;
}

int getSubvec(std::vector<double> &b)
{
    std::vector<double> bhat(b.begin() + 1, b.begin() + b.size() - 1);
    b = bhat;

    return 0;
}

int main(int argc, char *argv[])
{
    double s = 0.95;  //edit this for the value that you want to evaluate
    auto f = [](double x)
    { return 2.0; };        //edit this for your right-hand side

    auto a = [](double x)
    { return 1.0; };    //Any positive function or constant will do
    

    omp_set_num_threads(8);
    auto start = std::chrono::high_resolution_clock::now();
    int degree = 1;
    std::vector<std::vector<double>> m;
    std::vector<double> in = {0.0, 1.0};
    mesh ms = mesh_generator_1D(in, 1000);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time for Mesh Generation: " << duration.count() << " ms\n";
    

    
    start = std::chrono::high_resolution_clock::now();
    FE_matrix_1D(m, a, ms, degree, 1, 1);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Elapsed time for Matrix Generation: " << duration.count() << " ms\n";
    
    

    getSubmatrix(m);
    std::vector<double> b;
    FE_vec_1D(b, f, ms, degree, 0);
    getSubvec(b);


    start = std::chrono::high_resolution_clock::now();
    std::vector<double> x = conjugate_gradient(m, b);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Elapsed time for Linear Solver: " << duration.count() << " ms\n";

    std::cout << "Solution: " << FE_evaluate_1D(x, s, ms, degree, 0) << std::endl;

    return 0;
}
