#pragma once

#include"winbgi2.h"
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include"mesh.h"

///////////////////////////////////////////
//class Vector
template<class T>
class Vector
{
private:
    int size;
    std::vector<T> data;
public:
    ///////////////////////////////////////////
    //creating stuff
    Vector(const int& s = 0) : size(s) { data.resize(size);}
    ~Vector() = default;
    void Fill(T value)
    {
        for (int i = 0; i < size; i++)
            data[i] = value;
    }
    void FillHorizontal(T value)
    {
        for (int i = 0; i < size; i += 2)
            data[i] = value;
    }
    void FillVertical(T value)
    {
        for (int i = 1; i < size; i += 2)
            data[i] = value;
    }
    void ApplyBoundaryCondition(const Mesh<T>& mesh)
    {
        for (int i = 0; i < size / 2; i++)
            if (mesh.GetBoundaryCondition(i))
            {
                data[2 * i] = 0;
                data[2 * i + 1] = 0;
            }
    }
    void ApplyGroundCondition(const Mesh<T>& mesh)
    {
        int rear_wheel_id = 0, front_wheel_id = 0;
        for (int i = 0; i < size / 2; i++)
            if (mesh.GetGroundCondition(i))
                if (rear_wheel_id == 0)
                    rear_wheel_id = i;
                else
                    front_wheel_id = i;

        bool rw_onground = false, fw_onground = false;
        if (mesh.GetNodeY(rear_wheel_id) <= 1)
            rw_onground = true;
        if (mesh.GetNodeY(front_wheel_id) <= 1)
            fw_onground = true;

        double force = 0, momentum = 0;
        for (int i = 0; i < size / 2; i++)
        {
            momentum -= data[2 * i] * (mesh.GetNodeY(i) - 40);
            momentum += data[2 * i + 1] * (mesh.GetNodeX(i) - 446);
            force += data[2 * i + 1];
        }

        data[2 * front_wheel_id + 1] = (force * (mesh.GetNodeX(rear_wheel_id) - 446) - momentum) / (mesh.GetNodeX(front_wheel_id) - mesh.GetNodeX(rear_wheel_id));
        data[2 * rear_wheel_id + 1] = -force - data[2 * front_wheel_id + 1];

        if (data[2 * front_wheel_id + 1] < 0 || fw_onground == false)
            data[2 * front_wheel_id + 1] = 0;
        if (data[2 * rear_wheel_id + 1] < 0 || rw_onground == false)
            data[2 * rear_wheel_id + 1] = 0;
    }

    ///////////////////////////////////////////
    //calculating stuff
    T& CalculateSum()
    {
        T result = 0;
        for (auto it : data)
            result += it;
        return result;
    }
    T& CalculateAverageValue()
    {
        T result = CalculateSum();
        result /= size;
        return result;
    }

    ///////////////////////////////////////////
    //getting access / printing
    void Print() const
    {
        for (int i=0; i<size; i++)
            std::cout << data[i] << std::endl;
    }
    int GetSize() const { return size; }
    const std::vector<T>& GetData() const { return data; }

    ///////////////////////////////////////////
    //operators
	T& operator()(const int& i) { return data[i]; }
    const T& operator()(const int& i) const { return data[i]; }
    Vector& operator+=(const Vector& b)
    {
        try
        {
            if (size != b.size)
                throw(0);
            
            for (int i = 0; i < size; i++)
                data[i] += b.data[i];
            return *this;
        }
        catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
    }
    Vector& operator-=(const Vector& b)
    {
        try
        {
            if (size != b.size)
                throw(0);

            for (int i = 0; i < size; i++)
                data[i] -= b.data[i];
            return *this;
        }
        catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
    }
    Vector& operator*=(const T& b)
    {
        for (int i = 0; i < size; i++)
            data[i] *= b;
        return *this;
    }
    Vector& operator=(const Vector& b)
    {
        try
        {
            if (size != b.size)
                throw(0);
            for (int i = 0; i < size; i++)
                data[i] = b.data[i];
            return *this;
        }
        catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
    }
    template<class T> friend Vector<T> operator*(const Vector<T>&, const T&);
    template<class T> friend T& operator*(const Vector<T>&, const Vector<T>&);
    template<class T> friend Vector<T> operator-(const Vector<T>&, const Vector<T>&);
    template<class T> friend Vector<T> operator+(const Vector<T>&, const Vector<T>&);
};

template<class T>
Vector<T> operator*(const Vector<T>& v1, const T& a)
{
    Vector<T> result(v1.size);
    for (int i = 0; i < v1.size; i++)
        result(i) = v1(i) * a;
    return result;
}

template<class T>
T& operator*(const Vector<T>& v1, const Vector<T>& v2)
{
    try
    {
        if (v1.size != v2.size)
            throw(0);
        T result = 0;
        for (int i = 0; i < v1.size; i++)
            result += v1.data[i] * v2.data[i];
        return result;
    }
    catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
}

template<class T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2)
{
    try
    {
        if (v1.size != v2.size)
            throw(0);
        Vector<T> result(v1.size);
        for (int i = 0; i < v1.size; i++)
            result(i) = v1(i) - v2(i);
        return result;
    }
    catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
}

template<class T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2)
{
    try
    {
        if (v1.size != v2.size)
            throw(0);
        Vector<T> result(v1.size);
        for (int i = 0; i < v1.size; i++)
            result(i) = v1(i) + v2(i);
        return result;
    }
    catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
}

///////////////////////////////////////////
//class Matrix
template<class T>
class Matrix
{
protected:
    std::vector<T> data;
    std::vector<int> columns;
    std::vector<int> rows;
    int size;
    bool Assembled;
public:
    ///////////////////////////////////////////
    //creating stuff
    Matrix() = default;
    Matrix(const int& s) : size(s), Assembled(0) 
    { 
        data.resize(s);
        columns.resize(s);
        rows.resize(s + 1); 
        for (int i = 0; i < size + 1; i++)
            rows[i] = i;
        for (int i = 0; i < size; i++)
            columns[i] = i;
    }
    virtual ~Matrix() = default;
    virtual void AssembleMatrix(const Mesh<T>&) = 0;
    const int& FindPlace(const int& i, const int& j, bool& ifnew)
    {
        int place = rows[i];
        while (columns[place] < j && place < rows[i + 1] - 1)
            place++;
        if (columns[place] == j)
            ifnew = false;
        else 
            ifnew = true;
        return place;
    }
    void ApplyBoundaryCondition(const Mesh<T>& mesh)
    {
        for (int i = 0; i < size; i++)
            if (mesh.GetBoundaryCondition(i/2) == true)
                for (int j = rows[i]; j < rows[i + 1]; j++)
                    if (columns[j] == i)
                        data[j] = 1;
                    else
                        data[j] = 0;
    }
    
    ///////////////////////////////////////////
    //drawing and printing stuff
    void DrawSparsityPattern(const int& problem_size)
    {
        clearscale();
        scale(0, 0, problem_size, problem_size);
        for (int i=0; i<size; i++)
            for(int j = rows[i]; j < rows[i+1]; j++)
                point(i, problem_size - columns[j]);
    }
    void PrintData()
    {
        for (auto it : data)
            std::cout << it << std::endl;
    }
    void PrintColumns()
    {
        for (auto it : columns)
            std::cout << it << std::endl;
    }
    void PrintRows()
    {
        for (auto it : rows)
            std::cout << it << std::endl;
    }

    ///////////////////////////////////////////
    //operators
    template<class T> friend Vector<T> operator*(const Matrix<T>& mat, const Vector<T>& vec);
    Matrix& operator*=(const T& a)
    {
        for (auto it : data)
            it *= a;
        return *this;
    }
    Matrix& operator+=(const Matrix<T>& mat)
    {
        try
        {
            if (size != mat.size)
                throw std::string("Matrixes do not have the same dimension \n");
            if (Assembled == false || mat.Assembled == false)
                throw std::string("Matrixes have not yet been assembled \n");

            int place;
            bool ifnew;
            for (int i = 0; i < size; i++)
                for (int j = mat.rows[i]; j < mat.rows[i + 1]; j++)
                {
                    place = FindPlace(i, j, ifnew);
                    if (ifnew == false)
                        data[place] += mat.data[j];
                    else
                    {
                        data.insert(data.begin() + place, mat.data[j]);
                        columns.insert(columns.begin() + place, mat.columns[j]);
                        for (int n = i + 1; n < size + 1; n++)
                            rows[n] ++;
                    }
                }
        }
        catch (std::string e) { std::cout << e; }

        return *this;
    }
    Matrix& operator=(const Matrix& mat) = default;
};

template<class T>
Vector<T> operator*(const Matrix<T>& mat, const Vector<T>& vec)
{
    try
    {
        if (mat.size != vec.GetSize())
            throw(0);
        Vector<T> result(vec.GetSize());
        for (int i = 0; i < mat.size; i++)
            for (int j = mat.rows[i]; j < mat.rows[i + 1]; j++)
                result(i) += mat.data[j] * vec(mat.columns[j]);
        return result;
    }
    catch (...) { std::cout << "vectors do not have the same dimension!\n"; exit(0); }
}

template<class T>
class ExactMassMatrix : public Matrix<T>
{
public:
    ExactMassMatrix() = default;
    ExactMassMatrix(const int& s) : Matrix<T>(s) {}
    virtual ~ExactMassMatrix() = default;
    virtual void AssembleMatrix(const Mesh<T>& mesh)
    {
        try
        {
            if (Matrix<T>::Assembled)
                throw std::string("Matrix already assembled \n");

            int place;
            bool ifnew;
            for (auto it = mesh.elements.begin(); it != mesh.elements.end(); it++)
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                    {
                        if (it->local_m[i][j] != 0)
                        {
                            place = Matrix<T>::FindPlace(it->GetDOF(i), it->GetDOF(j), ifnew);
                            if (ifnew == false)
                                Matrix<T>::data[place] += it->local_m[i][j];
                            else
                            {
                                Matrix<T>::data.insert(Matrix<T>::data.begin() + place, it->local_m[i][j]);
                                Matrix<T>::columns.insert(Matrix<T>::columns.begin() + place, it->GetDOF(j));
                                for (int n = it->GetDOF(i) + 1; n < Matrix<T>::size + 1; n++)
                                    Matrix<T>::rows[n] ++;
                            }
                        }
                    }
            Matrix<T>::Assembled = true;
        }
        catch (std::string e) { std::cout << e; }
    }
};

template <class T>
class ApproximateMassMatrix : public Matrix<T>
{
public:
    ApproximateMassMatrix() = default;
    ApproximateMassMatrix(const int& s) : Matrix<T>(s) {}
    virtual ~ApproximateMassMatrix() = default;
    virtual void AssembleMatrix(const Mesh<T> & mesh)
    {
        try
        {
            if (Matrix<T>::Assembled)
                throw std::string("Matrix already assembled \n");

            for (auto it = mesh.elements.begin(); it != mesh.elements.end(); it++)
                for (int i = 0; i < 6; i++)
                    Matrix<T>::data[it->GetDOF(i)] += it->GetMass() / 3;

            Matrix<T>::Assembled = true;
        }
        catch (std::string e) { std::cout << e; }
    }
};

template<class T>
class StiffnessMatrix : public Matrix<T>
{
public:
    StiffnessMatrix() = default;
    StiffnessMatrix(const int& s) : Matrix<T>(s) {}
    virtual ~StiffnessMatrix() = default;
    virtual void AssembleMatrix(const Mesh<T>& mesh)
    {
        try
        {
            if (Matrix<T>::Assembled)
                throw std::string("Matrix already assembled \n");

            int place;
            bool ifnew;
            for (auto it = mesh.elements.begin(); it != mesh.elements.end(); it++)
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                    {
                        if (it->local_k[i][j] != 0)
                        {
                            place = Matrix<T>::FindPlace(it->GetDOF(i), it->GetDOF(j), ifnew);
                            if (ifnew == false)
                                Matrix<T>::data[place] += it->local_k[i][j];
                            else
                            {
                                Matrix<T>::data.insert(Matrix<T>::data.begin() + place, it->local_k[i][j]);
                                Matrix<T>::columns.insert(Matrix<T>::columns.begin() + place, it->GetDOF(j));
                                for (int n = it->GetDOF(i) + 1; n < Matrix<T>::size + 1; n++)
                                    Matrix<T>::rows[n] ++;
                            }
                        }
                    }
            Matrix<T>::Assembled = true;
        }
        catch (std::string e) { std::cout << e; }
    }
};

///////////////////////////////////////////
//solving equations
template<class T>
void ConjugateGradient(const Matrix<T>& A, Vector<T>& x, const Vector<T>& b, const int& size)
{
    Vector<T> p(size);
    Vector<T> Ap(size);
    Vector<T> res(size);
    T norm_res;
    T alpha, beta = 1, gamma, delta;
    const double eps = 1e-3;
    const int max_iter = 10000;
    int iter = 0;
    
    res = b - (A*x);
    p = res;
    delta = res * res;
    for (;;)
    {
        Ap = A * p;
        alpha = delta / (p * Ap);
        x += p * alpha;
        res -= Ap * alpha;
        gamma = res * res;
        beta = gamma / delta;
        p *= beta;
        p += res;
        norm_res = sqrt(gamma);
        delta = gamma;
        if (norm_res < eps)	//checking converge condition
        {
            std::cout << "calculation converged! \n iteration performed: " << iter << "\n residuum obtained: " << norm_res << std::endl;
            break;
        }
        if (iter == max_iter - 1)	 //checking number of iterations condition
        {
            std::cout << "calculation did not converge! \n iteration performed: " << iter << "\n residuum obtained: " << norm_res << std::endl;
            break;
        }
        iter++;
    }
}

template<class T>
void Solve(Mesh<T>& mesh, const int& iter)
{
    const int problem_size = 2 * mesh.GetNoNodes();
    const T dt = 0.1;
    
    StiffnessMatrix<double> stiff(problem_size);
    stiff.AssembleMatrix(mesh);

    //ExactMassMatrix<double> mass(problem_size);
    ApproximateMassMatrix<double> mass(problem_size);
    mass.AssembleMatrix(mesh);

    Vector<double> x(problem_size);
    Vector<double> v(problem_size);
    v.FillHorizontal(100);
    v.ApplyBoundaryCondition(mesh);
    Vector<double> f(problem_size);
    f.FillVertical(-5);
    f = mass * f;
    f.ApplyBoundaryCondition(mesh);

    StiffnessMatrix<double> M(problem_size);
    M = stiff;
    M *= dt * dt * 0.01;
    M += mass;
    M.ApplyBoundaryCondition(mesh);
    Vector<double> b(problem_size);

    std::ofstream out_displacements, out_stresses;
    try
    {
        out_displacements.open("displacements.txt");
        out_stresses.open("stresses.txt");
        if (out_displacements.is_open() == false)
            throw std::string("displacements.txt");
        if (out_stresses.is_open() == false)
            throw std::string("stresses.txt");
        for (int i = 0; i < iter; i++)
        {
            f.ApplyGroundCondition(mesh);
            b = f - stiff * x;
            b *= dt;
            b += mass * v;
            ConjugateGradient(M, v, b, problem_size);
            x = v * dt;

            for (int j = 0; j < problem_size; j += 2)
                out_displacements << x(j) << "\t" << x(j + 1) << std::endl;

            mesh.CalculateStresses(x.GetData());
            for (int j = 0; j < mesh.GetNoElements(); j++)
                out_stresses << mesh.GetElementStress(j) << std::endl;

        }
        out_displacements.close();
        out_stresses.close();
    }
    catch (std::string file_name)
    {
        std::cout << "file " << file_name << " not opened\n";
        exit(0);
    }
}

template<class T>
void Animate(Mesh<T> mesh, const int& iter)
{
    const int problem_size = 2 * mesh.GetNoNodes();
    T stress, max_stress, min_stress;
    Vector<T> x(problem_size);
    std::ifstream in_displacements, in_stresses;
    //finding min and max stress for scale
    try
    {
        in_stresses.open("stresses.txt");
        if (in_stresses.is_open() == false)
            throw std::string("stresses.txt");
        in_stresses >> stress;
        max_stress = stress;
        min_stress = stress;
        for (int i = 1; i < iter * mesh.GetNoElements(); i++)
        {
            in_stresses >> stress;
            if (stress > max_stress)
                max_stress = stress;
            if (stress < min_stress)
                min_stress = stress;
        }
        in_stresses.close();
    }
    catch (std::string file_name)
    {
        std::cout << "file " << file_name << " not opened\n";
        exit(0);
    }

    //drawing
    Vector<double> start(problem_size);
    start.FillHorizontal(-500);
    mesh.AddDisplacement(start.GetData());
    Vector<double> speed(problem_size);
    speed.FillHorizontal(20);
    for (int i = 0; i < 500; i += 20)
    {
        mesh.AddDisplacement(speed.GetData());
        for (int j = 0; j < mesh.GetNoElements(); j++)
            mesh.AddElementStress(1, j);
        mesh.ColorMesh(1, 9e9);
        animate(2);
    }
    try
    {
        in_displacements.open("displacements.txt");
        in_stresses.open("stresses.txt");
        if (in_displacements.is_open() == false)
            throw std::string("displacements.txt");
        if (in_stresses.is_open() == false)
            throw std::string("stresses.txt");

        for (int i = 0; i < iter; i++)
        {
            for (int j = 0; j < problem_size; j += 2)
                in_displacements >> x(j) >> x(j + 1);
            for (int j = 0; j < mesh.GetNoElements(); j++)
            {
                in_stresses >> stress;
                mesh.AddElementStress(stress, j);
            }
            mesh.AddDisplacement(x.GetData());
            mesh.ColorMesh(min_stress, max_stress);
            animate(2);
        }

        in_displacements.close();
        in_stresses.close();
    }
    catch (std::string file_name)
    {
        std::cout << "file " << file_name << " not opened\n";
        exit(0);
    }
}