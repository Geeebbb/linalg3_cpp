#include<iostream>
#include<string.h>
#include<cmath>
#include <iomanip>
#include <sstream>
#include<algorithm>
#include <stdexcept>
#include <type_traits>
#include <initializer_list>
#include <stdexcept>
#ifndef MATRIX_H
#define MATRIX_H

using namespace std;
namespace linalg

{
    template<typename T>
    class Matrix
    {
    private:
        T* m_ptr;
        int m_rows = 0;
        int m_columns = 0;
        int m_capacity; //Вместимость выделенной памяти для матрицы
    public:
        T* begin() {
            return m_ptr;
        }

        T* end() {
            return m_ptr + m_rows * m_columns;
        }

        const T* begin() const {
            return m_ptr;
        }

        const T* end() const {
            return m_ptr + m_rows * m_columns;
        }

        template <typename T2>
        friend class Matrix;
        int rows() const noexcept {
            return m_rows;
        }

        int columns() const noexcept
        {
            return m_columns;
        }

        int capacity() const noexcept
        {
            return m_capacity;
        }

        bool empty() const
        {
            if ((m_rows == 0) || (m_columns == 0)) {
                return true;
            }
            else {
                return false;
            }
        }
        void clear();
        void reserve(int n);
        void shrink_to_fit();
        void reshape(int rows, int columns);
        // ƒефолтный конструктор
        Matrix() noexcept : m_ptr(nullptr), m_rows(0), m_columns(0), m_capacity(0) {}
        //  онструкторы с параметрами:
        Matrix(int rows);
        Matrix(int rows, int columns);
        template<typename T2>
        Matrix(initializer_list<T2>); // Конструктор с списком инициализации для столбцов
        template<typename T2>
        Matrix(initializer_list<initializer_list<T2>> lst); // Конструктор с списком инициализации для матрицы*/

        template<typename T2>
        Matrix(const Matrix<T2>& M);
        Matrix(const Matrix<T>& M);// Конструктор копирования
        Matrix(Matrix<T>&& moved) noexcept;

        Matrix<T>& operator = (const Matrix<T>& other);

        template <typename T2>
        Matrix<T>& operator = (const Matrix<T2>& M);

        template <typename T2>
        Matrix<T>& operator = (Matrix<T2>&& M) noexcept;

        T& operator ()(int i, int j); // Оператор доступа к элементу матрицы (i, j)
        T operator ()(int i, int j) const; // Константный оператор доступа к элементу матрицы (i, j)

        //оператор вывода
        template <typename T2>
        friend ostream& operator << (ostream& os, const Matrix<T2>& M);

        //поэлементное сложение матриц
        template <typename T2>
        Matrix<common_type_t<T, T2>> operator + (const Matrix<T2>& M) const; // Оператор сложения
        template <typename T2>
        Matrix<T>& operator += (const Matrix<T2>& M);

        //поэлементное вычитание матриц
        template <typename T2>
        Matrix<common_type_t<T, T2>> operator - (const Matrix<T2>& M) const;
        template <typename T2>
        Matrix<T>& operator -= (const Matrix<T2>& M);

        template <typename T2>
        Matrix<decltype(T{} *T2{}) > operator*(const Matrix<T2>& other) const;

        template <typename T2>
        friend Matrix<decltype(double{} *T2{}) > operator*(double x, const Matrix<T2>& mat);

        Matrix<decltype(T{} *double{}) > operator*(double x) const;

        template <typename T2>
        Matrix<T>& operator*=(const Matrix<T2>& other);

        Matrix<T>& operator*=(double x) noexcept;




        // operators == and !=
        template <typename T2>
        friend bool operator ==(const Matrix<T>& M1, const Matrix<T2>& M2) noexcept
        {
            if (M1.rows() != M2.rows() || M1.columns() != M2.columns()) {
                return false;
            }
            for (int i = 0; i < M1.rows(); ++i) {
                for (int j = 0; j < M1.columns(); ++j) {
                    if ((M1(i, j) - M2(i, j)) != 0) {
                        return false;
                    }
                }
            }
            return true;

        }

        template <typename T2>
        friend bool operator!=(const Matrix<T>& M1, const Matrix<T2>& M2) noexcept
        {
            return !(M1 == M2);
        }

        //norm
        double norm();
        //trace
        T trace() const;
        //det
        T det()const;
        //прямой ход метода Гаусса
        Matrix<T>& gauss_forward();

        Matrix<T>& gauss_backward();

        int rank() const;
        template <typename T2, typename T3>
        friend Matrix<common_type_t<T3, T2>> concatenate(const Matrix<T2>& M1, const Matrix<T3>& M2); // Функция конкатенации матриц

        template <typename T2>
        friend Matrix<T2> transpose(const Matrix<T2>& M);

        friend Matrix<long double> invert(const Matrix<T>& M); // Обратная матрица

        friend Matrix<T> power(const Matrix<T>& M, int p); // Возведение матрицы в степень

        template < typename T2>
        friend Matrix<long double> solve(const Matrix<T2>& M1, const Matrix<T2>& M2); // Решение системы линейных уравнений

        template <typename T2>
        friend Matrix<T2> ToDiagonal(Matrix<T2>& M);
        //friend Matrix<T> load_matrix(const char* filename);//загрузка матрицы с файла
        ~Matrix() {
            delete[] m_ptr;
        }
    };

    template<typename T2> class Matrix;

    template<typename T2>
    Matrix<long double> invert(const Matrix<T2>& M);

    template<typename T2>
    Matrix<long double> power(const Matrix<T2>& M, int p);

    template<typename T1, typename T2>
    Matrix<long double> solve(const Matrix<T1>& A, const Matrix<T2>& f);

    template <typename T2, typename T3>
    Matrix<common_type_t<T3, T2>> concatenate(const Matrix<T2>& M1, const Matrix<T3>& M2);

    template <typename T, typename T2>
    bool operator==(const Matrix<T>& M1, const Matrix<T2>& M2) noexcept;

    template <typename T, typename T2>
    bool operator!=(const Matrix<T>& M1, const Matrix<T2>& M2) noexcept;

}

#endif
