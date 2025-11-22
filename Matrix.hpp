#include"Matrix.h"
#include <iostream>
#include <stdexcept>
#include <initializer_list>
#include <iomanip>

using namespace std;
namespace linalg
{
    //--------------------КОНСТРУКТОРЫ-------------------//
    template <typename T>
    Matrix<T>::Matrix(int rows, int columns)
    {
        m_columns = columns;
        m_rows = rows;
        m_capacity = columns * rows;
        m_ptr = new T[rows * columns];
    }
    template <typename T>
    Matrix<T>::Matrix(int rows) : Matrix(rows, 1) {}

    template <typename T>
    Matrix<T>::Matrix(const Matrix<T>& M) : m_rows(M.m_rows), m_columns(M.m_columns), m_capacity(M.m_capacity) {
        m_ptr = new T[m_rows * m_columns];
        for (int i = 0; i < m_rows * m_columns; ++i) {
            m_ptr[i] = M.m_ptr[i];
        }
    }


    template <typename T>
    template <typename T2>
    Matrix<T>::Matrix(const Matrix<T2>& M)
    {
        if (M.empty())
        {
            throw runtime_error("\nError! Matrix is empty.\n");
        }
        m_rows = M.rows();
        m_columns = M.columns();
        m_capacity = m_rows * m_columns;
        m_ptr = new T[m_capacity]; // Выделение новой памяти
        for (int i = 0; i < m_capacity; i++)
        {
            m_ptr[i] = static_cast<T>(M(i / M.columns(), i % M.columns()));
        }
    }




    template <typename T>
    template <typename T2>
    Matrix<T>::Matrix(initializer_list<T2> lst)
    {
        this->m_rows = lst.size();
        this->m_columns = 1;
        this->m_capacity = this->m_rows;
        this->m_ptr = new T[m_rows];
        for (int i = 0; i < m_rows; i++)
        {
            try
            {
                m_ptr[i] = static_cast<T>(*(lst.begin() + i));
            }
            catch (const exception& e)
            {
                throw runtime_error("Can't use this type");
            }
        }

    }
    template <typename T>
    template <typename T2>
    Matrix<T>::Matrix(initializer_list<initializer_list<T2>> lst)
    {
        int tmp = 0;
        auto check = *lst.begin();
        this->m_rows = (int)lst.size();
        for (auto i : lst)
        {
            if (i.size() != check.size())
            {
                throw runtime_error("\nError! Wrong sizes of matrix.\n");
            }
            this->m_columns = (int)i.size();
        }
        this->m_ptr = new T[this->m_rows * this->m_columns];
        auto it = lst.begin();

        for (int i = 0; i < m_rows * m_columns; i++)
        {
            m_ptr[i] = static_cast<T>(*(it->begin() + i - tmp));
            if ((i + 1) % m_columns == 0)
            {
                tmp += m_columns;
                it++;
            }
        }
        this->m_capacity = this->m_rows * this->m_columns;

    }

    template <typename T>
    Matrix<T>::Matrix(Matrix<T>&& moved) noexcept : m_rows(moved.m_rows), m_columns(moved.m_columns), m_ptr(moved.m_ptr), m_capacity(moved.m_capacity) {
        moved.m_ptr = nullptr;
        moved.m_rows = 0;
        moved.m_columns = 0;
        moved.m_capacity = 0;
    }


    template <typename T>
    void Matrix<T>::shrink_to_fit()
    {
        this->m_capacity = this->m_rows * this->m_columns;
        T* new_ptr = static_cast<T*>(realloc(this->m_ptr, this->m_capacity * sizeof(T)));
        if (new_ptr) {
            this->m_ptr = new_ptr;
        }
        else {

            throw std::bad_alloc();
        }
    }

    //--------------------------------------ОПЕРАТОРЫ-----------------------------------------//
    template <typename T>
    Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other)
    {
        if (this != &other) {
            if (other.m_rows == 0 || other.m_columns == 0) {
                throw std::runtime_error("\n\nError! Can't assign empty object.\n");
            }

            if (m_capacity != other.m_capacity) {
                delete[] m_ptr;
                m_ptr = new T[other.m_capacity];
                m_capacity = other.m_capacity;
            }

            m_rows = other.m_rows;
            m_columns = other.m_columns;
            for (int i = 0; i < m_capacity; ++i) {
                m_ptr[i] = other.m_ptr[i];
            }
        }
        return *this;
    }
    template<typename T>
    template<typename T2>
    Matrix<T>& Matrix<T>::operator=(const Matrix<T2>& M)
    {

        this->m_rows = M.rows();
        this->m_columns = M.columns();
        if (M.capacity() > this->capacity())
        {
            this->m_ptr = new T[this->m_rows * this->m_columns];
        }
        for (int i = 0; i < (this->m_rows * this->m_columns); i++)
        {
            this->m_ptr[i] = static_cast<T>(M(i / M.columns(), i % M.columns()));
        }
        this->shrink_to_fit();
        return *this;
    }


    template<typename T>
    template<typename T2>
    Matrix<T>& Matrix<T>::operator=(Matrix<T2>&& M) noexcept
    {
        // Освобождаем текущую память
        delete[] this->m_ptr;

        // Перемещаем данные
        this->m_rows = M.m_rows;
        this->m_columns = M.m_columns;
        this->m_capacity = M.m_capacity;

        // Создаем новый массив с преобразованием типов
        this->m_ptr = new T[this->m_rows * this->m_columns];
        for (int i = 0; i < (this->m_rows * this->m_columns); i++) {
            this->m_ptr[i] = static_cast<T>(M.m_ptr[i]);
        }

        // Обнуляем данные перемещаемого объекта
        M.m_rows = 0;
        M.m_columns = 0;
        M.m_capacity = 0;
        delete[] M.m_ptr;
        M.m_ptr = nullptr;

        return *this;
    }

    template<typename T>
    T& Matrix<T>::operator()(int i, int j) {
        return m_ptr[i * m_columns + j];
    }
    template<typename T>
    T Matrix<T>::operator()(int i, int j) const {
        return m_ptr[i * m_columns + j];
    }
    //-------------------------оператор вывода----------------------------
    template <typename T>
    ostream& operator<<(ostream& out, const Matrix<T>& M) {
        T max_element = M.m_ptr[0];
        for (int i = 0; i < M.m_rows; i++)
        {
            for (int j = 0; j < M.m_columns; j++)
                max_element = max(max_element, M.m_ptr[i * M.m_columns + j]);
        }
        int element_width = (int)to_string(max_element).length() + 1;
        for (int i = 0; i < M.m_rows; i++)
        {
            for (int j = 0; j < M.m_columns; j++)
            {
                out << setw(element_width) << M.m_ptr[i * M.m_columns + j] << " ";
            }
            out << endl;
        }
        return out;
    }

    template <typename T>
    template <typename T2>
    Matrix<common_type_t<T, T2>> Matrix<T>::operator+(const Matrix<T2>& M) const
    {
        if ((this->m_rows == M.rows()) && (this->m_columns == M.columns()))
        {
            Matrix<common_type_t<T, T2>> result(m_rows, m_columns);
            for (int i = 0; i < m_capacity; ++i)
                result.m_ptr[i] = m_ptr[i] + M.m_ptr[i];

            return result;
        }
        else throw runtime_error("impossible operation");
    }

    template<typename T>
    template<typename T2>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T2>& M)
    {
        if ((this->m_columns == M.m_columns) && (this->m_rows == M.m_rows))
        {
            for (int i = 0; i < M.m_rows * M.m_columns; i++)
            {
                m_ptr[i] += static_cast<T>(M(i / M.columns(), i % M.columns()));
            }
            return *this;

        }
        else
            throw runtime_error("impossible operation");
    }


    template <typename T>
    template <typename T2>
    Matrix<common_type_t<T, T2>> Matrix<T>::operator - (const Matrix<T2>& M) const
    {
        if (((*this).m_columns == M.m_columns) && ((*this).m_rows == M.m_rows))
        {
            Matrix<common_type_t<T, T2>> result(m_rows, m_columns);
            for (int i = 0; i < m_capacity; ++i)
                result.m_ptr[i] = m_ptr[i] - M.m_ptr[i];
            return result;
        }
        else
            throw runtime_error("impossible operation");
    }


    template<typename T>
    template<typename T2>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T2>& M)
    {
        if ((this->m_columns == M.m_columns) && (this->m_rows == M.m_rows))
        {
            for (int i = 0; i < M.m_rows * M.m_columns; i++)
            {
                m_ptr[i] -= static_cast<T>(M(i / M.columns(), i % M.columns()));
            }
            return *this;

        }
        else
            throw runtime_error("impossible operation");
    }

    template <typename T>
    template <typename T2>
    Matrix<decltype(T{} *T2{}) > Matrix<T>::operator*(const Matrix<T2>& other) const
    {
        if (this->m_columns != other.rows())
        {
            throw runtime_error("\nError! Matrices are not competible");
        }
        Matrix<decltype(T{} *T2{}) > tmp(this->m_rows, other.columns());
        int counter = -1;
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++)
        {
            tmp.m_ptr[i] = 0;
            if (!(i % tmp.m_columns))
                counter++;
            for (int j = 0; j < this->m_columns; j++)
            {
                tmp.m_ptr[i] += this->m_ptr[counter * this->m_columns + j] * other(j, i % other.columns());
            }
        }
        return tmp;
    }

    template <typename T>
    Matrix<decltype(double{} *T{}) > operator*(double x, const Matrix<T>& mat)
    {
        Matrix<decltype(double{} *T{}) > tmp(mat.rows(), mat.columns());
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++)
        {
            tmp.m_ptr[i] = x * mat.m_ptr[i];
        }
        return tmp;
    }

    template <typename T>
    Matrix<decltype(T{} *double{}) > Matrix<T>::operator*(double x) const
    {
        Matrix<decltype(T{} *double{}) > tmp(this->rows(), this->columns());
        for (int i = 0; i < tmp.columns() * tmp.rows(); i++)
        {
            tmp(i / tmp.columns(), i % tmp.columns()) = this->m_ptr[i] * x;
        }
        return tmp;
    }

    template <typename T>
    template <typename T2>
    Matrix<T>& Matrix<T>::operator*=(const Matrix<T2>& other)
    {
        if (this->m_columns != other.rows())
        {
            throw runtime_error("\nError! Matrices are not competible");
        }

        Matrix<T> result(this->m_rows, other.columns());

        for (int i = 0; i < this->m_rows; i++) {
            for (int j = 0; j < other.columns(); j++)
            {
                result(i, j) = 0;
                for (int k = 0; k < this->m_columns; k++)
                {
                    result(i, j) = result(i, j) + (*this)(i, k) * other(k, j);

                }
            }
        }
        *this = result;
        return *this;
    }


    template <typename T>
    Matrix<T>& Matrix<T>::operator*=(double x) noexcept
    {
        for (int i = 0; i < this->m_columns * this->m_rows; i++)
        {
            this->m_ptr[i] = this->m_ptr[i] * static_cast<T>(x);
        }
        return *this;
    }






    //--------------------------Norm---------------------------//
    template<typename T>
    double Matrix<T>::norm() {
        double sum = 0.0;
        for (const auto& elem : m_ptr) {
            sum += elem * elem;
        }
        return sqrt(sum);
    }

    //--------------------------Trace-------------------------//
    template<typename T>
    T Matrix<T>::trace() const {
        if (m_rows != m_columns) {
            throw invalid_argument("Matrix must be square.");
        }
        T sum = 0;
        for (int i = 0; i < m_rows; ++i) {
            sum += m_ptr[i * m_columns + i];
        }
        return sum;
    }

    //--------------------------Det---------------------------//
    template <typename T>
    T Matrix<T>::det()const
    {
        if (m_rows == m_columns) {
            for (int i = 0; i < m_columns; i++) {
                int nonZeroIndex = -1;
                for (int j = i; j < m_rows; j++) {
                    if ((*this)(j, i) != 0) {
                        nonZeroIndex = j;
                        break;
                    }
                }
                if (nonZeroIndex != -1) {
                    if (nonZeroIndex != i) {
                        for (int k = 0; k < m_columns; k++) {
                            swap((*this)(i, k), (*this)(nonZeroIndex, k));
                        }
                    }
                    double diagonalElement = (*this)(i, i);
                    for (int j = i + 1; j < m_rows; j++) {
                        double factor = (*this)(j, i) / diagonalElement;
                        for (int k = 0; k < m_columns; k++) {
                            (*this)(j, k) -= (*this)(i, k) * factor;
                        }
                    }
                }
            }
            double determinant = 1;
            for (int i = 0; i < m_rows; i++) {
                determinant *= (*this)(i, i);
            }
            return determinant;
        }
        else {
            throw runtime_error("impossible operation");
        }
    }
    //-------------------------GAUSS--------------------------//
    template <typename T>
    Matrix<T>& Matrix<T>::gauss_forward()
    {
        for (int i = 0; i < m_columns; i++)
        {
            // ѕоиск ненулевого элемента в текущем столбце, начина€ с текущей строки
            int nonzero_row = -1; // »ндекс строки с ненулевым элементом
            for (int j = i; j < m_rows; j++)
            {
                if ((*this)(j, i) != 0)
                {
                    nonzero_row = j;
                    break;
                }
            }

            if (nonzero_row != -1)
            {
                // ≈сли найден ненулевой элемент, мен€ем местами текущую строку с строкой, содержащей ненулевой элемент
                for (int k = 0; k < m_columns; k++)
                {
                    double temp = (*this)(i, k);
                    (*this)(i, k) = (*this)(nonzero_row, k);
                    (*this)(nonzero_row, k) = temp;
                }

                // ƒелаем текущий элемент равным 1 путем делени€ всей строки на него
                double number = (*this)(i, i);
                for (int j = 0; j < m_columns; j++)
                    (*this)(i, j) /= number;
                // ¬ычитаем текущую строку, умноженную на число, из всех следующих строк
                for (int j = i + 1; j < m_rows; j++)
                {
                    number = (*this)(j, i);
                    for (int k = 0; k < m_columns; k++)
                        (*this)(j, k) -= (*this)(i, k) * number;

                }
            }
        }
        return (*this);
    }
    //--------------------------------GAUSS--------------------------------//
    template <typename T>
    Matrix<T>& Matrix<T>::gauss_backward()
    {

        // мы начинаем с последней строки и идем вверх
        for (int i = m_columns - 1; i >= 0; i--) {

            // ƒл€ каждой строки ищем не нулевой элемент в этом столбце
            for (int j = 0; j < m_rows; j++)
            {
                if ((*this)(j, i) != 0)
                {
                    double temp;
                    // ћен€ем строки местами
                    for (int k = m_columns - 1; k >= 0; k--)
                    {
                        temp = (*this)(i, k);
                        (*this)(i, k) = (*this)(j, k);
                        (*this)(j, k) = temp;
                    }
                    break;
                }
            }

            // ≈сли элемент по диагонали не равен нулю
            if ((*this)(i, i) != 0) {
                double number = (*this)(i, i);

                // ƒелим все элементы строки на диагональный дл€ получени€ единицы в верхнетреугольной матрице
                for (int j = 0; j < m_columns; j++)
                    (*this)(i, j) /= number;

                // ¬ычитаем из верхних строк текущую строку (умноженную на соответствующий элемент), чтобы занулить остальные элементы в текущем столбце
                for (int j = i - 1; j >= 0; j--) {
                    number = (*this)(j, i);
                    for (int k = m_columns - 1; k >= 0; k--)
                        (*this)(j, k) -= (*this)(i, k) * number;
                }
            }
        }
        return (*this);
    }
    //-------------------------RANK-------------------------------//
    template <typename T>
    int Matrix<T>::rank() const
    {
        Matrix<T> Mcopy = *this;

        // ѕримен€ем метод √аусса
        Mcopy.gauss_backward();

        // —читаем количество ненулевых строк из преобразованной матрицы
        int rank = 0;
        for (int i = 0; i < Mcopy.m_rows; ++i) {
            bool allZeros = true;
            for (int j = 0; j < Mcopy.m_columns; ++j) {
                if (Mcopy(i, j) != 0.0) { // ѕредполагаетс€, что метод (i, j) возвращает значение €чейки
                    allZeros = false;
                    break;
                }
            }
            if (!allZeros) {
                rank++;
            }
        }

        return rank;
    }

    //-------------------------соединить матриц------------------//
    template <typename T>
    Matrix<T> concatenate(const Matrix<T>& M1, const Matrix<T>& M2)
    {
        if (M1.m_rows != M2.m_rows)
        {
            throw runtime_error("\nError! Matrices dont have the appropriate sizes.\n");
        }
        Matrix<T> tmp(M1.m_rows, M1.m_columns + M2.m_columns);
        int flag = 0, k = 1;
        tmp.m_ptr[0] = M1.m_ptr[0];
        for (int i = 1; i < tmp.m_columns * tmp.m_rows; i++)
        {
            if (k % M1.m_columns == 0 && (flag % 2) == 0)
            {
                k = 0;
                flag++;
            }
            else if (k % M2.m_columns == 0 && (flag % 2) != 0)
            {
                k = 0;
                flag++;
            }
            if (!(flag % 2))
            {
                tmp.m_ptr[i] = M1.m_ptr[k + (flag / 2) * M1.m_columns];
                k++;
            }
            else if (flag % 2)
            {
                tmp.m_ptr[i] = M2.m_ptr[k + ((flag - 1) / 2) * M2.m_columns];
                k++;
            }
        }
        return tmp;
    }
    //--------------------------транспонирование-------------------//
    template <typename T>
    Matrix<T> transpose(const Matrix<T>& M)

    {
        Matrix<T> res(M.columns(), M.rows());
        for (int i = 0; i < M.rows(); i++) {
            for (int j = 0; j < M.columns(); j++) {
                res(j, i) = M(i, j);
            }
        }
        return res;
    }
    //--------------------------Обратная матрица-------------------//
    template <typename T>
    Matrix<T> invert(const Matrix<T>& m)

    {
        int size = m.rows();
        // Проверяем, является ли матрица квадратной
        if (size != m.columns())
            throw std::runtime_error("Matrix must be square");

        // Создаем единичную матрицу того же размера
        Matrix<T> I(size, size);
        for (int i = 0; i < size; ++i)
            I(i, i) = 1.0;

        for (int i = 0; i < size; ++i) {
            double diag = m(i, i);
            for (int j = 0; j < size; ++j) {
                m(i, j) /= diag;
                I(i, j) /= diag;
            }

            for (int j = 0; j < size; ++j) {
                if (i != j) {
                    double ratio = m(j, i);

                    for (int k = 0; k < size; ++k) {
                        m(j, k) -= ratio * m(i, k);
                        I(j, k) -= ratio * I(i, k);
                    }
                }
            }
        }
        return I;
    }
    //------------------------RESHAPE----------------------------//
    template <typename T>
    void Matrix<T>::reshape(int rows, int columns)
    {
        if (this->m_rows * this->m_columns != rows * columns)
        {
            throw runtime_error("\nError Not correct dimensions.\n");
        }
        this->m_rows = rows;
        this->m_columns = columns;

        cout << *this;
    }

    //------------------------POWer------------------------------//

    template <typename T>
    Matrix<T> power(const Matrix<T>& m, int p)

    {
        if (m.rows() != m.columns()) {
            throw runtime_error("Operation is impossible");
        }

        Matrix<T> result(m.rows(), m.columns());
        if (p < 0) {
            m = invert(m);
            p = -p;
        }

        for (int i = 0; i < m.rows(); ++i) {
            result(i, i) = 1;
        }

        while (p) {
            if (p % 2 == 1)
                result = result * m;
            m = m * m;
            p >>= 1;
        }

        return result;
    }
    //------------------------SOLVE------------------------------//
    template <typename T>

    Matrix<long double> solve(const Matrix<T>& A, const Matrix<T>& f)
    {
        if (A.rows() != f.rows() || f.columns() > 1)
        {
            throw runtime_error("Invalid sizes of the input matrices.");
        }

        Matrix<T> P = concatenate(A, f);
        P.gauss_forward();

        Matrix<T> tmp_matr(A.rows(), 1);

        for (int i = 0; i < A.rows(); ++i) {

            tmp_matr(i, 0) = P(i, A.columns());

        }


        return tmp_matr;
    }
    template <typename T>

    Matrix<T> ToDiagonal(Matrix<T>& M)
    {
        if (M.rows() == M.columns())
        {
            Matrix<T> result = M;
            result = transpose(result.gauss_forward());
            return result.gauss_forward();

        }
        else
            throw runtime_error("Invalid sizes of the input matrices.");

    }


    template <typename T>
    void Matrix<T>::reserve(int n)
    {
        this->m_ptr = (T*)realloc(this->m_ptr, n);
        this->m_capacity = n;
    }



    template <typename T>
    void Matrix<T>::clear()
    {
        this->m_rows = 0;
        this->m_columns = 0;
    }

}





