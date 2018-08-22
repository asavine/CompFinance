#pragma once

#include <vector>
using namespace std;

//  Simple matrix class that wraps a vector,
//  See chapters 1 and 2

template <class T>
class matrix
{
    size_t      myRows;
    size_t      myCols;
    vector<T>   myVector;

public:

    //  Constructors
    matrix() : myRows(0), myCols(0) {}
    matrix(const size_t rows, const size_t cols) : myRows(rows), myCols(cols), myVector(rows*cols) {}

    //  Copy, assign
    matrix(const matrix& rhs) : myRows(rhs.myRows), myCols(rhs.myCols), myVector(rhs.myVector) {}
    matrix& operator=(const matrix& rhs)
    {
        if (this == &rhs) return *this;
        matrix<T> temp(rhs);
        swap(temp);
        return *this;
    }

    //  Copy, assign from different (convertible) type
    template <class U>
    matrix(const matrix<U>& rhs)
        : myRows(rhs.rows()), myCols(rhs.cols())
    {
        myVector.resize(rhs.rows() * rhs.cols());
        copy(rhs.begin(), rhs.end(), myVector.begin());
    }
    template <class U>
    matrix& operator=(const matrix<U>& rhs)
    {
        if (this == &rhs) return *this;
        matrix<T> temp(rhs);
        swap(temp);
        return *this;
    }

    //  Move, move assign
    matrix(matrix&& rhs) : myRows(rhs.myRows), myCols(rhs.myCols), myVector(move(rhs.myVector)) {}
    matrix& operator=(matrix&& rhs)
    {
        if (this == &rhs) return *this;
        matrix<T> temp(move(rhs));
        swap(temp);
        return *this;
    }

    //  Swapper
    void swap(matrix& rhs)
    {
        myVector.swap(rhs.myVector);
        ::swap(myRows, rhs.myRows);
        ::swap(myCols, rhs.myCols);
    }

    //  Resizer
    void resize(const size_t rows, const size_t cols)
    {
        myRows = rows;
        myCols = cols;
        if (myVector.size() < rows*cols) myVector = vector<T>(rows*cols);
    }

    //  Access
    size_t rows() const { return myRows; }
    size_t cols() const { return myCols; }
    //  So we can call matrix [i][j]
    T* operator[] (const size_t row) { return &myVector[row*myCols]; }
    const T* operator[] (const size_t row) const { return &myVector[row*myCols]; }
    bool empty() const { return myVector.empty(); }

    //  Iterators
    typedef typename vector<T>::iterator iterator;
    typedef typename vector<T>::const_iterator const_iterator;
    iterator begin() { return myVector.begin(); }
    iterator end() { return myVector.end(); }
    const_iterator begin() const { return myVector.begin(); }
    const_iterator end() const { return myVector.end(); }
};

template <class T>
inline matrix<T> transpose(const matrix<T>& mat)
{
    matrix<T> res(mat.cols(), mat.rows());
    for (size_t i = 0; i < res.rows(); ++i)
    {
        for (size_t j = 0; j < res.cols(); ++j)
        {
            res[i][j] = mat[j][i];
        }
    }

    return res;
}
