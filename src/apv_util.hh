#ifndef APV_UTIL_HH
#define APV_UTIL_HH

#include <cstddef>

#ifdef _OPENMP
  #include "omp.h"
#endif

struct apv_threads
{
  apv_threads() {}
  ~apv_threads() {}

#ifdef _OPENMP
  static int thread_id() {return omp_get_thread_num();}
  static int numthreads() {return omp_get_max_threads();}
  static double tic() {return omp_get_wtime();}
  static double toc(double t0_) {return omp_get_wtime()-t0_;}
#else
  static int thread_id() {return 0;}
  static int numthreads() {return 1;}
  static double tic() {return 0.0;}
  static double toc(double t0_) {return 0.0;}
#endif
};

struct apv_chunk
{
  apv_chunk(const bool owner_,size_t size_): owner(owner_), size(size_) {}
  ~apv_chunk() {}

  const bool owner;
  const size_t size;

};
  struct apv_ichunk : public apv_chunk
  {
    apv_ichunk(bool owner_, size_t size_, int *chunk_): apv_chunk(owner_,size_), chunk(chunk_) {}

      apv_ichunk(size_t size_): apv_ichunk(true,size_,new int[size_]) {}
      apv_ichunk(size_t size_,int *chunk_): apv_ichunk(false,size_,chunk_) {}

    ~apv_ichunk() {if (owner) delete [] chunk;}

    int * const chunk;
  };
  struct apv_dchunk : public apv_chunk
  {
    apv_dchunk(bool owner_, size_t size_, double *chunk_): apv_chunk(owner_,size_), chunk(chunk_) {}

      apv_dchunk(size_t size_): apv_dchunk(true,size_,new double[size_]) {}
      apv_dchunk(size_t size_,double *chunk_): apv_dchunk(false,size_,chunk_) {}

    ~apv_dchunk() {if (owner) delete [] chunk;}

    double * const chunk;
  };


struct apv_matrix
{
  apv_matrix(bool chunk_owner_, size_t rows_, size_t cols_):
    chunk_owner(chunk_owner_), rows(rows_), cols(cols_) {}
  ~apv_matrix() {}

  const bool chunk_owner;
  const size_t  rows,
                cols;
};
  struct apv_imatrix : public apv_matrix
  {

    apv_imatrix(bool chunk_owner_, size_t rows_, size_t cols_, apv_ichunk * const ichunk_):
      apv_matrix(chunk_owner_,rows_,cols_), ichunk(ichunk_) {}

      apv_imatrix(size_t rows_, size_t cols_): apv_imatrix(true,rows_,cols_,new apv_ichunk(rows_*cols_)) {}
      apv_imatrix(size_t rows_, size_t cols_,apv_ichunk *ichunk_):
        apv_imatrix(false,(ichunk_->size)/cols_,(ichunk_->size)/rows_,ichunk_) {}

    ~apv_imatrix() {if (chunk_owner) delete ichunk;}

    apv_ichunk * const ichunk;

  };
  struct apv_dmatrix : public apv_matrix
  {

    apv_dmatrix(bool chunk_owner_, size_t rows_, size_t cols_, apv_dchunk * const dchunk_):
      apv_matrix(chunk_owner_,rows_,cols_), dchunk(dchunk_) {}

      apv_dmatrix(size_t rows_, size_t cols_): apv_dmatrix(true,rows_,cols_,new apv_dchunk(rows_*cols_)) {}
      apv_dmatrix(size_t rows_, size_t cols_,apv_dchunk *dchunk_):
        apv_dmatrix(false,(dchunk_->size)/cols_,(dchunk_->size)/rows_,dchunk_) {}

    ~apv_dmatrix() {if (chunk_owner) delete dchunk;}

    apv_dchunk * const dchunk;

  };


template <typename T> T ** Tmatrix(int M_, int N_)
{
  T * chunk = new T[M_*N_],
    ** rows = new T*[M_];
  for (int i = 0,j=0; i < M_; i++,j+=N_)
    rows[i] = chunk+j;
  return rows;
}
template <typename T> void free_Tmatrix(T ** Tmat_)
{
  delete [] Tmat_[0];
  delete [] Tmat_;
}
template <typename T> T *** T3tensor(int L_, int M_, int N_)
{
  T * chunk = new T[L_*M_*N_],
    ** rows = new T*[L_*M_],
    *** mats = new T**[L_];
  for (size_t l = 0; l < L_; l++)
  {
    mats[l] = rows + l*M_;
    for (size_t m = 0; m < M_; m++)
      mats[l][m] = chunk + (l*M_*N_) + (m*N_);
  }
  return mats;
}
template <typename T> void free_T3tensor(T *** Ttens_)
{
  delete [] Ttens_[0][0];
  delete [] Ttens_[0];
  delete [] Ttens_;
}
template <typename T> T **** T4tensor(int K_, int L_, int M_, int N_)
{
  T * chunk = new T[K_*L_*M_*N_],
    ** rows = new T*[K_*L_*M_],
    *** mats = new T**[K_*L_],
    **** tsrs = new T***[K_];
  for (size_t k = 0; k < K_; k++)
  {
    tsrs[k] = mats + (k*L_);
    for (size_t l = 0; l < L_; l++)
    {
      tsrs[k][l] = rows + (k*L_*M_) + (l*M_);
      for (size_t m = 0; m < M_; m++)
        tsrs[k][l][m] = chunk + (k*L_*M_*N_) + (l*M_*N_) + (m*N_);
    }
  }
  return tsrs;
}
template <typename T> void free_T4tensor(T **** Ttens_)
{
  delete [] Ttens_[0][0][0];
  delete [] Ttens_[0][0];
  delete [] Ttens_[0];
  delete [] Ttens_;
}
template <typename T> T ** Tsym(int M_)
{
  T * chunk = new T[((M_+1)*(M_))/2],
    ** rows = new T*[M_];
  rows[0] = chunk;
  for (int i = 1, j=M_; i < M_; i++,j--) rows[i] = rows[i-1]+j;
  return rows;
}
template <typename T> T ** Tsym_lower(int M_)
{
  T * chunk = new T[((M_+1)*(M_))/2],
    ** rows = new T*[M_];
  rows[0] = chunk;
  for (int i = 1; i < M_; i++) rows[i] = rows[i-1]+i;
  return rows;
}

#endif
