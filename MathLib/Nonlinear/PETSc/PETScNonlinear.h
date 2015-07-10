/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_NONLINEAR_PETSCNONLINEAR_H_
#define MATHLIB_NONLINEAR_PETSCNONLINEAR_H_

#include <petscsnes.h>

#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

namespace MathLib
{

class IVector;

template <class T>
void fromPETScVec(Vec &src, T &dest)
{
    for (PetscInt i=0; i<(PetscInt)dest.size(); i++) {
        double v = .0;
        VecGetValues(src, 1, &i, &v);
        dest[i] = v;
    }
}

template <class T>
void toPETScVec(T &src, Vec &dest)
{
    for (PetscInt i=0; i<(PetscInt)src.size(); i++) {
        VecSetValue(dest, i, src[i], INSERT_VALUES);
    }
    VecAssemblyBegin(dest);
    VecAssemblyEnd(dest);
}

template <class T>
void toPETScMat(T &src, Mat &dest)
{
    for (PetscInt i=0; i<(PetscInt)src.getNRows(); i++) {
        for (PetscInt j=0; j<(PetscInt)src.getNCols(); j++) {
            if (src(i,j)==.0) continue;
            MatSetValue(dest, i, j, src(i,j), INSERT_VALUES);
        }
    }
    MatAssemblyBegin(dest, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dest, MAT_FINAL_ASSEMBLY);
}

template <class F_R, class F_J, class T_VEC, class T_MAT>
struct Wrapper
{
    Wrapper(F_R &fr, F_J &fj, T_VEC &x) : _f_r(fr), _f_j(fj), _x(x)
    {
        _r = new T_VEC(x);
        _jac = new T_MAT(x.size(), x.size());
    }

    ~Wrapper()
    {
        delete _r;
        delete _jac;
    }

    PetscErrorCode residual(SNES, Vec x_, Vec r_,void*)
    {
        fromPETScVec(x_, _x);
        _f_r(_x, *_r);
        toPETScVec(*_r, r_);
        return 0;
    }

    PetscErrorCode jacobian(SNES, Vec x_, Mat jac_, Mat /*B*/, void*)
    {
        fromPETScVec(x_, _x);
        _f_j(_x, *_jac);
        toPETScMat(*_jac, jac_);
        return 0;
    }

    F_R& _f_r;
    F_J& _f_j;
    T_VEC& _x;
    T_VEC* _r;
    T_MAT* _jac;
};

template <class F_R, class F_J>
struct Wrapper<F_R, F_J, MathLib::PETScVector, MathLib::PETScMatrix>
{
    Wrapper(F_R &fr, F_J &fj, MathLib::PETScVector &x) : _f_r(fr), _f_j(fj) {}

    ~Wrapper() {}

    PetscErrorCode residual(SNES, Vec x_, Vec r_,void*)
    {
        PETScVector x(x_);
        PETScVector r(r_);
        _f_r(x, r);
        return 0;
    }

    PetscErrorCode jacobian(SNES, Vec x_, Mat jac_, Mat /*B*/, void*)
    {
        PETScVector x(x_);
        PETScMatrix jac(jac_);
        _f_j(x, jac);
        return 0;
    }

    F_R& _f_r;
    F_J& _f_j;
};

template <class T>
PetscErrorCode wrap_residual(SNES snes,Vec x, Vec r, void* v)
{
    auto& w = *(T*)v;
    return w.residual(snes, x, r, nullptr);
}

template <class T>
#if (PETSC_VERSION_NUMBER >= 3050)
PetscErrorCode wrap_jacobian(SNES snes, Vec x, Mat jac, Mat B, void* v)
{
#else
PetscErrorCode wrap_jacobian(SNES snes, Vec x, Mat* jac_, Mat* B_, MatStructure* str, void* v)
{
    Mat& jac(*jac_);
    Mat& B(*B_);
    *str = SAME_NONZERO_PATTERN;
#endif
    auto& w = *(T*)v;
    return w.jacobian(snes, x, jac, B, nullptr);
}

class PETScNonlinear
{
#if (PETSC_VERSION_NUMBER >= 3050)
    typedef PetscErrorCode (*F_Residual)(SNES,Vec,Vec,void*);
    typedef PetscErrorCode (*F_Jacobian)(SNES,Vec,Mat,Mat,void*);
#else
    typedef PetscErrorCode (*F_Residual)(SNES,Vec,Vec,void*);
    typedef PetscErrorCode (*F_Jacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
#endif

public:
    PETScNonlinear()
    : _A(nullptr), _b(nullptr), _x(nullptr), _snes(nullptr), _external(false)
    {
        SNESCreate(PETSC_COMM_WORLD, &_snes);
    }


    PETScNonlinear(MathLib::IMatrix* J, MathLib::IVector* r)
    : _A(static_cast<MathLib::PETScMatrix*>(J)->getRawMatrix()),
      _b(static_cast<MathLib::PETScVector*>(r)->getData()),
      _x(nullptr), _snes(nullptr), _external(true)
    {
        SNESCreate(PETSC_COMM_WORLD, &_snes);
    }

    virtual ~PETScNonlinear()
    {
        if (!_external) {
            VecDestroy(&_b);
            MatDestroy(&_A);
        }
        VecDestroy(&_x);
        SNESDestroy(&_snes);
    }

    template <class F_R, class F_J, class T_V, class T_MAT>
    bool solve(F_R &f_residual, F_J &f_jacobian, T_V &x_)
    {
        typedef Wrapper<F_R, F_J, T_V, T_MAT> WF;
        WF wf(f_residual, f_jacobian, x_);

        if (_x==nullptr) {
            VecCreate(PETSC_COMM_WORLD, &_x);
            VecSetSizes(_x, PETSC_DECIDE, x_.size());
            VecSetUp(_x);
        }
        for (std::size_t i=0; i<x_.size(); i++)
            VecSetValue(_x, i, x_[i], INSERT_VALUES);

        bool status = solve(wrap_residual<WF>, wrap_jacobian<WF>, _x, &wf);

        for (PetscInt i=0; i<(PetscInt)x_.size(); i++) {
            double v = .0;
            VecGetValues(_x, 1, &i, &v);
            x_[i] = v;
        }

        return status;
    }

    template <class F_R, class F_J>
    bool solve(F_R &f_residual, F_J &f_jacobian, MathLib::PETScVector &x_)
    {
        typedef Wrapper<F_R, F_J, MathLib::PETScVector, MathLib::PETScMatrix> WF;
        WF wf(f_residual, f_jacobian, x_);

        bool status = solve(wrap_residual<WF>, wrap_jacobian<WF>, x_.getData(), &wf);

        return status;
    }

    bool solve(F_Residual f_residual, F_Jacobian f_jacobian, Vec &x, void* ctx)
    {
        PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------\n");
        PetscPrintf(PETSC_COMM_WORLD, "*** PETSc nonlinear solver\n");

        PetscErrorCode ierr;
        if (_b == nullptr) {
            PetscInt m = 0;
            VecGetSize(x, &m);

            VecDuplicate(x, &_b);
            MatCreate(PETSC_COMM_WORLD, &_A);
            ierr = MatSetSizes(_A,PETSC_DECIDE,PETSC_DECIDE,m,m);
            CHKERRCONTINUE(ierr);
#ifdef USE_MPI
            MatSetType(_A,MATMPIAIJ);
#else
            MatSetType(_A,MATSEQAIJ);
#endif
            MatSetOption(_A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); // for MatZeroRows()
            MatSetFromOptions(_A);
    //        MatMPIAIJSetPreallocation(_A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
    //        MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
            MatSetUp(_A);
        }

        SNESSetFunction(_snes, _b, f_residual, ctx);
        SNESSetJacobian(_snes, _A, _A, f_jacobian, ctx);
        SNESLineSearch ls;
        SNESGetLineSearch(_snes, &ls);
        SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);
        SNESLineSearchSetMonitor(ls, PETSC_TRUE);
        SNESMonitorSet(_snes, SNESMonitorDefault, nullptr, nullptr);
        SNESSetFromOptions(_snes);

        // -snes_monitor -snes_converged_reason -snes_view
        SNESSolve(_snes, NULL, x);

        SNESConvergedReason reason;
        SNESGetConvergedReason(_snes, &reason);
        PetscInt iter_nlin = 0;
        SNESGetIterationNumber(_snes, &iter_nlin);
        PetscReal xnorm, fnorm, ynorm;
        SNESLineSearchGetNorms(ls, &xnorm, &fnorm, &ynorm);
        PetscPrintf(PETSC_COMM_WORLD, "\nsummary\n");
        if (1 < reason && reason < 5)
            PetscPrintf(PETSC_COMM_WORLD, "\t status    : %s (%d)\n", "converged", reason);
        else
            PetscPrintf(PETSC_COMM_WORLD, "\t status    : %s (%d)\n", "diverged", reason);
        PetscPrintf(PETSC_COMM_WORLD, "\t iterations: %d\n", iter_nlin);
        PetscPrintf(PETSC_COMM_WORLD, "\t residual  : %e\n", fnorm);

        PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------\n");

        return true;
    }

    /// set the maximum number of iterations
    void setMaxIterations(std::size_t max_itr) {_max_itr = max_itr;}

    /// set the absolute residual tolerance used by the stopping criteria
    void setAbsResidualTolerance(double abs_tol) {_r_abs_tol = abs_tol;}

    /// set the relative residual tolerance used by the stopping criteria
    void setRelResidualTolerance(double rel_tol) {_r_rel_tol = rel_tol;}

    /// set the relative solution increment tolerance used by the stopping criteria
    void setRelDxTolerance(double rel_tol) {_dx_rel_tol = rel_tol;}


protected:
    Mat _A;
    Vec _b;
    Vec _x;
    SNES _snes;

    bool _external;

    /// absolute tolerance for residual
    double _r_abs_tol;
    /// relative tolerance for residual
    double _r_rel_tol;
    /// relative tolerance for dx
    double _dx_rel_tol;
    /// the maximum allowed iteration number
    std::size_t _max_itr;
};

} // MathLib


#endif // MATHLIB_NONLINEAR_PETSCNONLINEAR_H_
