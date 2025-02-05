#ifndef MFP_FLUX_H
#define MFP_FLUX_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"

using namespace amrex;

class State;

//================================================================================

class RiemannSolver
{
public:
    RiemannSolver();
    virtual ~RiemannSolver();

    int idx;
    State* istate;

    virtual void solve(Vector<Real> &L,
                               Vector<Real> &R,
                               Vector<Real> &F,
                               Real* shk);

    virtual const std::string& get_tag() const = 0;
    virtual int get_flux_size() const = 0;
    virtual bool requires_shock_detector() const {return false;}
};

template <typename D>
std::unique_ptr<RiemannSolver> RiemannSolverBuilder(const sol::table& def)
{
    int global_idx = def["global_idx"];

    if ((def["flux"] == D::tag) && (D::valid_state(global_idx))) {
        return std::unique_ptr<D>(new D(global_idx));
    } else {
        return nullptr;
    }
}

PhysicsFactory<RiemannSolver>& GetRiemannSolverFactory();

#endif // MFP_FLUX_H
