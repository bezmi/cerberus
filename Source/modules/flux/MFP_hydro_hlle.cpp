#include "MFP_hydro_hlle.H"

#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================
std::string HydroHLLE::tag = "HLLE";
bool HydroHLLE::registered = GetRiemannSolverFactory().Register(HydroHLLE::tag, RiemannSolverBuilder<HydroHLLE>);

HydroHLLE::HydroHLLE(){}
HydroHLLE::HydroHLLE(const int i)
{
    idx = i;
    istate = GD::get_state_ptr(i);

}

void HydroHLLE::solve(Vector<Real> &L,
                      Vector<Real> &R,
                      Vector<Real> &F,
                      Real* shk)
{
    BL_PROFILE("HydroHLLE::solve");

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroState::FluxIdx::Density];
    Real uL = L[+HydroState::FluxIdx::Xvel];
    Real vL = L[+HydroState::FluxIdx::Yvel];
    Real wL = L[+HydroState::FluxIdx::Zvel];
    Real pL = L[+HydroState::FluxIdx::Prs];

    Vector<Real> tLs;
    Real apL;
    for (int i = +HydroState::FluxIdx::Alpha; i < L.size(); ++i) {
        apL = L[+HydroState::FluxIdx::Alpha+i];
        tLs.push_back(apL * rhoL);

    }

    Real gamL = L[+HydroState::FluxIdx::Gamma];
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
    Real aL = sqrt(gamL*pL/rhoL);

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroState::FluxIdx::Density];
    Real uR = R[+HydroState::FluxIdx::Xvel];
    Real vR = R[+HydroState::FluxIdx::Yvel];
    Real wR = R[+HydroState::FluxIdx::Zvel];
    Real pR = R[+HydroState::FluxIdx::Prs];

    Vector<Real> tRs;
    Real apR;
    for (int i = +HydroState::FluxIdx::Alpha; i < R.size(); ++i) {
        apR = R[+HydroState::FluxIdx::Alpha];

        tRs.push_back(apR*rhoR);

    }

    Real gamR = R[+HydroState::FluxIdx::Gamma];
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real aR = sqrt(gamR*pR/rhoR);

    // speeds
    Real sL = std::min(uL - aL, uR - aR);
    Real sR = std::max(uL + aL, uR + aR);

    if (sL >= 0.0) {
        // flux vector L
        F[+HydroState::ConsIdx::Density]  = rhoL*uL;
        F[+HydroState::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        F[+HydroState::ConsIdx::Ymom]   = rhoL*uL*vL;
        F[+HydroState::ConsIdx::Zmom]   = rhoL*uL*wL;
        F[+HydroState::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i = 0; i < tLs.size(); ++i) {
            F[+HydroState::ConsIdx::Tracer + i] = tLs[i] * uL;
        }
    } else if ((sL <= 0.0) && (sR >= 0.0)) {
        int nc = F.size();
        Vector<Real> fvL(nc), fvR(nc), svL(nc), svR(nc);
        // flux vector L
        fvL[+HydroState::ConsIdx::Density]  = rhoL*uL;
        fvL[+HydroState::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        fvL[+HydroState::ConsIdx::Ymom]   = rhoL*uL*vL;
        fvL[+HydroState::ConsIdx::Zmom]   = rhoL*uL*wL;
        fvL[+HydroState::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i = 0; i < tLs.size(); ++i) {
            fvL[+HydroState::ConsIdx::Tracer + i] = tLs[i] * uL;
        }

        // flux vector R
        fvR[+HydroState::ConsIdx::Density]  = rhoR*uR;
        fvR[+HydroState::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        fvR[+HydroState::ConsIdx::Ymom]   = rhoR*uR*vR;
        fvR[+HydroState::ConsIdx::Zmom]   = rhoR*uR*wR;
        fvR[+HydroState::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i = 0; i < tRs.size(); ++i) {
            fvR[+HydroState::ConsIdx::Tracer + i] = tRs[i] * uR;
        }

        // state vector L
        svL[+HydroState::ConsIdx::Density]  = rhoL;
        svL[+HydroState::ConsIdx::Xmom]   = rhoL*uL;
        svL[+HydroState::ConsIdx::Ymom]   = rhoL*vL;
        svL[+HydroState::ConsIdx::Zmom]   = rhoL*wL;
        svL[+HydroState::ConsIdx::Eden]   = nrgL;
        for (int i = 0; i < tLs.size(); ++i) {
            svL[+HydroState::ConsIdx::Tracer + i] = tLs[i];
        }

        // state vector R
        svR[+HydroState::ConsIdx::Density]  = rhoR;
        svR[+HydroState::ConsIdx::Xmom]   = rhoR*uR;
        svR[+HydroState::ConsIdx::Ymom]   = rhoR*vR;
        svR[+HydroState::ConsIdx::Zmom]   = rhoR*wR;
        svR[+HydroState::ConsIdx::Eden]   = nrgR;
        for (int i = 0; i < tRs.size(); ++i) {
            svR[+HydroState::ConsIdx::Tracer + i] = tRs[i];
        }

        for (int i=0; i<nc; ++i) {
            F[i] = (sR*fvL[i] - sL*fvR[i] + sL*sR*(svR[i] - svL[i]))/(sR - sL);
        }
    } else {
        F[+HydroState::ConsIdx::Density]  = rhoR*uR;
        F[+HydroState::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        F[+HydroState::ConsIdx::Ymom]   = rhoR*uR*vR;
        F[+HydroState::ConsIdx::Zmom]   = rhoR*uR*wR;
        F[+HydroState::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i = 0; i < tRs.size(); ++i) {
            F[+HydroState::ConsIdx::Tracer + i] = tRs[i] * uR;
        }
    }
}

bool HydroHLLE::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}
