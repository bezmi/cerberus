#include "MFP_hllc.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================
std::string HydroHLLC::tag = "HLLC";
bool HydroHLLC::registered = GetRiemannSolverFactory().Register(HydroHLLC::tag, RiemannSolverBuilder<HydroHLLC>);

HydroHLLC::HydroHLLC(){}
HydroHLLC::HydroHLLC(const int i)
{
    idx = i;
    istate = GD::get_state_ptr(i);
}

void HydroHLLC::solve(Vector<Real> &L,
                      Vector<Real> &R,
                      Vector<Real> &F,
                      Real* shk)
{
    BL_PROFILE("HydroHLLC::solve");

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroState::FluxIdx::Density];
    Real uL = L[+HydroState::FluxIdx::Xvel];
    Real vL = L[+HydroState::FluxIdx::Yvel];
    Real wL = L[+HydroState::FluxIdx::Zvel];
    Real pL = L[+HydroState::FluxIdx::Prs];
    Real gamL = L[+HydroState::FluxIdx::Gamma];
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);

    Vector<Real> tLs;
    Real apL;
    for (int i = +HydroState::FluxIdx::Density + 1; i < L.size(); ++i) {
        apL = L[i];
        tLs.push_back(apL * rhoL);

    }

    Real aL = sqrt(gamL*pL/rhoL);

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroState::FluxIdx::Density];
    Real uR = R[+HydroState::FluxIdx::Xvel];
    Real vR = R[+HydroState::FluxIdx::Yvel];
    Real wR = R[+HydroState::FluxIdx::Zvel];
    Real pR = R[+HydroState::FluxIdx::Prs];
    Real gamR = R[+HydroState::FluxIdx::Gamma];
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real aR = sqrt(gamR*pR/rhoR);

    Vector<Real> tRs;
    Real apR;
    for (int i = +HydroState::FluxIdx::Density + 1; i < R.size(); ++i) {
        apR = R[i];

        tRs.push_back(apR*rhoR);

    }

    // Calculate wave speeds S_L, S_star and S_R

    Real rho_bar = 0.5*(rhoL + rhoR);
    Real a_bar = 0.5*(aL + aR);

    Real p_star = 0.5*(pL + pR) - 0.5*(uR - uL)*rho_bar*a_bar;
    Real S_star = 0.5*(uL + uR) - 0.5*(pR - pL)/(rho_bar*a_bar);

    Real qL;
    if (p_star <= pL) {
        qL = 1.0;
    } else {
        qL = std::sqrt(1.0 + ((gamL+1.0)/(2*gamL))*(p_star/pL - 1.0));
    }

    Real S_L = uL - aL*qL;

    Real qR;
    if (p_star <= pR) {
        qR = 1.0;
    } else {
        qR = std::sqrt(1.0 + ((gamR+1.0)/(2*gamR))*(p_star/pR - 1.0));
    }

    Real S_R = uR + aR*qR;
    int nc = F.size();


    if (S_L >= 0.0) {
        // flux vector L
        F[+HydroState::ConsIdx::Density]  = rhoL*uL;
        F[+HydroState::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        F[+HydroState::ConsIdx::Ymom]   = rhoL*uL*vL;
        F[+HydroState::ConsIdx::Zmom]   = rhoL*uL*wL;
        F[+HydroState::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i = 0; i < tLs.size(); ++i) {
            F[+HydroState::ConsIdx::Density + 1 + i] = tLs[i] * uL;
        }
        return;
    } else if ((S_L <= 0.0) && (0.0 <= S_star)) {
        Vector<Real> svLs(nc), fvL(nc), svL(nc);
        // flux vector L
        fvL[+HydroState::ConsIdx::Density]  = rhoL*uL;
        fvL[+HydroState::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        fvL[+HydroState::ConsIdx::Ymom]   = rhoL*uL*vL;
        fvL[+HydroState::ConsIdx::Zmom]   = rhoL*uL*wL;
        fvL[+HydroState::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i = 0; i < tLs.size(); ++i) {
            fvL[+HydroState::ConsIdx::Density + 1 + i] = tLs[i] * uL;
        }

        // state vector L
        svL[+HydroState::ConsIdx::Density]  = rhoL;
        svL[+HydroState::ConsIdx::Xmom]   = rhoL*uL;
        svL[+HydroState::ConsIdx::Ymom]   = rhoL*vL;
        svL[+HydroState::ConsIdx::Zmom]   = rhoL*wL;
        svL[+HydroState::ConsIdx::Eden]   = nrgL;
        for (int i = 0; i < tLs.size(); ++i) {
            svL[+HydroState::ConsIdx::Density + 1 + i] = tLs[i];
        }

        Real coeff = rhoL*((S_L - uL)/(S_L - S_star));

        svLs[+HydroState::ConsIdx::Density] = coeff;
        svLs[+HydroState::ConsIdx::Xmom] = coeff*S_star;
        svLs[+HydroState::ConsIdx::Ymom] = coeff*vL;
        svLs[+HydroState::ConsIdx::Zmom] = coeff*wL;
        svLs[+HydroState::ConsIdx::Eden] = coeff*(nrgL/rhoL + (S_star - uL)*(S_star + pL/(rhoL*(S_L - uL))));
        for (int i = 0; i < tLs.size(); ++i) {
            svLs[+HydroState::ConsIdx::Density + 1 + i] = tLs[i]*((S_L - uL)/(S_L - S_star));
        }

        for (int i=0; i<nc; ++i) {
            F[i] = fvL[i] + S_L*(svLs[i] - svL[i]);
        }
    } else if ((S_star <= 0.0) && (0.0 <= S_R)) {
        Vector<Real> svRs(nc), fvR(nc), svR(nc);
        // flux vector R
        fvR[+HydroState::ConsIdx::Density]  = rhoR*uR;
        fvR[+HydroState::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        fvR[+HydroState::ConsIdx::Ymom]   = rhoR*uR*vR;
        fvR[+HydroState::ConsIdx::Zmom]   = rhoR*uR*wR;
        fvR[+HydroState::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i = 0; i < tRs.size(); ++i) {
            fvR[+HydroState::ConsIdx::Density + 1 + i] = tRs[i] * uR;
        }

        // state vector R
        svR[+HydroState::ConsIdx::Density]  = rhoR;
        svR[+HydroState::ConsIdx::Xmom]   = rhoR*uR;
        svR[+HydroState::ConsIdx::Ymom]   = rhoR*vR;
        svR[+HydroState::ConsIdx::Zmom]   = rhoR*wR;
        svR[+HydroState::ConsIdx::Eden]   = nrgR;
        for (int i = 0; i < tRs.size(); ++i) {
            svR[+HydroState::ConsIdx::Density + 1 + i] = tRs[i];
        }

        Real coeff = rhoR*((S_R - uR)/(S_R - S_star));

        svRs[+HydroState::ConsIdx::Density] = coeff;
        svRs[+HydroState::ConsIdx::Xmom] = coeff*S_star;
        svRs[+HydroState::ConsIdx::Ymom] = coeff*vR;
        svRs[+HydroState::ConsIdx::Zmom] = coeff*wR;
        svRs[+HydroState::ConsIdx::Eden] = coeff*(nrgR/rhoR + (S_star - uR)*(S_star + pR/(rhoR*(S_R - uR))));
        for (int i = 0; i < tLs.size(); ++i) {
            svRs[+HydroState::ConsIdx::Density + 1 + i] = tRs[i]*((S_R - uR)/(S_R - S_star));
        }

        for (int i=0; i<nc; ++i) {
            F[i] = fvR[i] + S_R*(svRs[i] - svR[i]);
        }

    } else {
        // flux vector R
        F[+HydroState::ConsIdx::Density]  = rhoR*uR;
        F[+HydroState::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        F[+HydroState::ConsIdx::Ymom]   = rhoR*uR*vR;
        F[+HydroState::ConsIdx::Zmom]   = rhoR*uR*wR;
        F[+HydroState::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i = 0; i < tRs.size(); ++i) {
            F[+HydroState::ConsIdx::Density + 1 + i] = tRs[i] * uR;
        }
    }
}

bool HydroHLLC::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}
