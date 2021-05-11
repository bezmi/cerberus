#include "MFP_hydro.H"

#include <functional>
#include <string>

#include "MFP_optional_func.H"
#include "sol.hpp"
#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_lua.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_hydro_bc.H"
#include "Eigen"

#ifdef PYTHON
#include "MFP_diagnostics.H"
#endif

using GD = GlobalData;

//------------
// Hydro State
// TODO: get rid of the Density + 1 + a indexing

// TODO: make this a map so the order is not important?
Vector<set_bc> HydroState::bc_set = {
    &set_x_vel_bc, // Xvel
    &set_y_vel_bc, // Yvel
    &set_z_vel_bc, // Zvel
    &set_scalar_bc, // Prs
    &set_scalar_bc, // Temp
    &set_scalar_bc, // Density
};

Vector<std::string>  prim_names = {"x_vel", "y_vel", "z_vel", "p", "T", "rho"};
Vector<std::string>  cons_names = {"x_mom", "y_mom", "z_mom", "nrg", "rho"};

Vector<int> HydroState::flux_vector_idx = {+HydroState::FluxIdx::Xvel};
Vector<int> HydroState::cons_vector_idx = {+HydroState::ConsIdx::Xmom};
Vector<int> HydroState::prim_vector_idx = {+HydroState::PrimIdx::Xvel};

//-----------------------------------------------------------------------------

std::string HydroState::tag = "hydro";
bool HydroState::registered = GetStateFactory().Register(HydroState::tag, StateBuilder<HydroState>);

HydroState::HydroState() : State(){}
HydroState::HydroState(const sol::table& def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}
HydroState::~HydroState(){}

void HydroState::init_from_lua()
{
    BL_PROFILE("HydroState::init_from_lua");


    sol::state& lua = GD::lua;

    const sol::table state_def = lua["states"][name];

    // TODO: maybe move this to the factory class once the state has been registered?
    GD::num_fluid += 1;


    //
    // get mass, charge, and density
    //

    // mass = {1.0, 2.0};
    // charge = {1.0, 1.0};
    // gamma = {5.0/3.0, 5.0/3.0};
    set_values(state_def["mass"], mass);
    set_values(state_def["charge"], charge);
    set_values(state_def["gamma"], gamma);

    if (any_equal(mass.begin(), mass.end(), 0.0) or any_equal(gamma.begin(), gamma.end(), 0.0))
        Abort("State: "+name+"; mass and gamma cannot be 0");

    mass_const = all_equal(mass.begin(), mass.end(), mass[0]);
    charge_const = all_equal(charge.begin(), charge.end(), charge[0]);
    gamma_const = all_equal(gamma.begin(), gamma.end(), gamma[0]);

    if ((mass.size() != charge.size()) or (mass.size() != gamma.size()) or (charge.size() != gamma.size()))
        Abort("State: "+name+"; 'mass', 'charge' and 'gamma' must have the same number of components");


    n_components = mass.size();

    if (n_components >= 2)
        is_multicomponent = true;

    n_alphas = n_components - 1;

    for (int i = 0; i < n_alphas; ++i) {
        bc_set.push_back(&set_scalar_bc);
    }

    sum_cons.resize(n_cons());

    //
    // viscous terms coefficients
    //
    set_viscosity();

    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(n_prim());

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh=0; lh<2; ++lh) {

            std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]]["fill_hydro_bc"].get_or<std::string>("outflow");
            int i_side_bc = bc_names.at(side_bc);

            // loop over the named primitives and get any custom functions
            // TODO: write function to eliminate redundant code
            bool has_invalid_f = false;
            for (auto const& [prim_name, prim_idx] : singular_prim_names_map) {

                bs.set_phys_fill_bc(ax, i_side_bc, lh, prim_idx);

                const sol::object v = state_def["bc"][dir_name[ax]][side_name[lh]][prim_name].get_or(sol::object());
                Optional3D1VFunction f = get_udf(v);

                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid())
                    Abort("State "+name+": Setting 'fill_hydro_bc = inflow' requires all primitive variables to be defined, '" + prim_name + "' is not defined");

                bs.set(ax,prim_idx,lh,f);
            }


            if (is_multicomponent) {
                // handle alpha values
                const sol::table alpha_bc_defs = state_def["bc"][dir_name[ax]][side_name[lh]][multicomp_prim_name].get_or(sol::table());

                if (i_side_bc == PhysBCType::inflow and !alpha_bc_defs.valid())
                    Abort("State "+name+": Setting 'fill_hydro_bc = inflow' requires all primitive variables to be defined, '" + multicomp_prim_name + "' is not defined");

                if (alpha_bc_defs.valid())
                    if (alpha_bc_defs.size() < n_alphas)
                        Abort("State "+name+": Length of 'alpha' bc defs must be equal to (n_components - 1)");

                for (int i = 0; i < n_alphas; ++i) {
                    bs.set_phys_fill_bc(ax, i_side_bc, lh, +HydroState::PrimIdx::Density + 1 + i);
                    Optional3D1VFunction f;
                    if (alpha_bc_defs.valid()) {
                        f = get_udf(alpha_bc_defs[i+1]);

                    }

                    bs.set(ax,+HydroState::PrimIdx::Density + 1 + i,lh,f);

                }
            }

#ifdef AMREX_USE_EB
            bool is_symmetry = (i_side_bc == PhysBCType::symmetry) || (i_side_bc == PhysBCType::slipwall) || (i_side_bc == PhysBCType::noslipwall);
            if (lh==0) {
                bs.eb_bc.setLo(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            } else {
                bs.eb_bc.setHi(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            }
#endif
        }
    }

    // check validity of inflow bc
    boundary_conditions.post_init();


    //
    // shock detector threshold
    //
    set_shock_detector();

    //
    // particles
    //

    particle_init = state_def["particles"].get_or<std::string>("");


    //
    // positivity
    //
    enforce_positivity = state_def["enforce_positivity"].get_or(0);

    extra_slope_limits = state_def["extra_slope_limits"].get_or(1);

    State::init_from_lua();
}



#ifdef AMREX_USE_EB
void HydroState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == HydroSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroSlipWall(flux_solver.get())));
    } else if (bc_type == HydroNoSlipWall::tag) {
        if (!viscous) {
            Abort("Requested EB bc of type '" + bc_type + "' without defining 'viscosity' for state '" + name + "'");
        }
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroNoSlipWall(flux_solver.get(), viscous.get(), bc_def)));
    } else if (bc_type == DirichletWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new DirichletWall(flux_solver.get(), get_prim_names(), get_prim_vector_idx(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}
#endif

Real HydroState::init_from_number_density(std::map<std::string, Real> data)
{
    Real nd = nd_function(data);
    Real mass_inv = 0.0;
    if (mass_const)
        return nd * mass[0];

    Real alphai = 0.0;
    Real S_alphas = 0.0;
    for (int i=0; i < n_alphas; ++i) {
        alphai = functions[+HydroState::PrimIdx::Density + 1 + i](data);
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        mass_inv += alphai / mass[i];
    }
        mass_inv += (1.0-S_alphas) / mass[n_alphas]; // deduce alpha value for final component

    return nd / mass_inv;
}

void HydroState::set_udf()
{
    using namespace std::placeholders;

    sol::state& lua = GD::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());


    if (!value.valid())
        Abort("State "+name+" does not have 'value' defined for initial conditions");

    // get a list of any initialisation functions that need to be called during the run
    sol::table dynamic = state_def["dynamic"].get_or(sol::table());

    const auto ignore = init_with_value();

    bool check, success;

    for (const auto& [prim_name, prim_idx] : HydroState::singular_prim_names_map) {

        // is there a variable with this name?
        success = false;
        check = value[prim_name].valid();

        // doesn't exist, is there an alternative?
        if (!check) {
            if (prim_name == "rho") {

                check = value["nd"].valid();

                if (!check)
                    Abort("State "+name+": 'rho' or 'nd' not defined for initial conditions");

                Optional3D1VFunction nd;
                success = get_udf(value["nd"], nd, 0.0);

                nd_function = nd;

                Optional3D1VFunction rho;

                rho.set_func(std::bind(&HydroState::init_from_number_density, this, _1));

                functions[prim_idx] = rho;
            }
        }

        // TODO: success does abolutely nothing?
        if (!success) {

            Optional3D1VFunction v;

            success = get_udf(value[prim_name], v, 0.0);

            if (not success) {
                for (const auto &j : ignore) {
                    if (prim_idx == j.first) {
                        v.set_value(j.second);
                        success = true;
                        break;
                    }
                }
            }

            functions[prim_idx] = v;
        }

    }

    // TODO: use a function to reduce redundant code
    if (is_multicomponent) {
        sol::table alpha_defs = value[multicomp_prim_name].get_or(sol::table());

        if (!alpha_defs.valid())
            Abort("State "+name+": 'alpha' list must be defined if n_components >= 2 to avoid ambiguity");

        if (alpha_defs.size() < n_alphas or alpha_defs.size() > n_alphas)
            Abort("State "+name+": n_components >= 2 so 'value["+multicomp_prim_name+"]' must be of size (n_components - 1)");

        for (int i = 0; i < n_alphas; ++i) {
            Optional3D1VFunction alpha = get_udf(alpha_defs[i+1]);

            if (!alpha.is_valid())
                Abort("State "+name+", 'value["+multicomp_prim_name+"]["+std::to_string(i)+"]' must be either a lua function, or numerical value");

            functions[+HydroState::PrimIdx::Density + 1 + i] = alpha;

        }
    }

    if (dynamic.valid()) {
        for (const auto &d : dynamic) {
            std::string dynamic_prim_name = d.second.as<std::string>();

            // name is not valid
            if (singular_prim_names_map.count(dynamic_prim_name) == 0 or dynamic_prim_name != multicomp_prim_name)
                Abort("State "+name+": 'dynamic', '"+dynamic_prim_name+"' is not a valid primitive name");

            if (dynamic_prim_name == multicomp_prim_name) {
                if (!is_multicomponent)
                    Abort("State "+name+": 'alpha' is listed as dynamic, but this is not a multicomponent state");

                for (int i = 0; i < n_alphas; ++i) {
                    dynamic_functions[+HydroState::PrimIdx::Density + 1 +i] = &functions[+HydroState::PrimIdx::Density + 1 +i];
                }

            } else {
                dynamic_functions[singular_prim_names_map[dynamic_prim_name]] = &functions[singular_prim_names_map[dynamic_prim_name]];
            }
        }
    }

}

AssociatedType HydroState::get_association_type() const
{
    // TODO: this needs to be able to work with a combined ions/neutral species definition
    if (all_less(charge.begin(), charge.end(), 0.0)) {
        return AssociatedType::Electron;
    } else if (all_greater(charge.begin(), charge.end(), 0.0)) {
        return AssociatedType::Ion;
    } else if (all_equal(charge.begin(), charge.end(), 0.0)) {
        return AssociatedType::Neutral;
    } else {
        return AssociatedType::Mix;
    }
}

int HydroState::n_prim() const {return +PrimIdx::NUM + n_components - 1;} // there are n-1 alpha values
int HydroState::n_cons() const {return +ConsIdx::NUM + n_components - 1;}
int HydroState::n_flux() const {return +FluxIdx::NUM + n_components - 1;}

const Vector<std::string>& HydroState::get_cons_names() const
{
    return cons_names;
}

const Vector<std::string>& HydroState::get_prim_names() const
{
    return prim_names;
}

const std::string HydroState::get_prim_name(int idx) const
{
    if (idx > +HydroState::PrimIdx::Density && is_multicomponent) {
        std::string s = "alpha"+std::to_string(idx-+HydroState::PrimIdx::NUM);
        return s;
    } else
        return get_prim_names()[idx];
}

const std::string HydroState::get_cons_name(int idx) const
{
    if (idx > +HydroState::ConsIdx::Density && is_multicomponent) {
        std::string s = "tracer"+std::to_string(idx-+HydroState::ConsIdx::NUM);
        return s;
    } else
        return get_cons_names()[idx];
}

const Vector<set_bc>& HydroState::get_bc_set() const
{
    return bc_set;
}

// TODO: clamp alpha values!

Vector<Real> HydroState::get_alphas_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_alphas_from_cons");

    Real rho = U[+HydroState::ConsIdx::Density];
    Vector<Real> tracers = Vector<Real>(U.begin() + +HydroState::ConsIdx::Density+1, U.end());
    Vector<Real> alphas;

    for (auto tr : tracers) {
        alphas.push_back(tr / rho);
    }

    return alphas;
}

Real HydroState::get_mass_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("HydroState::get_mass_from_prim");

    if (mass_const) return mass[0];

    Real S_alphai_mi = 0.0;
    Real S_alphas = 0.0;

    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = Q[+PrimIdx::Density + 1 + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        S_alphai_mi += alphai / mass[i];
    }

    S_alphai_mi += (1.0-S_alphas) / mass[n_alphas];

    return 1.0 / S_alphai_mi;
}

Real HydroState::get_mass_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_mass_from_cons");

    if (mass_const) return mass[0];

    Real rho = U[+ConsIdx::Density];

    Real S_alphai_mi = 0.0;
    Real S_alphas = 0.0;

    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = U[+ConsIdx::Density + 1 + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        S_alphai_mi += alphai / mass[i];
    }

    S_alphai_mi += (1.0-S_alphas) / mass[n_alphas];

    return 1.0 / S_alphai_mi;
}

Real HydroState::get_charge_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("HydroState::get_charge_from_prim");

    if (charge_const) return charge[0];

    Real S_alphaiqi_mi = 0;
    Real S_alphai_mi = 0;
    Real S_alphas = 0;

    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = Q[+PrimIdx::Density + 1 + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphaiqi_mi += alphai * charge[i] / mass[i];
        S_alphai_mi += alphai / mass[i];
        S_alphas += alphai;
    }

    S_alphaiqi_mi += (1.0-S_alphas) * charge[n_alphas] / mass[n_alphas];
    S_alphai_mi += (1.0-S_alphas) / mass[n_alphas];

    return S_alphaiqi_mi / S_alphai_mi;
}

Real HydroState::get_charge_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_charge_from_cons");

    if (charge_const) return charge[0];

    Real rho = U[+ConsIdx::Density];

    Real S_alphaiqi_mi = 0;
    Real S_alphai_mi = 0;
    Real S_alphas = 0;

    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = U[+ConsIdx::Density + 1 + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphaiqi_mi += alphai * charge[i] / mass[i];
        S_alphai_mi += alphai / mass[i];
        S_alphas += alphai;
    }

    S_alphaiqi_mi += (1.0-S_alphas) * charge[n_alphas] / mass[n_alphas];
    S_alphai_mi += (1.0-S_alphas) / mass[n_alphas];

    return S_alphaiqi_mi / S_alphai_mi;
}

Real HydroState::get_gamma_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("HydroState::get_gamma_from_prim");

    if (gamma_const) return gamma[0];

    Real S_alphaicpi = 0.0;
    Real S_alphaicvi = 0.0;
    Real S_alphai = 0.0;

    Real cvi = 0.0;
    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = Q[+PrimIdx::Density + 1 + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        cvi = 1.0 / (mass[i] * (gamma[i] - 1.0));

        S_alphaicpi += alphai * cpi;
        S_alphaicvi += alphai * cvi;
    }

    cpi = gamma[n_alphas] / (mass[n_alphas] * (gamma[n_alphas] - 1.0));
    cvi = 1.0 / (mass[n_alphas] * (gamma[n_alphas] - 1.0));

    S_alphaicpi += (1.0 - S_alphai) * cpi;
    S_alphaicvi += (1.0 - S_alphai) * cvi;

    return S_alphaicpi / S_alphaicvi;
}

Real HydroState::get_gamma_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_gamma_from_cons");

    if (gamma_const) return gamma[0];

    Real rho = U[+ConsIdx::Density];

    Real S_alphaicpi = 0.0;
    Real S_alphaicvi = 0.0;
    Real S_alphai = 0.0;

    Real cvi = 0.0;
    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = U[+ConsIdx::Density + 1 + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        cvi = 1.0 / (mass[i] * (gamma[i] - 1.0));

        S_alphaicpi += alphai * cpi;
        S_alphaicvi += alphai * cvi;
    }

    cpi = gamma[n_alphas] / (mass[n_alphas] * (gamma[n_alphas] - 1.0));
    cvi = 1.0 / (mass[n_alphas] * (gamma[n_alphas] - 1.0));

    S_alphaicpi += (1.0 - S_alphai) * cpi;
    S_alphaicvi += (1.0 - S_alphai) * cvi;

    return S_alphaicpi / S_alphaicvi;
}

Real HydroState::get_cp_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("HydroState::get_cp_from_prim");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real S_alphai = 0.0;
    Real S_alphaicpi = 0.0;

    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = Q[+PrimIdx::Density + 1 + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        S_alphaicpi += alphai * cpi;
    }

    cpi = gamma[n_alphas] / (mass[n_alphas] * (gamma[n_alphas] - 1.0));
    S_alphaicpi += (1.0- S_alphai) * cpi;

    return S_alphaicpi;
}

Real HydroState::get_cp_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_cp_from_cons");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real rho = U[+ConsIdx::Density];

    Real S_alphai = 0.0;
    Real S_alphaicpi = 0.0;

    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_alphas; ++i) {
        alphai = U[+ConsIdx::Density + 1 + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        S_alphaicpi += alphai * cpi;
    }

    cpi = gamma[n_alphas] / (mass[n_alphas] * (gamma[n_alphas] - 1.0));
    S_alphaicpi += (1.0- S_alphai) * cpi;

    return S_alphaicpi;

}


// in place conversion from conserved to primitive
bool HydroState::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("HydroState::cons2prim");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real m = get_mass_from_cons(U);
    Real g = get_gamma_from_cons(U);
    Real p = (ed - ke)*(g - 1);
    Real T = p*rhoinv*m;

    Q[+PrimIdx::Density] = rho;
    Q[+PrimIdx::Xvel] = u;
    Q[+PrimIdx::Yvel] = v;
    Q[+PrimIdx::Zvel] = w;
    Q[+PrimIdx::Prs] = p;
    Q[+PrimIdx::Temp] = T;

    for (int i = 0; i < n_alphas; ++i) {
        Q[+PrimIdx::Density + 1 + i] = U[+ConsIdx::Density + 1 + i] * rhoinv;
    }

    return prim_valid(Q);
}

// in place conversion from conserved to primitive
bool HydroState::cons2prim(Vector<autodiff::dual>& U, Vector<autodiff::dual>& Q) const
{
    BL_PROFILE("HydroState::cons2prim");

    autodiff::dual rho = U[+ConsIdx::Density];
    autodiff::dual mx = U[+ConsIdx::Xmom];
    autodiff::dual my = U[+ConsIdx::Ymom];
    autodiff::dual mz = U[+ConsIdx::Zmom];
    autodiff::dual ed = U[+ConsIdx::Eden];

    autodiff::dual rhoinv = 1/rho;
    autodiff::dual u = mx*rhoinv;
    autodiff::dual v = my*rhoinv;
    autodiff::dual w = mz*rhoinv;
    autodiff::dual ke = 0.5*rho*(u*u + v*v + w*w);

    Vector<Real> U_vals;
    for (int i = 0; i < U.size(); ++i) {
        U_vals.push_back(U[i].val);
    }

    Real m = get_mass_from_cons(U_vals);
    Real g = get_gamma_from_cons(U_vals);
    autodiff::dual p = (ed - ke)*(g - 1);
    autodiff::dual T = p*rhoinv*m;

    Q[+PrimIdx::Density] = rho;
    Q[+PrimIdx::Xvel] = u;
    Q[+PrimIdx::Yvel] = v;
    Q[+PrimIdx::Zvel] = w;
    Q[+PrimIdx::Prs] = p;
    Q[+PrimIdx::Temp] = T;

    for (int i = 0; i < n_alphas; ++i) {
        Q[+PrimIdx::Density + 1 + i] = U[+ConsIdx::Density + 1 + i] * rhoinv;
    }

    if ((Q[+PrimIdx::Density].val <= 0.0) ||  (Q[+PrimIdx::Prs].val <= 0.0)) {
        return false;
    }
    return true;
}

// in-place conversion from primitive to conserved variables
void HydroState::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("HydroState::prim2cons");

    Real rho = Q[+PrimIdx::Density];
    Real u = Q[+PrimIdx::Xvel];
    Real v = Q[+PrimIdx::Yvel];
    Real w = Q[+PrimIdx::Zvel];
    Real p = Q[+PrimIdx::Prs];

    Real mx = u*rho;
    Real my = v*rho;
    Real mz = w*rho;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real g = get_gamma_from_prim(Q);
    Real ed = p/(g - 1) + ke;

    U[+ConsIdx::Density] = rho;
    U[+ConsIdx::Xmom] = mx;
    U[+ConsIdx::Ymom] = my;
    U[+ConsIdx::Zmom] = mz;
    U[+ConsIdx::Eden] = ed;

    for (int i = 0; i < n_alphas; ++i) {
        U[+ConsIdx::Density + 1 + i] = Q[+PrimIdx::Density + 1 + i] * rho;
    }

}

bool HydroState::prim_valid(Vector<Real> &Q) const
{
    if ((Q[+PrimIdx::Density] <= 0.0) ||  (Q[+PrimIdx::Prs] <= 0.0)) {
        return false;
    }
    return true;
}

bool HydroState::cons_valid(Vector<Real> &U) const
{
    if ((U[+ConsIdx::Density] <= 0.0) ||  (U[+ConsIdx::Eden] <= 0.0)) {
        return false;
    }
    return true;
}

Real HydroState::get_energy_from_cons(const Vector<Real> &U) const
{
    return U[+ConsIdx::Eden];
}

Real HydroState::get_temperature_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_temperature_from_cons");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];

    Real rhoinv = 1/rho;
    Real ke = 0.5*rhoinv*(mx*mx + my*my + mz*mz);
    Real g = get_gamma_from_cons(U);

    Real prs = (ed - ke)*(g - 1);

    Real m = get_mass_from_cons(U);
    return prs*rhoinv*m;

}

Real HydroState::get_temperature_from_prim(const Vector<Real> &Q) const
{
    return Q[+PrimIdx::Temp];
}

dual HydroState::get_temperature_from_prim(const Vector<dual> &Q) const
{
    dual T = Q[+PrimIdx::Temp];
    return T;
}

RealArray HydroState::get_speed_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("HydroState::get_speed_from_cons");

    Real rho = U[+ConsIdx::Density];
    Real mx = U[+ConsIdx::Xmom];
    Real my = U[+ConsIdx::Ymom];
    Real mz = U[+ConsIdx::Zmom];
    Real ed = U[+ConsIdx::Eden];

    Real rhoinv = 1/rho;

    Real ux = mx*rhoinv;
    Real uy = my*rhoinv;
    Real uz = mz*rhoinv;

    Real kineng = 0.5*rho*(ux*ux + uy*uy + uz*uz);
    Real g = get_gamma_from_cons(U);
    Real prs = (ed - kineng)*(g - 1);
    Real a = std::sqrt(g*prs*rhoinv);

    RealArray s = {AMREX_D_DECL(a + std::abs(ux), a + std::abs(uy), a + std::abs(uz))};

    return s;

}

RealArray HydroState::get_speed_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("HydroState::get_speed_from_prim");

    Real g = get_gamma_from_prim(Q);

    Real a = std::sqrt(g*Q[+PrimIdx::Prs]/Q[+PrimIdx::Density]);

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+PrimIdx::Xvel]),
                                a + std::abs(Q[+PrimIdx::Yvel]),
                                a + std::abs(Q[+PrimIdx::Zvel]))};


    return s;

}

RealArray HydroState::get_current_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("HydroState::get_current_from_cons");
    Real q = get_charge_from_cons(U);
    Real m = get_mass_from_cons(U);
    Real r = q/m;

    RealArray c = {AMREX_D_DECL(
                   r*U[+ConsIdx::Xmom],
                   r*U[+ConsIdx::Ymom],
                   r*U[+ConsIdx::Zmom]
                   )};

    return c;
}



void HydroState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     EB_OPTIONAL(,const EBCellFlagFab &flag)
                                     EB_OPTIONAL(,const FArrayBox &vfrac)
                                     ) const
{
    BL_PROFILE("HydroState::calc_reconstruction");
    // if we don't want to apply extra limiting on the slopes (forced to 2nd order)
    // we can use the default reconstruction scheme
    if (!extra_slope_limits) {
        State::calc_reconstruction(box, prim, rlo, rhi EB_OPTIONAL(,flag,vfrac));
    }

    // convert pressure
    const Box &pbox = prim.box();
    const Dim3 p_lo = amrex::lbound(pbox);
    const Dim3 p_hi = amrex::ubound(pbox);

    FArrayBox gamma_minus_one(pbox);
    Array4<Real> const& src4 = prim.array();
    Array4<Real> const& gam4 = gamma_minus_one.array();

#ifdef AMREX_USE_EB
    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);

    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    int N = prim.nComp();

    Vector<Real> stencil(reconstruction->stencil_length);
    int offset = reconstruction->stencil_length/2;
    Array<int,3> stencil_index;
    Vector<Real> cell_value(N), cell_slope(N);

    Real rho_lo, rho_hi;
    Real alpha_lo, alpha_hi;
    Real abs_phi, phi_scale, coeff_eps;
    Real gam_lo, gam_hi;
    Vector<Real> prims(N);
    Vector<Real> prims_hi(N);
    Vector<Real> prims_lo(N);

    Vector<Real> phis_alpha;

    // make sure our arrays for putting lo and hi reconstructed values into
    // are the corect size
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rlo[d].resize(box, N);
        rhi[d].resize(box, N);

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif
    }

    // change pressure to internal energy
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif
                for (int h = 0; h < N; ++h) {
                    prims[h] = src4(i,j,k,h);
                }

                gam4(i,j,k) = get_gamma_from_prim(prims) - 1.0;

                src4(i,j,k,+PrimIdx::Prs) /= gam4(i,j,k);

            }
        }
    }

    // now do reconstruction

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (check_eb) {

                        // covered cell doesn't need calculating
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }

                        // cell that references a covered cell doesn't need calculating
                        bool skip = false;
                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            // check if any of the stencil values are from a covered cell
                            if (f4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2]).isCovered()) {
                                skip = true;
                                break;
                            }
                        }

                        if (skip) {
                            continue;
                        }

                    }
#endif

                    // cycle over all components
                    for (int n = 0; n<N; ++n) {

                        // fill in the stencil along dimension index
                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        cell_slope[n] = reconstruction->get_slope(stencil);
                        cell_value[n] = stencil[offset];

                    }


                    // apply corrections to slopes
                    // J. Sci. Comput. (2014) 60:584-611
                    // Robust Finite Volume Schemes for Two-Fluid Plasma Equations

                    Real &rho     = cell_value[+PrimIdx::Density];
                    Real &phi_rho = cell_slope[+PrimIdx::Density];

                    Real &u     = cell_value[+PrimIdx::Xvel];
                    Real &phi_u = cell_slope[+PrimIdx::Xvel];

                    Real &v     = cell_value[+PrimIdx::Yvel];
                    Real &phi_v = cell_slope[+PrimIdx::Yvel];

                    Real &w     = cell_value[+PrimIdx::Zvel];
                    Real &phi_w = cell_slope[+PrimIdx::Zvel];

                    Real &eps     = cell_value[+PrimIdx::Prs];
                    Real &phi_eps = cell_slope[+PrimIdx::Prs];

                    if (is_multicomponent) {
                        for (int a = 0 ; a < n_alphas; ++a) {
                            Real &alpha = cell_value[+PrimIdx::Density + 1 + a];
                            Real &phi_alpha = cell_slope[+PrimIdx::Density + 1 + a];

                            alpha_lo = alpha - 0.5*phi_alpha;
                            prims_lo[+PrimIdx::Density + 1 + a] = alpha_lo;

                            alpha_hi = alpha + 0.5*phi_alpha;
                            prims_hi[+PrimIdx::Density + 1 + a] = alpha_hi;

                            lo4(i,j,k,+PrimIdx::Density + 1 + a) = alpha_lo;
                            hi4(i,j,k,+PrimIdx::Density + 1 + a) = alpha_hi;

                        }
                    }

                    // correct density slope
                    if (std::abs(phi_rho) > 2*rho) {
                        phi_rho = 2*sign(phi_rho, 0.0)*rho;
                    }

                    // get some face values
                    rho_lo = rho - 0.5*phi_rho;
                    rho_hi = rho + 0.5*phi_rho;

                    gam_lo = get_gamma_from_prim(prims_lo);
                    gam_hi = get_gamma_from_prim(prims_hi);

                    abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;

                    // correct velocity slope
                    Real eps_face = eps - 0.5*std::abs(phi_eps);

                    if (eps_face <= 0.0) {
                        // if the reconstructed face value goes non-physical
                        // just set back to first order with zero slope
                        phi_u = 0.0;
                        phi_v = 0.0;
                        phi_w = 0.0;
                        phi_eps = 0.0;
                    } else {
                        coeff_eps = (rho/(rho_lo*rho_hi))*eps_face;
                        if ((0.125*abs_phi) > coeff_eps) {
                            phi_scale = sqrt(abs_phi);
                            coeff_eps = sqrt(8*coeff_eps);
                            phi_u = (phi_u/phi_scale)*coeff_eps;
                            phi_v = (phi_v/phi_scale)*coeff_eps;
                            phi_w = (phi_w/phi_scale)*coeff_eps;
                        }
                        // update eps
                        abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;
                        eps -= (rho_lo*rho_hi/rho)*0.125*abs_phi;
                    }



                    // density
                    lo4(i,j,k,+PrimIdx::Density) = rho_lo;
                    hi4(i,j,k,+PrimIdx::Density) = rho_hi;

                    // x - velocity
                    lo4(i,j,k,+PrimIdx::Xvel) = u - 0.5*(rho_hi/rho)*phi_u;
                    hi4(i,j,k,+PrimIdx::Xvel) = u + 0.5*(rho_lo/rho)*phi_u;

                    // y - velocity
                    lo4(i,j,k,+PrimIdx::Yvel) = v - 0.5*(rho_hi/rho)*phi_v;
                    hi4(i,j,k,+PrimIdx::Yvel) = v + 0.5*(rho_lo/rho)*phi_v;

                    // z - velocity
                    lo4(i,j,k,+PrimIdx::Zvel) = w - 0.5*(rho_hi/rho)*phi_w;
                    hi4(i,j,k,+PrimIdx::Zvel) = w + 0.5*(rho_lo/rho)*phi_w;

                    // epsilon -> pressure
                    lo4(i,j,k,+PrimIdx::Prs) = (eps - 0.5*phi_eps)*(gam_lo - 1.0);
                    hi4(i,j,k,+PrimIdx::Prs) = (eps + 0.5*phi_eps)*(gam_hi - 1.0);


                    // Temperature (calculate from pressure and density)
                    lo4(i,j,k,+PrimIdx::Temp) = lo4(i,j,k,+PrimIdx::Prs)/(rho_lo/get_mass_from_prim(prims_lo));
                    hi4(i,j,k,+PrimIdx::Temp) = hi4(i,j,k,+PrimIdx::Prs)/(rho_hi/get_mass_from_prim(prims_hi));

                }
            }
        }
    }


    // convert back to pressure
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {
#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif
                src4(i,j,k,+PrimIdx::Prs) *= gam4(i,j,k);

            }
        }
    }

    return;
}

void HydroState::get_state_values(const Box& box,
                                  const FArrayBox& src,
                                  std::map<std::string,FArrayBox>& out,
                                  Vector<std::string>& updated
                                  EB_OPTIONAL(,const FArrayBox& vfrac)
                                  ) const
{
    BL_PROFILE("HydroState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<n_cons(); ++i) {
        const std::string s = get_cons_name(i);
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<n_prim(); ++i) {
        const std::string s = get_prim_name(i);
        if (s == get_cons_name(0)) continue;
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // additional variables

    Vector<std::string> other;

    const std::string charge_name = "charge-"+name;
    bool load_charge = out.find(charge_name) != out.end();
    if (load_charge) other.push_back(charge_name);

    const std::string mass_name = "mass-"+name;
    bool load_mass = out.find(mass_name) != out.end();
    if (load_mass) other.push_back(mass_name);

    const std::string gamma_name = "gamma-"+name;
    bool load_gamma = out.find(gamma_name) != out.end();
    if (load_gamma) other.push_back(gamma_name);

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
    if (load_vfrac) other.push_back(vfrac_name);
#endif

    updated.insert(updated.end(), other.begin(), other.end());

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out[s].setVal(0.0);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Vector<Real> S(n_cons()), Q(n_prim());

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    continue;
                }
#endif

                get_state_vector(src, i, j, k, S);

                if (load_charge) out4[charge_name](i,j,k) = get_charge_from_cons(S);
                if (load_mass)   out4[mass_name](i,j,k)   = get_mass_from_cons(S);
                if (load_gamma)  out4[gamma_name](i,j,k)  = get_gamma_from_cons(S);
            #ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
            #endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

                if (!prim_tags.empty()) {
                    cons2prim(S, Q);

                    for (const auto& var : prim_tags) {
                        out4[var.first](i,j,k) = Q[var.second];
                    }
                }
            }
        }
    }


    return;
}

void HydroState::calc_velocity(const Box& box,
                               FArrayBox& src,
                               FArrayBox& vel
                               EB_OPTIONAL(,const EBCellFlagFab& flag)
                               ) const
{
    BL_PROFILE("HydroState::calc_velocity");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real> const& src4 = src.array();
    Array4<Real> const& vel4 = vel.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    int nc = vel.nComp();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    for (int n=0; n<nc; ++n) {
                        vel4(i,j,k,n) = 0.0;
                    }
                    continue;
                }
#endif

                Real invrho = 1.0/src4(i,j,k,+ConsIdx::Density);

                for (int n=0; n<nc; ++n) {
                    vel4(i,j,k,n) = src4(i,j,k,+ConsIdx::Xmom+n)*invrho;
                }
            }
        }
    }

    return;
}

// given all of the available face values load the ones expected by the flux calc into a vector
void HydroState::load_state_for_flux(Vector<Array4<const Real>> &face,
                                               int i, int j, int k, Vector<Real> &S) const
{
    BL_PROFILE("HydroState::load_state_for_flux");


    // first get the primitives of this state
    Array4<const Real> const &f4 = face[global_idx];
    Vector<Real> prims(n_prim());

    S[+FluxIdx::Density] = f4(i,j,k,+PrimIdx::Density);
    S[+FluxIdx::Xvel] = f4(i,j,k,+PrimIdx::Xvel);
    S[+FluxIdx::Yvel] = f4(i,j,k,+PrimIdx::Yvel);
    S[+FluxIdx::Zvel] = f4(i,j,k,+PrimIdx::Zvel);
    S[+FluxIdx::Prs] = f4(i,j,k,+PrimIdx::Prs);

    if (is_multicomponent) {
        for (int a = 0; a < n_alphas; ++a) {
            prims[+PrimIdx::Density + 1 + a] = f4(i,j,k,+PrimIdx::Density + 1 + a);
            S[+FluxIdx::Density + 1 + a] = f4(i,j,k,+PrimIdx::Density + 1 + a);
        }
    }
    S[+FluxIdx::Gamma] = get_gamma_from_prim(prims);

    return;
}

void HydroState::calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox> &prim,
                                     #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                     #endif
                                     const Real* dx) const
{
    BL_PROFILE("HydroState::calc_viscous_fluxes");
    // now calculate viscous fluxes and load them into the flux arrays
    if (viscous) {

        Viscous &V = *viscous;
        switch (V.get_type()) {
            case Viscous::Neutral :

//                plot_FAB_2d(flag, "flag", false);
//                plot_FAB_2d(prim[global_idx], 0, "prim density", false, true);

                calc_neutral_viscous_fluxes(box,
                                            fluxes,
                                            pbox,
                                            prim[global_idx],
                                            EB_OPTIONAL(flag,)
                                            dx);
                break;
            case Viscous::Ion :
                calc_ion_viscous_fluxes(box,
                                        fluxes,
                                        pbox,
                                        prim,
                                        EB_OPTIONAL(flag,)
                                        dx);
                break;
            case Viscous::Electron :
                calc_electron_viscous_fluxes(box,
                                             fluxes,
                                             pbox,
                                             prim,
                                             EB_OPTIONAL(flag,)
                                             dx);
                break;
            default :
                break;
        }
    }

}

void HydroState::calc_neutral_diffusion_terms(const Box& box,
                                              const FArrayBox& prim,
                                              FArrayBox& diff
                                              EB_OPTIONAL(,const EBCellFlagFab& flag)
                                              ) const
{
    BL_PROFILE("HydroState::calc_neutral_diffusion_terms");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& prim4 = prim.array();
    Array4<Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Viscous &V = *viscous;

    int np = n_prim();

    Vector<Real> Q(np);
    Real T, mu, kappa;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                for (int n=0; n<np; ++n) {
                    Q[n] = prim4(i,j,k,n);
                }

                V.get_neutral_coeffs(Q, T, mu, kappa);


                d4(i,j,k,Viscous::NeutralTemp) = T;
                d4(i,j,k,Viscous::NeutralKappa) = kappa;
                d4(i,j,k,Viscous::NeutralMu) = mu;
            }
        }
    }

    return;
}

#ifdef AMREX_USE_EB
void HydroState::calc_neutral_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                                AMREX_SPACEDIM> &fluxes,
                                                const Box& pbox,
                                                const FArrayBox &prim,
                                                EB_OPTIONAL(const EBCellFlagFab& flag,)
                                                const Real* dx) const
{
    BL_PROFILE("HydroState::calc_neutral_viscous_fluxes_eb");
    FArrayBox diff(pbox, Viscous::NUM_NEUTRAL_DIFF_COEFFS);
    calc_neutral_diffusion_terms(pbox,
                                 prim,
                                 diff
                                 EB_OPTIONAL(,flag)
                                 );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

    Array4<const EBCellFlag> const& f4 = flag.array();

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    const Array<Real,3>  weights = {0.0, 1.0, 0.5};
    Real whi, wlo;

    Real tauxx, tauxy, tauxz, dTdx, muf;

    int iTemp = Viscous::NeutralTemp;
    int iMu = Viscous::NeutralMu;
    int iKappa = Viscous::NeutralKappa;

    Vector<int> prim_vel_id = get_prim_vector_idx();

    int Xvel = prim_vel_id[0] + 0;
    int Yvel = prim_vel_id[0] + 1;
    int Zvel = prim_vel_id[0] + 2;

    Vector<int> cons_vel_id = get_cons_vector_idx();

    int Xmom = cons_vel_id[0] + 0;
    int Ymom = cons_vel_id[0] + 1;
    int Zmom = cons_vel_id[0] + 2;

    Vector<int> nrg_id = get_nrg_idx();
    int Eden = nrg_id[0];


    // X - direction
    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y-AMREX_D_PICK(0,1,1); j <= hi.y+AMREX_D_PICK(0,1,1); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(-1,0,0);
                bool other_covered = f4(i-1,j,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdx = (d4(i,j,k,iTemp) - d4(i-1,j,k,iTemp))*dxinv[0];

                dudx = (p4(i,j,k,Xvel) - p4(i-1,j,k,Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,Yvel) - p4(i-1,j,k,Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,Zvel) - p4(i-1,j,k,Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2

                const int jhip = j + (int)f4(i,j,k).isConnected(0, 1,0);
                const int jhim = j - (int)f4(i,j,k).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i-1,j,k).isConnected(0, 1,0);
                const int jlom = j - (int)f4(i-1,j,k).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];
                dudy = (0.5*dxinv[1]) * ((p4(i  ,jhip,k,Xvel)-p4(i  ,jhim,k,Xvel))*whi+(p4(i-1,jlop,k,Xvel)-p4(i-1,jlom,k,Xvel))*wlo);
                dvdy = (0.50*dxinv[1]) * ((p4(i  ,jhip,k,Yvel)-p4(i  ,jhim,k,Yvel))*whi+(p4(i-1,jlop,k,Yvel)-p4(i-1,jlom,k,Yvel))*wlo);

#endif
#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i-1,j,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i-1,j,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];
                dudz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,Xvel)-p4(i  ,j,khim,Xvel))*whi + (p4(i-1,j,klop,Xvel)-p4(i-1,j,klom,Xvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,Zvel)-p4(i  ,j,khim,Zvel))*whi + (p4(i-1,j,klop,Zvel)-p4(i-1,j,klom,Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,iMu)+d4(i-1,j,k,iMu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,Xmom) -= tauxx;
                fluxX(i,j,k,Ymom) -= tauxy;
                fluxX(i,j,k,Zmom) -= tauxz;
                fluxX(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel) +  p4(i-1,j,k,Xvel))*tauxx
                                          +(p4(i,j,k,Yvel) + p4(i-1,j,k,Yvel))*tauxy
                                          +(p4(i,j,k,Zvel) + p4(i-1,j,k,Zvel))*tauxz
                                          +(d4(i,j,k,iKappa)+d4(i-1,j,k,iKappa))*dTdx);
            }
        }
    }

    // Y - direction
#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,-1,0);
                bool other_covered = f4(i,j-1,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdy = (d4(i,j,k,iTemp)-d4(i,j-1,k,iTemp))*dxinv[1];
                dudy = (p4(i,j,k,Xvel)-p4(i,j-1,k,Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,Yvel)-p4(i,j-1,k,Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,Zvel)-p4(i,j-1,k,Zvel))*dxinv[1];

                const int ihip = i + (int)f4(i,j,  k).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,  k).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j-1,k).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j-1,k).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,Xvel)-p4(ihim,j  ,k,Xvel))*whi + (p4(ilop,j-1,k,Xvel)-p4(ilom,j-1,k,Xvel))*wlo);
                dvdx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,Yvel)-p4(ihim,j  ,k,Yvel))*whi + (p4(ilop,j-1,k,Yvel)-p4(ilom,j-1,k,Yvel))*wlo);

#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,  k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,  k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i,j-1,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i,j-1,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];

                dvdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,Yvel)-p4(i,j  ,khim,Yvel))*whi + (p4(i,j-1,klop,Yvel)-p4(i,j-1,klom,Yvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,Zvel)-p4(i,j  ,khim,Zvel))*whi + (p4(i,j-1,klop,Zvel)-p4(i,j-1,klom,Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j-1,k,iMu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,Xmom) -= tauxy;
                fluxY(i,j,k,Ymom) -= tauyy;
                fluxY(i,j,k,Zmom) -= tauyz;
                fluxY(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j-1,k,Xvel))*tauxy
                                          +(p4(i,j,k,Yvel)+p4(i,j-1,k,Yvel))*tauyy
                                          +(p4(i,j,k,Zvel)+p4(i,j-1,k,Zvel))*tauyz
                                          +(d4(i,j,k,iKappa) + d4(i,j-1,k,iKappa))*dTdy);

            }
        }
    }
#endif

    // Z - direction
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z+1; ++k) {
        for   (int j = lo.y-1; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,0,-1);
                bool other_covered = f4(i,j,k-1).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];

                const int ihip = i + (int)f4(i,j,k  ).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,k  ).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j,k-1).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j,k-1).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,Xvel)-p4(ihim,j,k  ,Xvel))*whi + (p4(ilop,j,k-1,Xvel)-p4(ilom,j,k-1,Xvel))*wlo);
                dwdx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,Zvel)-p4(ihim,j,k  ,Zvel))*whi + (p4(ilop,j,k-1,Zvel)-p4(ilom,j,k-1,Zvel))*wlo);

                const int jhip = j + (int)f4(i,j,k  ).isConnected(0 ,1,0);
                const int jhim = j - (int)f4(i,j,k  ).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i,j,k-1).isConnected(0 ,1,0);
                const int jlom = j - (int)f4(i,j,k-1).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];

                dvdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,Yvel)-p4(i,jhim,k  ,Yvel))*whi + (p4(i,jlop,k-1,Yvel)-p4(i,jlom,k-1,Yvel))*wlo);
                dwdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,Zvel)-p4(i,jhim,k  ,Zvel))*whi + (p4(i,jlop,k-1,Zvel)-p4(i,jlom,k-1,Zvel))*wlo);

                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j,k-1,iMu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,Xmom) -= tauxz;
                fluxZ(i,j,k,Ymom) -= tauyz;
                fluxZ(i,j,k,Zmom) -= tauzz;
                fluxZ(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j,k-1,Xvel))*tauxz
                                          +(p4(i,j,k,Yvel)+p4(i,j,k-1,Yvel))*tauyz
                                          +(p4(i,j,k,Zvel)+p4(i,j,k-1,Zvel))*tauzz
                                          +(d4(i,j,k,iKappa) +d4(i,j,k-1,iKappa))*dTdz);

            }
        }
    }

#endif

}

#endif

void HydroState::calc_neutral_viscous_fluxes(const Box& box, Array<FArrayBox,
                                             AMREX_SPACEDIM> &fluxes,
                                             const Box& pbox,
                                             const FArrayBox &prim,
                                             EB_OPTIONAL(const EBCellFlagFab& flag,)
                                             const Real* dx) const
{
    BL_PROFILE("HydroState::calc_neutral_viscous_fluxes");
#ifdef AMREX_USE_EB
    if (flag.getType() != FabType::regular) {
        calc_neutral_viscous_fluxes_eb(box, fluxes, pbox, prim, flag, dx);
        return;
    }
#endif

    FArrayBox diff(pbox, Viscous::NUM_NEUTRAL_DIFF_COEFFS);
    calc_neutral_diffusion_terms(pbox,
                                 prim,
                                 diff
                             #ifdef AMREX_USE_EB
                                 ,flag
                             #endif
                                 );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    Real tauxx, tauxy, tauxz, dTdx, muf;

    int iTemp = Viscous::NeutralTemp;
    int iMu = Viscous::NeutralMu;
    int iKappa = Viscous::NeutralKappa;

    Vector<int> prim_vel_id = get_prim_vector_idx();

    int Xvel = prim_vel_id[0] + 0;
    int Yvel = prim_vel_id[0] + 1;
    int Zvel = prim_vel_id[0] + 2;

    Vector<int> cons_vel_id = get_cons_vector_idx();

    int Xmom = cons_vel_id[0] + 0;
    int Ymom = cons_vel_id[0] + 1;
    int Zmom = cons_vel_id[0] + 2;

    Vector<int> nrg_id = get_nrg_idx();
    int Eden = nrg_id[0];

    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdx = (d4(i,j,k,iTemp) - d4(i-1,j,k,iTemp))*dxinv[0];

                dudx = (p4(i,j,k,Xvel) - p4(i-1,j,k,Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,Yvel) - p4(i-1,j,k,Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,Zvel) - p4(i-1,j,k,Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2
                dudy = (p4(i,j+1,k,Xvel)+p4(i-1,j+1,k,Xvel)-p4(i,j-1,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[1]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i-1,j+1,k,Yvel)-p4(i,j-1,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[1]);
#endif
#if AMREX_SPACEDIM == 3
                dudz = (p4(i,j,k+1,Xvel)+p4(i-1,j,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i-1,j,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,iMu)+d4(i-1,j,k,iMu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,Xmom) -= tauxx;
                fluxX(i,j,k,Ymom) -= tauxy;
                fluxX(i,j,k,Zmom) -= tauxz;
                fluxX(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel) +  p4(i-1,j,k,Xvel))*tauxx
                                          +(p4(i,j,k,Yvel) + p4(i-1,j,k,Yvel))*tauxy
                                          +(p4(i,j,k,Zvel) + p4(i-1,j,k,Zvel))*tauxz
                                          +(d4(i,j,k,iKappa)+d4(i-1,j,k,iKappa))*dTdx);

            }
        }
    }

#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdy = (d4(i,j,k,iTemp)-d4(i,j-1,k,iTemp))*dxinv[1];
                dudy = (p4(i,j,k,Xvel)-p4(i,j-1,k,Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,Yvel)-p4(i,j-1,k,Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,Zvel)-p4(i,j-1,k,Zvel))*dxinv[1];
                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j-1,k,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[0]);
                dvdx = (p4(i+1,j,k,Yvel)+p4(i+1,j-1,k,Yvel)-p4(i-1,j,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[0]);
#if AMREX_SPACEDIM == 3
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i,j-1,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j-1,k,iMu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,Xmom) -= tauxy;
                fluxY(i,j,k,Ymom) -= tauyy;
                fluxY(i,j,k,Zmom) -= tauyz;
                fluxY(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j-1,k,Xvel))*tauxy
                                          +(p4(i,j,k,Yvel)+p4(i,j-1,k,Yvel))*tauyy
                                          +(p4(i,j,k,Zvel)+p4(i,j-1,k,Zvel))*tauyz
                                          +(d4(i,j,k,iKappa) + d4(i,j-1,k,iKappa))*dTdy);

            }
        }
    }


#endif
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];
                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j,k-1,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[0]);
                dwdx = (p4(i+1,j,k,Zvel)+p4(i+1,j,k-1,Zvel)-p4(i-1,j,k,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[0]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i,j+1,k-1,Yvel)-p4(i,j-1,k,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[1]);
                dwdy = (p4(i,j+1,k,Zvel)+p4(i,j+1,k-1,Zvel)-p4(i,j-1,k,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[1]);
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j,k-1,iMu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,Xmom) -= tauxz;
                fluxZ(i,j,k,Ymom) -= tauyz;
                fluxZ(i,j,k,Zmom) -= tauzz;
                fluxZ(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j,k-1,Xvel))*tauxz
                                          +(p4(i,j,k,Yvel)+p4(i,j,k-1,Yvel))*tauyz
                                          +(p4(i,j,k,Zvel)+p4(i,j,k-1,Zvel))*tauzz
                                          +(d4(i,j,k,iKappa) +d4(i,j,k-1,iKappa))*dTdz);

            }
        }
    }

#endif


    return;
}

// ====================================================================================
void HydroState::calc_ion_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                          State& EMstate,Array4<const Real> const& prim_EM4,
                                          State& ELEstate,Array4<const Real> const& prim_ELE4,
                                          FArrayBox& diff
                                          EB_OPTIONAL(,const EBCellFlagFab& flag)
                                          ) const {

    BL_PROFILE("HydroState::calc_ion_diffusion_terms");
    return;
}

void HydroState::calc_ion_viscous_fluxes(const Box& box, 
                                         Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                         const Box& pbox, const Vector<FArrayBox>& prim,
                                         EB_OPTIONAL(const EBCellFlagFab& flag,)
                                         const Real* dx) const {

    BL_PROFILE("HydroState::calc_ion_viscous_fluxes");
    return;
}

// ====================================================================================

void HydroState::calc_electron_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                               State& EMstate,
                                               Array4<const Real> const& prim_EM4,
                                               State& IONstate,
                                               Array4<const Real> const& prim_ION4,
                                               FArrayBox& diff
                                               EB_OPTIONAL(,const EBCellFlagFab& flag)
                                               ) const {
    BL_PROFILE("HydroState::calc_electron_diffusion_terms");
    return;
}

void HydroState::calc_electron_viscous_fluxes(const Box& box, 
                                              Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                              const Box& pbox, const Vector<FArrayBox>& prim,
                                              EB_OPTIONAL(const EBCellFlagFab& flag,)
                                              const Real* dx) const {
    BL_PROFILE("HydroState::calc_electron_viscous_fluxes");
    return;
}

// ====================================================================================

void HydroState::calc_charged_viscous_fluxes(int passed_idx,
                                             int ion_idx,
                                             int electron_idx,
                                             int em_idx,
                                             const Box& box,
                                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                             const Box& pbox,
                                             const Vector<FArrayBox>& prim,
                                             EB_OPTIONAL(const EBCellFlagFab& flag,)
                                             const Real* dx, FArrayBox& diff) const {
    BL_PROFILE("HydroState::calc_charged_viscous_fluxes");
    return;
}


// ====================================================================================
void HydroState::write_info(nlohmann::json& js) const
{

    State::write_info(js);

    js["mass"] = mass;
    js["charge"] = charge;
    js["gamma"] = gamma;

    if (viscous) {

        auto& grp = js["viscosity"];

        int tp = viscous->get_type();
        grp["type"] = tp;

        const auto coeffs = viscous->get_refs();

        for (const auto& cf : coeffs) {
            grp[cf.first] = cf.second;
        }
    }
}

std::string HydroState::str() const
{
    std::stringstream msg;

    msg << State::str();

    if (viscous) {

        msg << "    viscosity : " << viscous->str() << "\n";

    }

    return msg.str();
}
