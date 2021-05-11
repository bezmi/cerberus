#ifndef MFP_HYDRO_H
#define MFP_HYDRO_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include "MFP_bc.H"
#include "MFP_utility.H"
#include "MFP_state.H"

using namespace amrex;

struct HydroState : public State
{
public:
    HydroState();
    HydroState(const sol::table& def);
    ~HydroState();

    enum class ConsIdx : int {
        Xmom=0,
        Ymom,
        Zmom,
        Eden,
        Density,
        //Tracer,
        NUM,
    };

    enum class PrimIdx : int {
        Xvel=0,
        Yvel,
        Zvel,
        Prs,
        Temp,
        Density,
        //Alpha,
        NUM
    };

    enum class FluxIdx : int {
        Xvel=0,
        Yvel,
        Zvel,
        Prs,
        Gamma,
        Density,
        //Alpha,
        NUM
    };

    const std::string multicomp_prim_name = "alpha";
    const std::string multicomp_cons_name = "tracer";


// set cons/prim_names based on ordering of Cond/PrimIdx
    std::map<std::string, int> singular_cons_names_map = {{"rho", +HydroState::ConsIdx::Density},
                                                            {"x_mom", +HydroState::ConsIdx::Xmom},
                                                            {"y_mom", +HydroState::ConsIdx::Ymom},
                                                            {"z_mom", +HydroState::ConsIdx::Zmom},
                                                            {"nrg", +HydroState::ConsIdx::Eden}};

    std::map<std::string, int> singular_prim_names_map = {{"rho", +HydroState::PrimIdx::Density},
                                                            {"x_vel", +HydroState::PrimIdx::Xvel},
                                                            {"y_vel", +HydroState::PrimIdx::Yvel},
                                                            {"z_vel", +HydroState::PrimIdx::Zvel},
                                                            {"p", +HydroState::PrimIdx::Prs},
                                                            {"T", +HydroState::PrimIdx::Temp}};
    Optional3D1VFunction nd_function;


    static Vector<int> flux_vector_idx;
    static Vector<int> cons_vector_idx;
    static Vector<int> prim_vector_idx;

    const std::map<std::string, int> bc_names = {{"interior",  PhysBCType::interior},
                                                 {"inflow",    PhysBCType::inflow},
                                                 {"outflow",   PhysBCType::outflow},
                                                 {"symmetry",  PhysBCType::symmetry},
                                                 {"slipwall",  PhysBCType::slipwall},
                                                 {"noslipwall",PhysBCType::noslipwall}};

    static bool registered;
    static std::string tag;

    bool is_multicomponent = false;
    bool is_two_component = false;

    static Vector<set_bc> bc_set;

    bool enforce_positivity = false;
    // const Vector<int> enforce_positivity_idx = {+ConsIdx::Density, +ConsIdx::Eden};

    bool extra_slope_limits = true;

    int n_components = 0;
    int n_alphas = 0;

    virtual void init_from_lua() override;
#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

    Real init_from_number_density(std::map<std::string, Real> data);
    virtual void set_udf() override;

    virtual const std::string& get_tag() const override {return tag;}
    virtual int get_type() const override {return +StateType::isHydro;}
    virtual AssociatedType get_association_type() const override;
    virtual int n_prim() const override;
    virtual int n_cons() const override;
    virtual int n_flux() const override;

    virtual bool get_enforce_positivity() const override {return enforce_positivity;}
    // virtual Vector<int> get_enforce_positivity_idx() const override {return enforce_positivity_idx;}

    virtual const Vector<std::string>& get_cons_names() const override;
    virtual const Vector<std::string>& get_prim_names() const override;
    virtual const std::string get_prim_name(int idx) const override;
    virtual const std::string get_cons_name(int idx) const override;
    virtual const Vector<set_bc>& get_bc_set() const override;

    Vector<Real> get_alphas_from_cons(const Vector<Real> &U) const;

    virtual Real get_mass_from_cons(const Vector<Real> &U) const override;
    virtual Real get_charge_from_cons(const Vector<Real> &U) const override;
    Real get_gamma_from_cons(const Vector<Real> &U) const;
    Real get_cp_from_cons(const Vector<Real> &U) const;

    Real get_mass_from_prim(const Vector<Real> &Q) const;
    Real get_charge_from_prim(const Vector<Real> &Q) const;
    Real get_gamma_from_prim(const Vector<Real> &Q) const;
    Real get_cp_from_prim(const Vector<Real> &Q) const;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual bool cons2prim(Vector<autodiff::dual>& U, Vector<autodiff::dual>& Q) const override;

    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;
    virtual bool prim_valid(Vector<Real>& Q) const override;
    virtual bool cons_valid(Vector<Real>& U) const override;

    virtual int get_cons_density_idx() const override {
        return +ConsIdx::Density;
    }
    virtual const Vector<int> get_prim_vector_idx() const override {
        return prim_vector_idx;
    }

    virtual const Vector<int> get_flux_vector_idx() const override {
        return flux_vector_idx;
    }

    virtual const Vector<int> get_cons_vector_idx() const override {
        return cons_vector_idx;
    }

    virtual const Vector<int> get_prs_idx() const override {
        return {+PrimIdx::Prs};
    }
    virtual const Vector<int> get_nrg_idx() const override {
        return {+ConsIdx::Eden};
    }

    virtual Real get_density_from_cons(const Vector<Real> &U) const override
    {
        return U[+ConsIdx::Density];
    }
    virtual Real get_density_from_prim(const Vector<Real> &Q) const override
    {
        return Q[+PrimIdx::Density];
    }
    virtual dual get_density_from_prim(const Vector<dual> &Q) const override
    {
        return Q[+PrimIdx::Density];
    }

    virtual Real get_energy_from_cons(const Vector<Real> &U) const override;

    virtual Real get_temperature_from_cons(const Vector<Real> &U) const override;
    virtual Real get_temperature_from_prim(const Vector<Real> &Q) const override;
    virtual dual get_temperature_from_prim(const Vector<dual> &Q) const override;


    virtual RealArray get_speed_from_cons(const Vector<Real> &U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real> &Q) const override;

    virtual RealArray get_current_from_cons(const Vector<Real> &U) const override;

    void calc_reconstruction(const Box& box,
                             FArrayBox &prim,
                             Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                             Array<FArrayBox, AMREX_SPACEDIM> &rhi
                             EB_OPTIONAL(,const EBCellFlagFab& flag)
                             EB_OPTIONAL(,const FArrayBox &vfrac)
                             ) const override;

    void get_state_values(const Box& box,
                                  const FArrayBox& src,
                                  std::map<std::string,FArrayBox>& out,
                                  Vector<std::string>& updated
                                  EB_OPTIONAL(,const FArrayBox& vfrac)
                                  ) const override;


    void calc_velocity(const Box& box,
                       FArrayBox& src,
                       FArrayBox& vel
                       EB_OPTIONAL(,const EBCellFlagFab& flag)
                       ) const override;

    virtual void load_state_for_flux(Vector<Array4<const Real> > &face,
                                             int i, int j, int k, Vector<Real> &S) const override;

    virtual void calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox> &prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx) const override;

    void calc_neutral_diffusion_terms(const Box& box,
                                      const FArrayBox &prim,
                                      FArrayBox& diff
                                      EB_OPTIONAL(,const EBCellFlagFab& flag)
                                      ) const;

    void calc_neutral_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box &pbox,
                                     const FArrayBox &prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx) const;

    void calc_charged_viscous_fluxes(int passed_idx,
                                     int ion_idx,
                                     int electron_idx,
                                     int em_idx,
                                     const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox>& prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx,
                                     FArrayBox& diff) const;

    void calc_ion_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                  State& EMstate,Array4<const Real> const& prim_EM4,
                                  State& ELEstate,Array4<const Real> const& prim_ELE4,
                                  FArrayBox& diff
                                  EB_OPTIONAL(,const EBCellFlagFab& flag)
                                  ) const ;

    void calc_ion_viscous_fluxes(const Box& box,
                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                 const Box &pbox,
                                 const Vector<FArrayBox>& prim,
                                 EB_OPTIONAL(const EBCellFlagFab& flag,)
                                 const Real* dx) const;

    void calc_electron_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                       State& EMstate,Array4<const Real> const& prim_EM4,
                                       State& ELEstate,Array4<const Real> const& prim_ELE4,
                                       FArrayBox& diff
                                       EB_OPTIONAL(,const EBCellFlagFab& flag)
                                       ) const ;
    
    void calc_electron_viscous_fluxes(const Box& box, Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                      const Box &pbox, const Vector<FArrayBox>& prim,
                                      EB_OPTIONAL(const EBCellFlagFab& flag,)
                                      const Real* dx) const;

#ifdef AMREX_USE_EB
    void calc_neutral_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                        AMREX_SPACEDIM> &fluxes,
                                        const Box& pbox,
                                        const FArrayBox &prim,
                                        const EBCellFlagFab& flag,
                                        const Real* dx) const;
#endif

    virtual void write_info(nlohmann::json &js) const override;

    virtual std::string str() const override;
};

#endif // MFP_HYDRO_H
