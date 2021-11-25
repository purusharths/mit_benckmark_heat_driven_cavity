// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

/// \author Philipp Gerstner

// System includes.
#include "hiflow.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType DataType;
typedef LAD::VectorType VectorType;
typedef LAD::MatrixType MatrixType;

// DIM of the problem.
const int DIM = 2;

typedef Vec<DIM, DataType> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Functor used to impose u(x) = c on the boundary.
struct VelocityDirichletBC
{
  VelocityDirichletBC(DataType time,
                      int top_mat_number,
                      int bottom_mat_number)
      : time_(time),
        top_mat_number_(top_mat_number),
        bottom_mat_number_(bottom_mat_number) {}

  void evaluate(const mesh::Entity &face,
                const Vec<DIM, DataType> &pt_coord,
                std::vector<DataType> &vals) const
  {
    const int material_number = face.get_material_number();
    vals.clear();

    // **********************************************
    // TODO exercise A
    // DataType time_factor = 1.;
    if (material_number == top_mat_number_ || material_number == bottom_mat_number_)
    {
      vals.resize(DIM, 0.);
    }
    // END Exercise A
    // *********************************************
  }

  size_t nb_comp() const
  {
    return DIM;
  }

  size_t nb_func() const
  {
    return 1;
  }

  size_t iv2ind(size_t j, size_t v) const
  {
    return v;
  }

  int top_mat_number_;
  int bottom_mat_number_;
  DataType inflow_vel_x_;
  DataType time_;
};
// make changes
// @@@ NOTE: should return/modify temp not velocity
struct TemperatureDirichletBC
{
  TemperatureDirichletBC(DataType time,
                         int lhs_mat_number,
                         int rhs_mat_number)
      : time_(time),
        lhs_mat_number_(lhs_mat_number),
        rhs_mat_number_(rhs_mat_number)
  {
  }

  void evaluate(const mesh::Entity &face,
                const Vec<DIM, DataType> &pt_coord,
                std::vector<DataType> &vals) const
  {
    const int material_number = face.get_material_number();
    vals.clear();

    // **********************************************
    // TODO exercise A
    // DataType time_factor = 1.;
    if (material_number == lhs_mat_number_)
    {
      // vals.resize(1, 0.);
      vals.resize(DIM, 0.);
      vals[0] = 0.5; //  .5 to temp
    }

    if (material_number == rhs_mat_number_)
    {
      // vals.resize(1, 0.);
      vals.resize(DIM, 0.);
      vals[0] = -0.5; // -.5 to temp
    }
    // END Exercise A
    // *********************************************
  }

  size_t nb_comp() const
  {
    return 1; // add 1 for temp
  }

  size_t nb_func() const
  {
    return 1;
  }

  size_t iv2ind(size_t j, size_t v) const
  {
    return v;
  }

  int rhs_mat_number_;
  int lhs_mat_number_;
  // int inflow_mat_num_;
  // int outflow_mat_num_;
  // int foil1_mat_num_;
  // int foil2_mat_num_;
  // DataType inflow_vel_x_;
  DataType time_;
};

// Functor used for the local assembly of the stiffness matrix and load vector.
class LocalFlowAssembler : private AssemblyAssistant<DIM, DataType>
{
public:
  void set_parameters(DataType theta, DataType dt, DataType Ra, DataType Pr)
  {
    this->Ra = Ra;
    this->Pr = Pr;
    this->dt_ = dt;
    this->theta_ = theta;
  }

  void set_newton_solution(VectorType const *newton_sol)
  {
    prev_newton_sol_ = newton_sol;
  }

  void set_previous_solution(VectorType const *prev_sol)
  {
    prev_time_sol_ = prev_sol;
  }

  // compute local matrix
  // [in]  element:    contains information about current cell
  // [in]  quadrature: quadrature rule to be used for approximating the integrals
  // [out] lm: contribution of the current cell to the global system matrix

  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature,
                  LocalMatrix &lm)
  {
    const bool need_basis_hessians = false;

    // AssemblyAssistant sets up the local FE basis functions for the current cell
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, need_basis_hessians);

    // evaluate last newton iterate at all quadrature points
    this->evaluate_last_newton_and_time_iterate();

    // number of degrees of freedom on current cell
    const size_t num_dof = this->num_dofs_total();

    // number of quadrature points
    const int num_q = this->num_quadrature_points();
    DataType c1 = std::sqrt(Pr / Ra); // @TODO: crosscheck if the division is floating pt.
    DataType c2 = 1.0 / std::sqrt(Pr * Ra);

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      // quadrature weight
      const DataType wq = w(q);

      // volume element of cell transformation
      const DataType dJ = std::abs(this->detJ(q));

      // previous Newton iterate
      // store velocity, velocity jacobian and pressure at current quad point
      Vec<DIM, DataType> vk;               //-> u_k
      Mat<DIM, DIM, DataType> Dvk;         // -> grad(u_k)
      DataType pk = this->sol_ns_[DIM][q]; // -> p_k
      

      Vec<DIM, DataType> Dtk = this->grad_sol_ns_[DIM+1][q];

      for (int v = 0; v != DIM; ++v)
      {
        vk[v] = this->sol_ns_[v][q];
        for (int d = 0; d != DIM; ++d)
        {
          Dvk(v, d) = this->grad_sol_ns_[v][q][d];
        }
      }

      // previous time step
      Vec<DIM, DataType> vn;
      Mat<DIM, DIM, DataType> Dvn;
      for (int v = 0; v != DIM; ++v)
      {
        vn[v] = this->sol_ts_[v][q];
        for (int d = 0; d != DIM; ++d)
        {
          Dvn(v, d) = this->grad_sol_ts_[v][q][d];
        }
      }

      // loop over test DOFs <-> test function v
      for (int i = 0; i < num_dof; ++i)
      {
        // get test function values for flux variables
        Vec<DIM, DataType> phiV_i;
        Mat<DIM, DIM, DataType> DphiV_i;

        for (size_t var = 0; var < DIM; ++var)
        {
          phiV_i[var] = this->Phi(i, q, var);
          for (int d = 0; d != DIM; ++d)
          {
            DphiV_i(var, d) = this->grad_Phi(i, q, var)[d];
          }
        }
        // get test function values for pressure variable
        DataType phiP_i = this->Phi(i, q, DIM);
        Vec<DIM, DataType> DphiP_i = this->grad_Phi(i, q, DIM);

        // get test function values for temperature variable
        // @TODO change dim for vec<DIM, ... >??
        DataType phiT_i = this->Phi(i, q, DIM + 1);
        Vec<DIM, DataType> DphiT_i = this->grad_Phi(i, q, DIM + 1);

        // loop over trrial DOFs <-> trial function u
        for (int j = 0; j < num_dof; ++j)
        {
          // get ansatz function values for velocity variables
          Vec<DIM, DataType> phiV_j;
          Mat<DIM, DIM, DataType> DphiV_j;

          for (size_t var = 0; var < DIM; ++var)
          {
            phiV_j[var] = this->Phi(j, q, var);
            for (int d = 0; d != DIM; ++d)
            {
              DphiV_j(var, d) = this->grad_Phi(j, q, var)[d];
            }
          }
          // get ansatz function values for pressure variable
          DataType phiP_j = this->Phi(j, q, DIM);
          Vec<DIM, DataType> DphiP_j = this->grad_Phi(j, q, DIM);

          // get ansatz function values for Temperature variable
          DataType phiT_j = this->Phi(j, q, DIM + 1);
          Vec<DIM, DataType> DphiT_j = this->grad_Phi(j, q, DIM + 1);

          // --------------------------------------------------------------------------
          // ----- start assembly of individual terms in variational formulation ------

          // begin with velocity - velocity part
          if (this->first_dof_for_var(0) <= i && i < this->last_dof_for_var(DIM - 1) && this->first_dof_for_var(0) <= j && j < this->last_dof_for_var(DIM - 1))
          {
            DataType l0 = 0.;
            DataType l1 = 0.;
            DataType l2 = 0.;
            DataType l3 = 0.;

            // TODO EXERCISE A
            // Navier Stokes

            l0 = dot(phiV_j, phiV_i);         // check
            l1 = -c1 * dot(DphiV_j, DphiV_i); // check

            for (int v = 0; v != DIM; ++v)
            {
              for (int d = 0; d != DIM; ++d)
              {
                // v_k * grad(v)^T * w
                l2 += vk[d] * DphiV_j(v, d) * phiV_i[v]; // check

                // v * grad(v_k)^T * w
                l3 += phiV_j[d] * Dvk(v, d) * phiV_i[v]; // check
              }
            }

            lm(i, j) += wq * ((l0 / dt_) + theta_ * (l1 + l2 + l3)) * dJ;

            // END EXERCSIE A
          }

          // velocity - pressure part
          if (this->first_dof_for_var(0) <= i && i < this->last_dof_for_var(DIM - 1) && this->first_dof_for_var(DIM) <= j && j < this->last_dof_for_var(DIM))
          {
            // TODO EXERCISE A
            // - p * div(v)
            DataType l4 = -phiP_j * trace(DphiV_i); // check
            lm(i, j) += wq * l4 * dJ;
            // END EXERCISE A
          }

          // pressure - velocity part
          if (this->first_dof_for_var(DIM) <= i && i < this->last_dof_for_var(DIM) && this->first_dof_for_var(0) <= j && j < this->last_dof_for_var(DIM - 1))
          {
            // TODO EXERCISE A
            // div(u) * q
            DataType l5 = phiP_i * trace(DphiV_j); // check
            lm(i, j) += wq * l5 * dJ;
            // END EXERCISE A
          }

          // velocity (test) - temperature
          if (this->first_dof_for_var(0) <= i && i < this->last_dof_for_var(DIM - 1) && this->first_dof_for_var(DIM + 1) <= j && j < this->last_dof_for_var(DIM + 1))
          {
            // @TODO
            Vec<DIM, DataType> grav;
            grav[0] = 0;
            grav[1] = -1;
            // -(theta j, w)
            DataType l6 = -dot(phiT_j * grav, phiV_i);
            lm(i, j) += wq * l6 * dJ;
          }

          // temperature - temperature
          if (this->first_dof_for_var(DIM + 1) <= i && i < this->last_dof_for_var(DIM + 1) && this->first_dof_for_var(DIM + 1) <= j && j < this->last_dof_for_var(DIM + 1))
          {
            // @TODO
            // 1/k (theta, z)
            // DataType l7 = dot(phiT_j, phiT_i);
            DataType l7 = phiT_j * phiT_i;
            // c2 (div theta, div z)
            DataType l8 = c2 * dot(DphiT_j, DphiT_i);
            lm(i, j) += wq * ((l7 / dt_) + l8) * dJ;
          }

          // temperature - velocity

          if (this->first_dof_for_var(DIM + 1) <= i && i < this->last_dof_for_var(DIM + 1) && this->first_dof_for_var(0) <= j && j < this->last_dof_for_var(DIM - 1))
          {
            // @TODO
            DataType l9 = 0.;
            DataType l10 = 0.;

            // for (int v = 0; v != DIM; ++v)
            // {
              for (int d = 0; d != DIM; ++d)
              { // change DPhiV_j
                // v_k * grad(theta)^T * z
                l9 += vk[d] * DphiT_j[d] * phiT_i; // check

                // v * grad(theta_k)^T * z
                l10 += phiV_j[d] * Dtk[d] * phiT_i; // check
              }
            // }

            lm(i, j) += wq * (l9+l10) * dJ;
          }
        }
      }
    }
  }

  // compute local right hand side vector
  // [in]  element:    contains information about current cell
  // [in]  quadrature: quadrature rule to be used for approximating the integrals
  // [out] lv: contribution of the current cell to the global system right hand side
  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature,
                  LocalVector &lv)
  {
    const bool need_basis_hessians = false;

    // AssemblyAssistant sets up the local FE basis functions for the current cell
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, need_basis_hessians);

    // evaluate last newton iterate at all quadrature points
    this->evaluate_last_newton_and_time_iterate();

    // number of degrees of freedom on current cell
    const size_t num_dof = this->num_dofs_total();

    // number of quadrature points
    const int num_q = this->num_quadrature_points();
    DataType c1 = std::sqrt(Pr / Ra);
    DataType c2 = 1.0 / std::sqrt(Pr * Ra);

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      // quadrature weight
      const DataType wq = w(q);

      // volume element of cell transformation
      const DataType dJ = std::abs(this->detJ(q));

      // previous Newton iterate
      // store velocity, velocity jacobian and pressure at current quad point
      Vec<DIM, DataType> vk;               //-> u_k
      Mat<DIM, DIM, DataType> Dvk;         // -> grad(u_k)
      DataType pk = this->sol_ns_[DIM][q]; // -> p_k
      // temperature
      DataType theta_k = this->sol_ns_[DIM+1][q];
      
      Vec<DIM, DataType> Dtk = this->grad_sol_ns_[DIM+1][q];



      for (int v = 0; v != DIM; ++v)
      {
        vk[v] = this->sol_ns_[v][q];
        for (int d = 0; d != DIM; ++d)
        {
          Dvk(v, d) = this->grad_sol_ns_[v][q][d];
        }
      }

      // previous time step
      Vec<DIM, DataType> vn;
      Mat<DIM, DIM, DataType> Dvn;
      // temperature
      DataType theta_ts = this-> sol_ts_[DIM+1][q];

      for (int v = 0; v != DIM; ++v)
      {
        vn[v] = this->sol_ts_[v][q];
        for (int d = 0; d != DIM; ++d)
        {
          Dvn(v, d) = this->grad_sol_ts_[v][q][d];
        }
      }

      // loop over test DOFs <-> test function v
      for (int i = 0; i < num_dof; ++i)
      {
        // get test function values for flux variables
        Vec<DIM, DataType> phiV_i;
        Mat<DIM, DIM, DataType> DphiV_i;

        for (size_t var = 0; var < DIM; ++var)
        {
          phiV_i[var] = this->Phi(i, q, var);
          for (int d = 0; d != DIM; ++d)
          {
            DphiV_i(var, d) = this->grad_Phi(i, q, var)[d];
          }
        }
        // get test function values for pressure variable
        DataType phiP_i = this->Phi(i, q, DIM);
        Vec<DIM, DataType> DphiP_i = this->grad_Phi(i, q, DIM);

        // get test function values for temperature variable
        DataType phiT_i = this->Phi(i, q, DIM + 1);
        Vec<DIM, DataType> DphiT_i = this->grad_Phi(i, q, DIM + 1);

        // momentum equation
        if (this->first_dof_for_var(0) <= i && i < this->last_dof_for_var(DIM - 1))
        {
          DataType l0 = 0.;
          DataType l1_n = 0.;
          DataType l1_k = 0.;
          // DataType phiT_j = this->Phi(j, q, DIM+1);
          DataType l2_n = 0.;
          DataType l2_k = 0.;
          DataType lf = 0.;
          DataType lp = 0.;


          DataType t0 = 0.;
          
          // TODO EXERCISE B
          l0 = dot(vk - vn, phiV_i);
          l1_k = c1 * dot(Dvk, DphiV_i);
          // forcing function
          Vec<DIM, DataType> grav; // make grave class var
          grav[0] = 0;
          grav[1] = -1;
          lf = -dot(theta_k*grav, phiV_i);
          // pressure
          lp = -pk * trace(DphiV_i);

          for (int v = 0; v != DIM; ++v)
          {
            for (int d = 0; d != DIM; ++d)
            {
              // u_k * grad(u_k)^T * v
              l2_k += vk[d] * Dvk(v, d) * phiV_i[v];

            }
          }



          // lt = dot(phiT_j, phiV_i);

          lv[i] += wq * ((l0/dt_) + theta_ * (l1_k + l2_k) + dt_ * (lp + lf)) * dJ;
          
          // Crank Nikolson :
          // l1_n = c1 * dot(Dvn, DphiV_i);
          // l2_n += vn[d] * Dvn(v, d) * phiV_i[v];
          // dt_ * (1. - theta_) * (l1_n + l2_n) +
          // END EXERCISE B
        }

        // mass equation
        if (this->first_dof_for_var(DIM) <= i && i < this->last_dof_for_var(DIM))
        {
          // TODO EXERCISE B
          lv[i] += wq * phiP_i * trace(Dvk) * dJ;
          // END EXERCISE B
        }

        // temperature equation
        if (this->first_dof_for_var(DIM + 1) <= i && i < this->last_dof_for_var(DIM + 1))
        {
          // @TODO
          // DataType lt_0 = - c2 * (DphiT_j, DphiT_i);
          // lv[i] += wq * lt_0 * dJ;
          DataType lt0 = 0.;
          DataType lt1 = 0.;
          DataType lt2 = 0.;

          lt0 = (theta_k - theta_ts) * phiT_i;
          lt1 = c2 * dot(Dtk, DphiT_i);

            for (int d = 0; d != DIM; ++d)
            {
              // u_k * grad(u_k)^T * v
              lt2 += vk[d] * Dtk[d] * phiT_i;

              // lt_n += 
            }
            lv[i]+= wq * ((lt0/dt_) + lt1 +lt2) * dJ;          

          
        }
      }
    }
  }

  Vec<DIM, DataType> force(Vec<DIM, DataType> pt)
  {
    Vec<DIM, DataType> f;
    f[DIM - 1] = this->f_;
    return f;
  }

  // evaluate FE function corresponding to vector prev_newton_sol
  // at all quadrature points
  void evaluate_last_newton_and_time_iterate()
  {
    // assert (prev_newton_sol_ != nullptr);
    // assert (prev_time_sol_ != nullptr);
    
    // compute velocity solution values of previous time step u_(n-1) at each quadrature point xq:
    // sol_ts_[v][q] = u_(n-1)_v (xq)
    // v = velocity component
    for (int v = 0; v != DIM + 2; ++v) // @TODO change to DIM+2
    {
      sol_ns_[v].clear();
      grad_sol_ns_[v].clear();

      this->evaluate_fe_function(*prev_newton_sol_, v, sol_ns_[v]);
      this->evaluate_fe_function_gradients(*prev_newton_sol_, v, grad_sol_ns_[v]);
    }
    for (int v = 0; v != DIM + 2; ++v) // @TODO change to DIM+2
    {
      this->sol_ts_[v].clear();
      this->grad_sol_ts_[v].clear();
      this->evaluate_fe_function(*prev_time_sol_, v, sol_ts_[v]);
      this->evaluate_fe_function_gradients(*prev_time_sol_, v, grad_sol_ts_[v]);
    }
  }
  // TODO check if required to change these as well?
  FunctionValues<DataType> sol_ns_[DIM + 2];                // solution at previous newton step
  FunctionValues<Vec<DIM, DataType>> grad_sol_ns_[DIM + 2]; // gradient of solution at previous newton step

  VectorType const *prev_newton_sol_;

  VectorType const *prev_time_sol_;
  FunctionValues<DataType> sol_ts_[DIM+2];                // velocity solution at previous time step
  FunctionValues<Vec<DIM, DataType>> grad_sol_ts_[DIM+2]; // gradient of velocity solution at previous time step

  DataType Ra;
  DataType Pr;
  DataType theta_;
  DataType dt_;
  DataType nu_;
  DataType f_;
};

// @TODO Remove Force integral
template <int DIM, class LAD>
class ForceIntegral : private AssemblyAssistant<DIM, typename LAD::DataType>
{
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

public:
  ForceIntegral(DataType mu,
                const VectorType *sol,
                int force_direction,
                int obstacle_mat_num)
      : mu_(mu), force_dir_(force_direction), sol_(sol), obstacle_mat_num_(obstacle_mat_num)
  {
  }

  // compute local scalar contributions from boundary facets
  // [in]  element:    contains information about current cell
  // [in]  facet_number: local index of element facet
  // [in]  quadrature: quadrature rule to be used for approximating the integrals
  // [out] ls[facet_number]: contribution of the current facet to the global scalar
  void operator()(const Element<DataType, DIM> &element,
                  int facet_number,
                  const Quadrature<DataType> &quadrature,
                  std::vector<DataType> &ls)
  {
    const bool need_basis_hessians = false;

    AssemblyAssistant<DIM, DataType>::initialize_for_facet(element,
                                                           quadrature,
                                                           facet_number,
                                                           need_basis_hessians);

    this->recompute_function_values();

    // get material number of current facet
    IncidentEntityIterator facet = element.get_cell().begin_incident(DIM - 1);
    for (int i = 0; i < facet_number; ++i, ++facet)
    {
    }

    const int material_number = facet->get_material_number();

    // number of quadrature points
    const int num_q = this->num_quadrature_points();

    // this variable stores the integral over the current facet
    DataType facet_integral = 0.;

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      // quadrature weight
      const DataType wq = this->w(q);

      // surface element of cell transformation
      const DataType dS = std::abs(this->ds(q));

      // surface normal
      // multiply by -1,
      // because we want to have the outward normal w.r.t. the obstacle and not w.r.t. the compuational domain
      const Vec<DIM, DataType> nq = (-1.) * this->n(q);

      if (material_number != obstacle_mat_num_)
      {
        // consider only obstacle surface
        continue;
      }

      // force vector
      // t = sigma * n
      // sigma = mu * (gradv + grad v ^T) - p I

      Vec<DIM, DataType> t;
      Mat<DIM, DIM, DataType> sigma;
      for (int v = 0; v != DIM; ++v)
      {
        for (int d = 0; d != DIM; ++d)
        {
          sigma(v, d) += mu_ * (grad_vel_[v][q][d] + grad_vel_[d][q][v]);
        }
        sigma(v, v) -= press_[q];
      }
      sigma.VectorMult(nq, t);

      facet_integral += wq * t[force_dir_] * dS;
    }

    ls[facet_number] = facet_integral;
  }

private:
  void recompute_function_values()
  {
    for (int d = 0; d < DIM; ++d)
    {
      vel_[d].clear();
      grad_vel_[d].clear();

      this->evaluate_fe_function(*sol_, d, vel_[d]);
      this->evaluate_fe_function_gradients(*sol_, d, grad_vel_[d]);
    }
    press_.clear();
    this->evaluate_fe_function(*sol_, DIM, press_);
  }

  DataType mu_;
  int obstacle_mat_num_, force_dir_;
  FunctionValues<DataType> vel_[DIM], press_;
  FunctionValues<Vec<DIM, DataType>> grad_vel_[DIM];
  const VectorType *sol_;
};

template <int DIM, class LAD>
class DivergenceIntegral : private AssemblyAssistant<DIM, typename LAD::DataType>
{
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

public:
  DivergenceIntegral(const VectorType &sol)
      : sol_(sol)
  {
  }

  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature,
                  DataType &l2_div)
  {
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, false);
    for (int v = 0; v != DIM; ++v)
    {
      grad_vel_[v].clear();
      this->evaluate_fe_function_gradients(sol_, v, grad_vel_[v]);
    }

    const int num_q = this->num_quadrature_points();

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      const DataType wq = this->w(q);
      DataType dJ = std::abs(this->detJ(q));
      DataType div = 0.;
      for (int v = 0; v != DIM; ++v)
      {
        div += grad_vel_[v][q][v];
      }
      l2_div += wq * div * div * dJ;
    }
  }

  const VectorType &sol_;
  FunctionValues<Vec<DIM, DataType>> grad_vel_[DIM];
};

template <int DIM, class LAD>
class PressureIntegral : private AssemblyAssistant<DIM, typename LAD::DataType>
{
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

public:
  PressureIntegral(const VectorType &sol)
      : sol_(sol)
  {
  }

  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature,
                  DataType &pressure)
  {
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, false);
    this->evaluate_fe_function(sol_, DIM, p_);

    const int num_q = this->num_quadrature_points();

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      const DataType wq = this->w(q);
      DataType dJ = std::abs(this->detJ(q));
      pressure += wq * p_[q] * dJ;
    }
  }
  const VectorType &sol_;
  FunctionValues<DataType> p_;
};

template <int DIM, class LAD>
class VolumeIntegral : private AssemblyAssistant<DIM, typename LAD::DataType>
{
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

public:
  VolumeIntegral()
  {
  }

  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature,
                  DataType &vol)
  {
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, false);
    const int num_q = this->num_quadrature_points();

    // loop over quadrature points
    for (int q = 0; q < num_q; ++q)
    {
      const DataType wq = this->w(q);
      vol += wq * std::abs(this->detJ(q));
    }
  }
};
