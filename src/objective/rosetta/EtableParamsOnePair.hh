#ifndef INCLUDED_objective_rosetta_EtableParamsOnePair_hh
#define INCLUDED_objective_rosetta_EtableParamsOnePair_hh

#include "objective/rosetta/CubicPolynomial.hh"

namespace objective { 
namespace rosetta {


struct ExtraQuadraticRepulsion
{
	double xlo;
	double xhi;
	double slope;
	double extrapolated_slope;
	double ylo;
	ExtraQuadraticRepulsion() : xlo(0),xhi(0),slope(0),extrapolated_slope(0),ylo(0) {}
	ExtraQuadraticRepulsion(double a, double b, double c, double d, double e)
	 : xlo(a),xhi(b),slope(c),extrapolated_slope(d),ylo(e)
	{}
};


/// @brief %EtableParamsOnePair describes all of the parameters for a particular
/// pair of atom types necessary to evaluate the Lennard-Jones and LK solvation
/// energies.
struct EtableParamsOnePair
{
	typedef double Real;
	typedef CubicPolynomialParamsBase<Real> CubicPolynomial;

	Real maxd2;
	bool hydrogen_interaction;
	Real ljrep_linear_ramp_d2_cutoff;
	Real lj_r6_coeff;
	Real lj_r12_coeff;
	Real lj_switch_intercept;
	Real lj_switch_slope;
	Real lj_minimum;
	Real lj_val_at_minimum;
	Real ljatr_cubic_poly_xlo;
	Real ljatr_cubic_poly_xhi;
	CubicPolynomial ljatr_cubic_poly_parameters;
	ExtraQuadraticRepulsion ljrep_extra_repulsion;
	bool ljrep_from_negcrossing;

	Real lk_coeff1;
	Real lk_coeff2;
	Real lk_min_dis2sigma_value;
	Real fasol_cubic_poly_close_start;
	Real fasol_cubic_poly_close_end;

	Real fasol_cubic_poly_far_start;
	Real fasol_cubic_poly_far_end;

	Real lj_radius_1; // lj radii of atoms stored for analytic lk evaluation
	Real lj_radius_2; 

	Real lk_inv_lambda2_1; // lk params of atoms stored for analytic lk evaluation
	Real lk_inv_lambda2_2;

	CubicPolynomial fasol_cubic_poly_close;  // mututal desolvation of atoms 1 and 2
	Real fasol_cubic_poly_close_flat; // the fixed value used for distances beneath fasol_cubic_poly_close_start

	CubicPolynomial fasol_cubic_poly_far;    // mututal desolvation of atoms 1 and 2

	CubicPolynomial fasol_cubic_poly1_close; // desolvation of atom 1 by atom 2
	Real fasol_cubic_poly1_close_flat;       // the fixed value used for distances beneath fasol_cubic_poly_close_start
	CubicPolynomial fasol_cubic_poly1_far;   // desolvation of atom 1 by atom 2

	CubicPolynomial fasol_cubic_poly2_close; // desolvation of atom 2 by atom 1
	Real fasol_cubic_poly2_close_flat;       // the fixed value used for distances beneath fasol_cubic_poly_close_start
	CubicPolynomial fasol_cubic_poly2_far;   // desolvation of atom 2 by atom 1
	Real ljatr_final_weight;
	Real fasol_final_weight;

	EtableParamsOnePair() :

		maxd2(35.9999),
		hydrogen_interaction(0),
		ljrep_linear_ramp_d2_cutoff(0),
		lj_r6_coeff(0),
		lj_r12_coeff(0),
		lj_switch_intercept(0),
		lj_switch_slope(0),
		lj_minimum(0),
		lj_val_at_minimum(0),
		ljatr_cubic_poly_xlo(0),
		ljatr_cubic_poly_xhi(0),
		ljatr_cubic_poly_parameters(0,0,0,0),
		ljrep_extra_repulsion(0,0,0,0,0),
		ljrep_from_negcrossing(0),

		lk_coeff1(0),
		lk_coeff2(0),
		lk_min_dis2sigma_value(0),
		fasol_cubic_poly_close_start(0),
		fasol_cubic_poly_close_end(0),

		fasol_cubic_poly_far_start(0),
		fasol_cubic_poly_far_end(0),

		lj_radius_1(0),
		lj_radius_2(0),

		lk_inv_lambda2_1(0),
		lk_inv_lambda2_2(0),

		fasol_cubic_poly_close(0,0,0,0),
		fasol_cubic_poly_close_flat(0),

		fasol_cubic_poly_far(0,0,0,0),

		fasol_cubic_poly1_close(0,0,0,0),
		fasol_cubic_poly1_close_flat(0),
		fasol_cubic_poly1_far(0,0,0,0),

		fasol_cubic_poly2_close(0,0,0,0),
		fasol_cubic_poly2_close_flat(0),
		fasol_cubic_poly2_far(0,0,0,0),
		ljatr_final_weight(1),
		fasol_final_weight(1)
	{}


};

}
}

#endif

