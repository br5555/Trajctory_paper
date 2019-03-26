#include "Trajectory_without_drag.h"
#include <math.h>



#include <boost/math/tools/roots.hpp>
using boost::math::policies::policy;
using boost::math::tools::newton_raphson_iterate;
using boost::math::tools::halley_iterate;
using boost::math::tools::eps_tolerance; // Binary functor for specified number of bits.
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::toms748_solve;

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h> 

#include <boost/math/special_functions/next.hpp>

#include <tuple>
#include <utility> // pair, make_pair

//] [/root_finding_headers]

#include <iostream>
using std::cout; using std::endl;
#include <iomanip>
using std::setw; using std::setprecision;
#include <limits>
using std::numeric_limits;





/// ###################################################################

template <class T>
struct fifth_functor_2deriv
{ // Functor returning both 1st and 2nd derivatives.
	fifth_functor_2deriv(T const& a1, T const& a2, T const& a3, T const& a4, T const& a5)
		: m_a5(a5), m_a4(a4), m_a3(a3), m_a2(a2), m_a1(a1)
	{ // Constructor stores value to find root of, for example:
		std::cout << " a1 " << m_a1 << std::endl;
		std::cout << " a2 " << m_a2 << std::endl;
		std::cout << " a3 " << m_a3 << std::endl;
		std::cout << " a4 " << m_a4 << std::endl;
		std::cout << " a5 " << m_a5 << std::endl;

	}

	// using boost::math::tuple; // to return three values.
	std::tuple<T, T, T> operator()(T const& x)
	{ // Return both f(x) and f'(x) and f''(x).
		T fx = m_a5 + m_a4 * x + m_a3 * x *x + m_a2 * x *x *x + m_a1 * x*x *x*x; // Difference (estimate x^3 - value).

		T dx = m_a4 + 2 * m_a3  *x + 3 * m_a2  *x *x + 4 * m_a1 *x *x*x; // 1st derivative = 5x^4.
		T d2x = 2 * m_a3 + 3 * 2 * m_a2   *x + 4 * 3 * m_a1  *x*x; // 2nd derivative = 20 x^3
		return std::make_tuple(fx, dx, d2x); // 'return' fx, dx and d2x.
	}
private:
	T  m_a5, m_a4, m_a3, m_a2, m_a1; // to be 'fifth_rooted'.
}; // struct fifth_functor_2deriv

//] [/root_finding_fifth_functor_2deriv]


/*`Our fifth function is now:*/

//[root_finding_fifth_2deriv

template <class T>
T fifth_2deriv(T a1, T a2, T a3, T a4, T a5)
{ // return fifth root of x using 1st and 2nd derivatives and Halley.
	using namespace std;  // Help ADL of std functions.
	using namespace boost::math; // halley_iterate

	int exponent;
	srand(time(NULL));



	frexp(-1, &exponent); // Get exponent of z (ignore mantissa).
	T guess = 10.0;// ldexp(1., exponent / 5); // Rough guess is to divide the exponent by three.
	T min = ldexp(-1, exponent / 5); // Minimum possible value is half our guess.
	T max = ldexp(200., exponent / 5); // Maximum possible value is twice our guess.

	std::cout << "min " << min << std::endl;
	std::cout << "max " << max << std::endl;

	int digits = std::numeric_limits<T>::digits / 2; // Half maximum possible binary digits accuracy for type T.
	const boost::uintmax_t maxit = 500;
	boost::uintmax_t it = maxit;
	T result = halley_iterate(fifth_functor_2deriv<T>(a1, a2, a3, a4, a5), guess, min, max, digits, it);
	// Can show how many iterations (updated by halley_iterate).
	cout << it << " iterations (from max of " << maxit << ")" << endl;

	return result;
}


Trajectory_without_drag::Trajectory_without_drag()
{
	m_jerk_max = 5;
	m_jerk_min = -5;
	m_acc_max = 5.0;
	m_acc_min = -5.0;
	m_vel_max = 10;
	m_vel_min = -10;
	m_vel_start = 0.001;
	m_vel_end = 0;
	m_p_start = 0;
	m_p_end = 10.3;
	m_acc_start = 0;
	m_acc_end = 0;
}

Trajectory_without_drag::Trajectory_without_drag(double jerk_max, double  jerk_min, double  acc_max, double  acc_min, double vel_max, double  vel_min, double  vel_start, double  vel_end, double  p_start, double  p_end
	, double  acc_start, double  acc_end) : m_jerk_max(jerk_max), m_jerk_min(jerk_min), m_acc_max(acc_max), m_acc_min(acc_min), m_vel_max(vel_max), m_vel_min(vel_min),
	m_vel_start(vel_start), m_vel_end(vel_end), m_p_start(p_start), m_p_end(p_end), m_acc_start(acc_start), m_acc_end(acc_end) {}

Trajectory_without_drag::~Trajectory_without_drag()
{
}

bool Trajectory_without_drag::calculate_trajectory_case1(double t[]) {


	double t1, t2, t3, t4, t5, t6, t7;

	t1 = (m_acc_max - m_acc_start) / m_jerk_max;

	t2 = (-(m_acc_max*m_acc_max)*m_jerk_min + (m_acc_max*m_acc_max)*m_jerk_max + (m_acc_start*m_acc_start)*m_jerk_min + m_jerk_min * m_jerk_max*m_vel_max*2.0 - m_jerk_min * m_jerk_max*m_vel_start*2.0) / (m_acc_max*m_jerk_min*m_jerk_max*2.0);

	t3 = -m_acc_max / m_jerk_min;

	double A = m_acc_min * (m_acc_max*m_acc_max*m_acc_max*m_acc_max)*(m_jerk_min*m_jerk_min) - (m_acc_min*m_acc_min*m_acc_min*m_acc_min)*m_acc_max*(m_jerk_min*m_jerk_min) - m_acc_min * (m_acc_max*m_acc_max*m_acc_max*m_acc_max)*(m_jerk_max*m_jerk_max) + (m_acc_min*m_acc_min*m_acc_min*m_acc_min)*m_acc_max*(m_jerk_max*m_jerk_max) - m_acc_min * (m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*3.0;

	double B = (m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_max*(m_jerk_min*m_jerk_min)*3.0 - (m_acc_end*m_acc_end*m_acc_end)*m_acc_min*m_acc_max*(m_jerk_min*m_jerk_min)*8.0 + m_acc_min * m_acc_max*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*8.0;

	double C = (m_acc_end*m_acc_end)*(m_acc_min*m_acc_min)*m_acc_max*(m_jerk_min*m_jerk_min)*6.0 - m_acc_min * (m_acc_max*m_acc_max)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*6.0 + m_acc_min * (m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_max*m_vel_max)*1.2E1;

	double D = m_acc_max * (m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1 - m_acc_max * (m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_max*m_vel_max)*1.2E1 - m_acc_min * (m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1;

	double E = m_acc_min * m_acc_max*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*-2.4E1 + m_acc_min * m_acc_max*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1;

	double F = m_acc_min * (m_acc_max*m_acc_max)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_max*-1.2E1 + (m_acc_min*m_acc_min)*m_acc_max*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_max*1.2E1;

	double G = (m_acc_min*m_acc_min)*m_acc_max*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*-1.2E1 + m_acc_min * (m_acc_max*m_acc_max)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.2E1;

	double H = (m_acc_end*m_acc_end)*m_acc_max*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*-1.2E1 + m_acc_min * (m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.2E1;

	double I = m_acc_end * m_acc_min*m_acc_max*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1 - m_acc_min * m_acc_max*m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1;

	double J = m_acc_min * m_acc_max*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_max*2.4E1;

	t4 = -(A + B + C + D + E + F + G + H + I) / J;

	t5 = m_acc_min / m_jerk_min;

	t6 = (((m_acc_end*m_acc_end)*m_jerk_min - (m_acc_min*m_acc_min)*m_jerk_min + (m_acc_min*m_acc_min)*m_jerk_max - m_jerk_min * m_jerk_max*m_vel_end*2.0 + m_jerk_min * m_jerk_max*m_vel_max*2.0)*(-1.0 / 2.0)) / (m_acc_min*m_jerk_min*m_jerk_max);

	t7 = (m_acc_end - m_acc_min) / m_jerk_max;


	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case2(t);
	}

	return true;
}

bool Trajectory_without_drag::calculate_trajectory_case2(double t[]) {

	double t1, t2, t3, t4, t5, t6, t7;
	double v1_tmp;
	double v1;
	double p1;
	double v6;
	double p6;
	double t6_tmp_tmp;
	double b_t6_tmp_tmp;
	double c_t6_tmp_tmp;
	double d_t6_tmp_tmp;
	double t6_tmp_tmp_tmp;
	double e_t6_tmp_tmp;
	double f_t6_tmp_tmp;
	double g_t6_tmp_tmp;
	double h_t6_tmp_tmp;
	double i_t6_tmp_tmp;
	double j_t6_tmp_tmp;
	double k_t6_tmp_tmp;
	double l_t6_tmp_tmp;
	double m_t6_tmp_tmp;
	double b_t6_tmp_tmp_tmp;
	double n_t6_tmp_tmp;
	double o_t6_tmp_tmp;
	double p_t6_tmp_tmp;
	double q_t6_tmp_tmp;
	double r_t6_tmp_tmp;
	t1 = (m_acc_max - m_acc_start) / m_jerk_max;
	t4 = 0.0;
	t3 = -m_acc_max / m_jerk_min;
	t7 = (m_acc_end - m_acc_min) / m_jerk_max;
	t5 = m_acc_min / m_jerk_min;
	v1_tmp = 0.5 * (t1 * t1);
	v1 = (m_vel_start + t1 * m_acc_start) + v1_tmp * m_jerk_max;
	p1 = 0.5 * (t7 * t7);
	v6 = (m_vel_end - t7 * m_acc_min) - p1 * m_jerk_max;
	p6 = ((m_p_end - t7 * v6) - p1 * m_acc_min) - 0.16666666666666666 * pow
	(t7, 3.0) * m_jerk_max;
	p1 = ((m_p_start + t1 * m_vel_start) + v1_tmp * m_acc_start) + 0.16666666666666666 *
		pow(t1, 3.0) * m_jerk_max;
	v1_tmp = m_acc_max * m_acc_max;
	t6_tmp_tmp = 4.0 * v1_tmp;
	b_t6_tmp_tmp = m_acc_min * m_acc_min;
	c_t6_tmp_tmp = t6_tmp_tmp * b_t6_tmp_tmp;
	d_t6_tmp_tmp = 8.0 * m_acc_max * b_t6_tmp_tmp;
	t6_tmp_tmp_tmp = 8.0 * v1_tmp;
	e_t6_tmp_tmp = t6_tmp_tmp_tmp * m_acc_min;
	f_t6_tmp_tmp = v1 * v1;
	g_t6_tmp_tmp = 4.0 * m_acc_max * m_acc_min;
	h_t6_tmp_tmp = v6 * v6;
	i_t6_tmp_tmp = t3 * t3;
	j_t6_tmp_tmp = pow(t3, 3.0);
	k_t6_tmp_tmp = 4.0 * m_acc_max * b_t6_tmp_tmp * m_jerk_min;
	l_t6_tmp_tmp = e_t6_tmp_tmp * m_jerk_min;
	m_t6_tmp_tmp = pow(t5, 3.0);
	b_t6_tmp_tmp_tmp = m_jerk_min * m_jerk_min;
	n_t6_tmp_tmp = m_acc_max * m_acc_min * b_t6_tmp_tmp_tmp;
	o_t6_tmp_tmp = t5 * t5;
	p_t6_tmp_tmp = 4.0 * m_acc_max * m_acc_max;
	q_t6_tmp_tmp = p_t6_tmp_tmp * m_acc_min * m_jerk_min;
	r_t6_tmp_tmp = 2.0 * m_acc_max * m_acc_min;
	t6_tmp_tmp = ((((((((((((((((((((((((4.0 * b_t6_tmp_tmp * f_t6_tmp_tmp +
		t6_tmp_tmp * h_t6_tmp_tmp) + c_t6_tmp_tmp * i_t6_tmp_tmp) + c_t6_tmp_tmp *
		o_t6_tmp_tmp) - d_t6_tmp_tmp * p1) + e_t6_tmp_tmp * p1) + d_t6_tmp_tmp * p6)
		- e_t6_tmp_tmp * p6) - g_t6_tmp_tmp * f_t6_tmp_tmp) - g_t6_tmp_tmp *
		h_t6_tmp_tmp) - p_t6_tmp_tmp * b_t6_tmp_tmp * i_t6_tmp_tmp) + 4.0 * m_acc_max *
		v1_tmp * m_acc_min * i_t6_tmp_tmp) - t6_tmp_tmp * m_acc_max * m_acc_min *
		i_t6_tmp_tmp) - k_t6_tmp_tmp * j_t6_tmp_tmp / 3.0) - l_t6_tmp_tmp *
		j_t6_tmp_tmp / 3.0) + n_t6_tmp_tmp * pow(t3, 4.0)) - k_t6_tmp_tmp *
		m_t6_tmp_tmp / 3.0) - l_t6_tmp_tmp * m_t6_tmp_tmp / 3.0)
		+ n_t6_tmp_tmp * pow(t5, 4.0)) + t6_tmp_tmp_tmp *
		b_t6_tmp_tmp * t3 * t5) - k_t6_tmp_tmp * i_t6_tmp_tmp *
		t5) - t6_tmp_tmp * m_acc_min * m_jerk_min * t3 * o_t6_tmp_tmp)
		+ r_t6_tmp_tmp * b_t6_tmp_tmp_tmp * i_t6_tmp_tmp *
		o_t6_tmp_tmp) + q_t6_tmp_tmp * j_t6_tmp_tmp) - 8.0 * m_acc_max *
		m_acc_max * b_t6_tmp_tmp * t3 * t5) + q_t6_tmp_tmp * t3 *
		o_t6_tmp_tmp;
	if (t6_tmp_tmp < 0.0) {
		t1 = -1.0;
		t2 = -1.0;
		t3 = -1.0;
		t4 = -1.0;
		t5 = -1.0;
		t6 = -1.0;
		t7 = -1.0;
		return	calculate_trajectory_case3(t);
	}
	else {
		c_t6_tmp_tmp = r_t6_tmp_tmp * t3;
		d_t6_tmp_tmp = m_acc_min * m_jerk_min;
		t6_tmp_tmp_tmp = 2.0 * m_acc_max * v6 - 2.0 * m_acc_min * v6;
		e_t6_tmp_tmp = std::sqrt(t6_tmp_tmp);
		p1 = r_t6_tmp_tmp * t5;
		v1_tmp = d_t6_tmp_tmp * i_t6_tmp_tmp;
		d_t6_tmp_tmp *= o_t6_tmp_tmp;
		b_t6_tmp_tmp = 2.0 * m_acc_min * (m_acc_max - m_acc_min);
		t6 = ((((((t6_tmp_tmp_tmp + e_t6_tmp_tmp) - c_t6_tmp_tmp) + c_t6_tmp_tmp) -
			p1) + v1_tmp) + d_t6_tmp_tmp) / b_t6_tmp_tmp;
		if (t6 < 0.0) {
			t6 = ((((((t6_tmp_tmp_tmp - e_t6_tmp_tmp) - c_t6_tmp_tmp) + c_t6_tmp_tmp)
				- p1) + v1_tmp) + d_t6_tmp_tmp) / b_t6_tmp_tmp;
		}

		if (t6_tmp_tmp < 0.0) {
			t1 = -1.0;
			t2 = -1.0;
			t3 = -1.0;
			t4 = -1.0;
			t5 = -1.0;
			t6 = -1.0;
			t7 = -1.0;
			return calculate_trajectory_case3(t);
		}
		else {
			t2 = (((((v6 - t6 * m_acc_min) - 0.5 * o_t6_tmp_tmp * m_jerk_min) - t3 *
				m_acc_max) - 0.5 * i_t6_tmp_tmp * m_jerk_min) - v1) / m_acc_max;
		}
	}

	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case3(t);
	}

	return true;

}

bool Trajectory_without_drag::calculate_trajectory_case3(double t[]) {
	double t1, t2, t3, t4, t5, t6, t7;
	double v_6_tmp;
	double v_6;
	double v_5_tmp;
	double v_5;
	double b_t1;
	double a_t1;
	double t1_tmp;
	t2 = 0.0;
	t5 = m_acc_min / m_jerk_min;
	t7 = (m_acc_end - m_acc_min) / m_jerk_max;
	v_6_tmp = 0.5 * (t7 * t7);
	v_6 = (m_vel_end - t7 * m_acc_min) - v_6_tmp * m_jerk_max;
	v_5_tmp = 0.5 * (t5 * t5);
	v_5 = (m_vel_max + t5 * 0.0) + v_5_tmp * m_jerk_min;
	t6 = (v_6 - v_5) / m_acc_min;
	b_t1 = m_acc_start - m_acc_start * (m_jerk_max / m_jerk_min);
	a_t1 = 0.5 * m_jerk_max - 0.5 * (m_jerk_max * m_jerk_max / m_jerk_min);
	t1_tmp = b_t1 * b_t1 - 4.0 * a_t1 * ((m_vel_start - m_vel_max) - 0.5 * (m_acc_start *
		m_acc_start / m_jerk_min));
	if (t1_tmp < 0.0) {
		t1 = -1.0;
		t2 = -1.0;
		t3 = -1.0;
		t4 = -1.0;
		t5 = -1.0;
		t6 = -1.0;
		t7 = -1.0;
		return calculate_trajectory_case4(t);
	}
	else {
		t1_tmp = std::sqrt(t1_tmp);
		t1 = (-b_t1 + t1_tmp) / (2.0 * a_t1);
		if (t1 < 0.0) {
			t1 = (-b_t1 - t1_tmp) / (2.0 * a_t1);
		}

		b_t1 = m_acc_start + t1 * m_jerk_max;
		t3 = -b_t1 / m_jerk_min;
		a_t1 = 0.5 * (t1 * t1);
		t4 = (((((((((m_p_end - t7 * v_6) - v_6_tmp * m_acc_min) - 0.16666666666666666
			* pow(t7, 3.0) * m_jerk_max) - t6 * v_5) - 0.5 * (t6 *
				t6) * m_acc_min) - t5 * m_vel_max) - v_5_tmp * 0.0) -
			0.16666666666666666 * pow(t5, 3.0) * m_jerk_min) -
			((((((m_p_start + t1 * m_vel_start) + a_t1 * m_acc_start) +
				0.16666666666666666 * pow(t1, 3.0) * m_jerk_max) + t3 *
				((m_vel_start + t1 * m_acc_start) + a_t1 * m_jerk_max)) + 0.5 * (t3 *
					t3) * b_t1) + 0.16666666666666666 * pow(t3, 3.0) *
				m_jerk_min)) / m_vel_max;
	}

	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case4(t);
	}

	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case4(double t[]) {
	double t1, t2, t3, t4, t5, t6, t7;

	double v_2_tmp;
	double v_2;
	double v_1_tmp;
	double v_1;
	double a_t7;
	double b_t7;
	double t7_tmp;
	t6 = 0.0;
	t1 = (m_acc_max - m_acc_start) / m_jerk_max;
	t3 = -m_acc_max / m_jerk_min;
	v_2_tmp = 0.5 * (t3 * t3);
	v_2 = (m_vel_max - t3 * m_acc_max) - v_2_tmp * m_jerk_min;
	v_1_tmp = 0.5 * (t1 * t1);
	v_1 = (m_vel_start + t1 * m_acc_start) + v_1_tmp * m_jerk_max;
	t2 = (v_2 - v_1) / m_acc_max;
	a_t7 = m_jerk_max * m_jerk_max / (2.0 * m_jerk_min) - 0.5 * m_jerk_max;
	b_t7 = m_acc_end - m_acc_end * m_jerk_max / m_jerk_min;
	t7_tmp = b_t7 * b_t7 - 4.0 * a_t7 * ((m_vel_max - m_vel_end) + 0.5 * (m_acc_end *
		m_acc_end) / m_jerk_min);
	if (t7_tmp < 0.0) {
		t1 = -1.0;
		t2 = -1.0;
		t3 = -1.0;
		t4 = -1.0;
		t5 = -1.0;
		t6 = -1.0;
		t7 = -1.0;
		return calculate_trajectory_case5(t);
	}
	else {
		t7_tmp = std::sqrt(t7_tmp);
		t7 = (-b_t7 + t7_tmp) / (2.0 * a_t7);
		if (t7 < 0.0) {
			t7 = (-b_t7 - t7_tmp) / (2.0 * a_t7);
		}

		a_t7 = m_acc_end - t7 * m_jerk_max;
		t5 = a_t7 / m_jerk_min;
		b_t7 = 0.5 * (t7 * t7);
		t7_tmp = (m_vel_end - t7 * a_t7) - b_t7 * m_jerk_max;
		t4 = (((((((((m_p_end - t7 * t7_tmp) - b_t7 * a_t7) - 0.16666666666666666 *
			pow(t7, 3.0) * m_jerk_max) - 0.0 * t7_tmp) - 0.0 * a_t7)
			- t5 * m_vel_max) - 0.5 * pow(t5, 5.0) * 0.0) -
			0.16666666666666666 * pow(t5, 3.0) * m_jerk_min) -
			((((((((m_p_start + t1 * m_vel_start) + v_1_tmp * m_acc_start) +
				0.16666666666666666 * pow(t1, 3.0) * m_jerk_max) + t2 *
				v_1) + 0.5 * (t2 * t2) * m_acc_max) + t3 * v_2) + v_2_tmp *
				m_acc_max) + 0.16666666666666666 * pow(t3, 3.0) * m_jerk_min))
			/ m_vel_max;
	}
	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case5(t);
	}

	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case5(double t[]) {
	double t1, t2, t3, t4, t5, t6, t7;
	t2 = 0;
	t6 = 0;
	double a_3 = 0;
	double a_4 = a_3;
	double v_3 = m_vel_max;
	double v_4 = v_3;

	double a_t7 = m_jerk_max * (-1.0 / 2.0) + (m_jerk_max*m_jerk_max) / (m_jerk_min*2.0);
	double b_t7 = m_acc_end - (m_acc_end*m_jerk_max) / (m_jerk_min);
	double c_t7 = m_vel_max - m_vel_end + (m_acc_end*m_acc_end) / (m_jerk_min*2.0);

	double t7_tmp = (m_jerk_max*2.0 - ((m_jerk_max*m_jerk_max)*2.0) / m_jerk_min)*(-m_vel_end + m_vel_max + (m_acc_end*m_acc_end) / (m_jerk_min*2.0)) + pow(m_acc_end - (m_acc_end*m_jerk_max) / m_jerk_min, 2.0);

	if (t7_tmp < 0)
	{
		t1 = -1;
		t2 = -1;
		t3 = -1;
		t4 = -1;
		t5 = -1;
		t6 = -1;
		t7 = -1;
		return calculate_trajectory_case6(t);
	}

	t7 = ((b_t7 - sqrt(t7_tmp))*(-1.0 / 2.0)) / a_t7;
	if (t7 < 0)
		t7 = ((b_t7 + sqrt(t7_tmp))*(-1.0 / 2.0)) / a_t7;



	double a_5 = m_acc_end - t7 * m_jerk_max;
	t5 = (a_5) / (m_jerk_min);
	double a_6 = a_5;
	double v_6 = m_vel_end - t7 * a_6 - 0.5*(t7*t7)*m_jerk_max;
	double v_5 = v_6;

	double c_t1 = (m_vel_start - m_vel_max - 0.5*((m_acc_start*m_acc_start) / m_jerk_min));
	double b_t1 = m_acc_start - m_acc_start * (m_jerk_max / m_jerk_min);
	double a_t1 = 0.5*m_jerk_max - 0.5*((m_jerk_max*m_jerk_max) / m_jerk_min);

	double t1_tmp = b_t1 * b_t1 - 4.0*a_t1*c_t1;

	if (t1_tmp < 0)
	{
		t1 = -1;
		t2 = -1;
		t3 = -1;
		t4 = -1;
		t5 = -1;
		t6 = -1;
		t7 = -1;
		return calculate_trajectory_case6(t);
	}

	t1 = (-b_t1 + sqrt(t1_tmp)) / (2.0*a_t1);
	if (t1 < 0)
		t1 = (-b_t1 - sqrt(t1_tmp)) / (2.0*a_t1);




	double a_1 = m_acc_start + t1 * m_jerk_max;
	double a_2 = a_1;
	t3 = -a_2 / m_jerk_min;
	double p1 = m_p_start + t1 * m_vel_start + 0.5*(t1*t1)*m_acc_start + (1.0 / 6.0)*(t1*t1*t1)*m_jerk_max;
	double p_2 = p1;
	double v1 = m_vel_start + t1 * m_acc_start + 0.5*(t1*t1)*m_jerk_max;
	double v_2 = v1;

	double p_3 = p_2 + t3 * v_2 + 0.5*(t3*t3)*a_2 + (1.0 / 6.0)*(t3*t3*t3)*m_jerk_min;


	double p_6 = m_p_end - t7 * v_6 - 0.5*(t7*t7)*a_6 - (1.0 / 6.0)*(t7*t7*t7)*m_jerk_max;
	double p_5 = p_6 - t6 * v_5 - 0.5*(t6*t6)*a_5;
	double p_4 = p_5 - t5 * v_4 - 0.5*(t5*t5)*a_4 - (1.0 / 6.0)*(t5*t5*t5)*m_jerk_min;

	t4 = (p_4 - p_3) / v_3;


	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case6(t);
	}

	return true;

}
bool Trajectory_without_drag::calculate_trajectory_case6(double t[]) {
	
	double t1, t2, t3, t4, t5, t6, t7;
	double t1_sqrt;
	double v6;
	double p6;
	double t1_tmp_tmp_tmp;
	double t1_tmp_tmp;
	double b_t1_tmp_tmp;
	double c_t1_tmp_tmp;
	double t1_tmp;
	double t1_tmp2_tmp;
	double b_t1_tmp2_tmp;
	double c_t1_tmp2_tmp;
	double t1_tmp2_tmp_tmp;
	double b_t1_tmp2_tmp_tmp;
	double d_t1_tmp2_tmp;
	double e_t1_tmp2_tmp;
	double f_t1_tmp2_tmp;
	double c_t1_tmp2_tmp_tmp;
	double t1_sqrt3;
	double g_t1_tmp2_tmp;
	double h_t1_tmp2_tmp;
	double i_t1_tmp2_tmp;
	double t1_tmp2;
	double t1_sqrt2;
	t2 = 0.0;
	t4 = 0.0;
	t5 = m_acc_min / m_jerk_min;
	t7 = (m_acc_end - m_acc_min) / m_jerk_max;
	t1_sqrt = 0.5 * (t7 * t7);
	v6 = (m_vel_end - t7 * m_acc_min) - t1_sqrt * m_jerk_max;
	p6 = ((m_p_end - t7 * v6) - t1_sqrt * m_acc_min) - 0.16666666666666666 *
		pow(t7, 3.0) * m_jerk_max;
	t1_tmp_tmp_tmp = m_jerk_min * m_jerk_min;
	t1_sqrt = 9.0 * m_jerk_min * m_jerk_max - 9.0 * t1_tmp_tmp_tmp;
	t1_tmp_tmp = m_acc_start * m_acc_start;
	b_t1_tmp_tmp = t5 * t5;
	c_t1_tmp_tmp = 4.0 * m_jerk_min * m_jerk_max;
	t1_tmp = ((c_t1_tmp_tmp * b_t1_tmp_tmp / t1_sqrt - 4.0 * t1_tmp_tmp / t1_sqrt)
		- 8.0 * m_jerk_max * v6 / t1_sqrt) + 8.0 * m_jerk_max * m_vel_start /
		t1_sqrt;
	if (t1_tmp < 0.0) {
		t1 = -1.0;
		t2 = -1.0;
		t3 = -1.0;
		t4 = -1.0;
		t5 = -1.0;
		t6 = -1.0;
		t7 = -1.0;
		return calculate_trajectory_case7(t);
	}
	else {
		t1_sqrt = std::sqrt(t1_tmp);
		t1_tmp2_tmp = m_jerk_max * m_jerk_max;
		b_t1_tmp2_tmp = 8.0 * m_acc_min * t1_tmp2_tmp;
		t1_tmp = 4.0 * m_acc_min * m_jerk_min;
		c_t1_tmp2_tmp = 8.0 * m_acc_min * m_jerk_min;
		t1_tmp2_tmp_tmp = 4.0 * t1_tmp2_tmp;
		b_t1_tmp2_tmp_tmp = pow(m_acc_start, 3.0);
		d_t1_tmp2_tmp = ((((t1_tmp2_tmp_tmp * (v6 * v6) - 8.0 * m_acc_min *
			b_t1_tmp2_tmp_tmp / 3.0) + b_t1_tmp2_tmp * p6) -
			b_t1_tmp2_tmp * m_p_start) + 8.0 * m_acc_min * m_acc_start *
			m_jerk_max * m_vel_start) - b_t1_tmp2_tmp * t5 * v6;
		e_t1_tmp2_tmp = 4.0 * m_acc_min * t1_tmp_tmp * m_jerk_min * t1_sqrt;
		f_t1_tmp2_tmp = 2.0 * m_acc_min * t1_tmp_tmp * m_jerk_max * t1_sqrt;
		c_t1_tmp2_tmp_tmp = pow(t5, 3.0);
		t1_sqrt3 = c_t1_tmp2_tmp * t1_tmp2_tmp * c_t1_tmp2_tmp_tmp / 3.0;
		b_t1_tmp2_tmp = b_t1_tmp2_tmp * v6 * t1_sqrt;
		g_t1_tmp2_tmp = 4.0 * m_acc_min * t1_tmp2_tmp * m_vel_start * t1_sqrt;
		h_t1_tmp2_tmp = t1_tmp * t1_tmp2_tmp * b_t1_tmp_tmp * t1_sqrt;
		i_t1_tmp2_tmp = 2.0 * m_acc_min * t1_tmp_tmp_tmp * m_jerk_max * b_t1_tmp_tmp *
			t1_sqrt;
		t1_tmp = t1_tmp * m_jerk_max * v6 * t1_sqrt;
		c_t1_tmp2_tmp = c_t1_tmp2_tmp * m_jerk_max * m_vel_start * t1_sqrt;
		t1_tmp2 = ((((((((d_t1_tmp2_tmp - e_t1_tmp2_tmp) + f_t1_tmp2_tmp) + t1_sqrt3)
			- b_t1_tmp2_tmp) - g_t1_tmp2_tmp) + h_t1_tmp2_tmp) -
			i_t1_tmp2_tmp) + t1_tmp) + c_t1_tmp2_tmp;
		if (t1_tmp2 < 0.0) {
			t1 = -1.0;
			t2 = -1.0;
			t3 = -1.0;
			t4 = -1.0;
			t5 = -1.0;
			t6 = -1.0;
			t7 = -1.0;
			return calculate_trajectory_case7(t);
		}
		else {
			t1_sqrt2 = std::sqrt(t1_tmp2);
			t1_sqrt = ((((((((d_t1_tmp2_tmp + e_t1_tmp2_tmp) - f_t1_tmp2_tmp) +
				t1_sqrt3) + b_t1_tmp2_tmp) + g_t1_tmp2_tmp) -
				h_t1_tmp2_tmp) + i_t1_tmp2_tmp) - t1_tmp) - c_t1_tmp2_tmp;
			if (t1_sqrt < 0.0) {
				t1 = -1.0;
				t2 = -1.0;
				t3 = -1.0;
				t4 = -1.0;
				t5 = -1.0;
				t6 = -1.0;
				t7 = -1.0;
				return calculate_trajectory_case7(t);
			}
			else {
				t1_sqrt3 = std::sqrt(t1_sqrt);
				c_t1_tmp2_tmp = 2.0 * m_jerk_max * v6;
				f_t1_tmp2_tmp = c_t1_tmp2_tmp - t1_sqrt2;
				t1_tmp = 2.0 * m_acc_start * m_jerk_min;
				d_t1_tmp2_tmp = 6.0 * m_jerk_min * m_jerk_max;
				b_t1_tmp2_tmp = ((((((2.0 * t1_tmp_tmp * m_jerk_min - t1_tmp_tmp *
					m_jerk_max) + t1_tmp2_tmp_tmp * v6) + 2.0 *
					t1_tmp2_tmp * m_vel_start) - 2.0 * m_jerk_min *
					t1_tmp2_tmp * b_t1_tmp_tmp) + t1_tmp_tmp_tmp *
					m_jerk_max * b_t1_tmp_tmp) - 2.0 * m_jerk_min * m_jerk_max *
					v6) - c_t1_tmp_tmp * m_vel_start;
				t1_sqrt = 4.0 * m_acc_min * m_jerk_max * b_t1_tmp2_tmp;
				t1_tmp = ((((((((((m_acc_start * t1_tmp_tmp_tmp * b_t1_tmp_tmp -
					b_t1_tmp2_tmp_tmp) + 2.0 * t1_tmp_tmp_tmp * m_jerk_max *
					c_t1_tmp2_tmp_tmp) - t1_tmp * v6) + 4.0 * m_acc_start *
					m_jerk_max * v6) + t1_tmp * m_vel_start) + 2.0 * m_acc_start *
					m_jerk_max * m_vel_start) + d_t1_tmp2_tmp * p6) -
					d_t1_tmp2_tmp * m_p_start) - d_t1_tmp2_tmp * t5 * v6) -
					t1_tmp * m_jerk_max * b_t1_tmp_tmp) / b_t1_tmp2_tmp;
				d_t1_tmp2_tmp = 3.0 * m_jerk_min * v6;
				e_t1_tmp2_tmp = m_acc_min * b_t1_tmp2_tmp;
				t1 = (3.0 * m_jerk_min * (f_t1_tmp2_tmp * f_t1_tmp2_tmp) / t1_sqrt -
					t1_tmp) - d_t1_tmp2_tmp * f_t1_tmp2_tmp / e_t1_tmp2_tmp;

				// t1_1
				if (t1 < 0.0) {
					t1_tmp2 = c_t1_tmp2_tmp - t1_sqrt3;
					t1 = (3.0 * m_jerk_min * (t1_tmp2 * t1_tmp2) / t1_sqrt - t1_tmp) -
						d_t1_tmp2_tmp * t1_tmp2 / e_t1_tmp2_tmp;
				}

				if (t1 < 0.0) {
					t1_tmp2 = c_t1_tmp2_tmp + t1_sqrt2;
					t1 = (3.0 * m_jerk_min * (t1_tmp2 * t1_tmp2) / t1_sqrt - t1_tmp) -
						d_t1_tmp2_tmp * t1_tmp2 / e_t1_tmp2_tmp;
				}

				if (t1 < 0.0) {
					t1_tmp2 = c_t1_tmp2_tmp + t1_sqrt3;
					t1 = (3.0 * m_jerk_min * (t1_tmp2 * t1_tmp2) / t1_sqrt - t1_tmp) -
						d_t1_tmp2_tmp * t1_tmp2 / e_t1_tmp2_tmp;
				}

				// t3_1
				t1_tmp = 3.0 * t1_tmp2_tmp;
				c_t1_tmp2_tmp = 4.0 * m_acc_min * b_t1_tmp2_tmp;
				t1_tmp = 2.0 * (((((t1_tmp * p6 - t1_tmp * m_p_start) - b_t1_tmp2_tmp_tmp)
					- t1_tmp * t5 * v6) + m_jerk_min * t1_tmp2_tmp *
					c_t1_tmp2_tmp_tmp) + 3.0 * m_acc_start * m_jerk_max *
					m_vel_start) / b_t1_tmp2_tmp;
				t1_tmp2 = 3.0 * m_jerk_max * v6;
				t3 = (t1_tmp - 3.0 * (f_t1_tmp2_tmp * f_t1_tmp2_tmp) / c_t1_tmp2_tmp) +
					t1_tmp2 * f_t1_tmp2_tmp / e_t1_tmp2_tmp;
				if (t3 < 0.0) {
					t1_sqrt = 2.0 * m_jerk_max * v6 - t1_sqrt3;
					t3 = (t1_tmp - 3.0 * (t1_sqrt * t1_sqrt) / c_t1_tmp2_tmp) + t1_tmp2 *
						(2.0 * m_jerk_max * v6 - t1_sqrt3) / e_t1_tmp2_tmp;
				}

				if (t3 < 0.0) {
					t1_sqrt = 2.0 * m_jerk_max * v6 + t1_sqrt2;
					t3 = (t1_tmp - 3.0 * (t1_sqrt * t1_sqrt) / c_t1_tmp2_tmp) + t1_tmp2 *
						(2.0 * m_jerk_max * v6 + t1_sqrt2) / e_t1_tmp2_tmp;
				}

				if (t3 < 0.0) {
					t1_sqrt = 2.0 * m_jerk_max * v6 + t1_sqrt3;
					t3 = (t1_tmp - 3.0 * (t1_sqrt * t1_sqrt) / c_t1_tmp2_tmp) + t1_tmp2 *
						(2.0 * m_jerk_max * v6 + t1_sqrt3) / e_t1_tmp2_tmp;
				}

				// t6_1
				t1_sqrt = 2.0 * m_acc_min * m_jerk_max;
				t6 = -f_t1_tmp2_tmp / t1_sqrt;
				if (t6 < 0.0) {
					t6 = -(2.0 * m_jerk_max * v6 - t1_sqrt3) / t1_sqrt;
				}

				if (t6 < 0.0) {
					t6 = -(2.0 * m_jerk_max * v6 + t1_sqrt2) / t1_sqrt;
				}

				if (t6 < 0.0) {
					t6 = -(2.0 * m_jerk_max * v6 + t1_sqrt3) / t1_sqrt;
				}
			}
		}
	}
	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return calculate_trajectory_case7(t);
	}
	
	return true;
}
bool Trajectory_without_drag::calculate_trajectory_case7(double t[]) {
	

	double t1, t2, t3, t4, t5, t6, t7;

	t1 = (m_acc_max - m_acc_start)/m_jerk_max;
	t3 = -m_acc_max/m_jerk_min;
	t6 = 0.0;
	t4 = 0.0;


	double a1 = (m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max)*3.0+(m_acc_max*m_acc_max*m_acc_max*m_acc_max)*(m_jerk_min*m_jerk_min)-(m_acc_max*m_acc_max*m_acc_max*m_acc_max)*(m_jerk_max*m_jerk_max)-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*3.0-(m_acc_end*m_acc_end*m_acc_end)*m_acc_max*(m_jerk_max*m_jerk_max)*8.0+m_acc_max*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*8.0+(m_acc_end*m_acc_end)*(m_acc_max*m_acc_max)*(m_jerk_max*m_jerk_max)*6.0-(m_acc_max*m_acc_max)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*6.0+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-m_acc_max*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1+m_acc_max*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*1.2E1-(m_acc_max*m_acc_max)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*1.2E1+(m_acc_max*m_acc_max)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.2E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.2E1+m_acc_end*m_acc_max*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*2.4E1-m_acc_max*m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1;
	double a2 = (m_jerk_max*m_jerk_max)*(m_acc_end-m_acc_max)*(m_jerk_min-m_jerk_max)*(m_acc_end*m_acc_max+m_jerk_min*m_vel_end*2.0-m_acc_end*m_acc_end)*-1.2E1;
	double a3 = (m_jerk_max*m_jerk_max)*(m_jerk_min-m_jerk_max)*((m_acc_end*m_acc_end)*m_jerk_min*2.0-(m_acc_end*m_acc_end)*m_jerk_max*3.0-(m_acc_max*m_acc_max)*m_jerk_max-m_acc_end*m_acc_max*m_jerk_min*2.0+m_acc_end*m_acc_max*m_jerk_max*4.0+m_jerk_min*m_jerk_max*m_vel_end*2.0)*6.0;
	double a4 = (m_jerk_max*m_jerk_max*m_jerk_max)*(m_jerk_min-m_jerk_max)*(m_acc_end*m_jerk_min*3.0-m_acc_end*m_jerk_max*3.0-m_acc_max*m_jerk_min+m_acc_max*m_jerk_max*2.0)*-4.0;
	double a5 = (m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*pow(m_jerk_min-m_jerk_max,2.0)*3.0;

	t7 = fifth_2deriv(a1, a2, a3, a4, a5);
	t5 = (m_acc_end-m_jerk_max*t7)/m_jerk_min;
	t2 = ((m_jerk_max*m_jerk_max)*(t7*t7)+m_jerk_max*m_vel_end*2.0-m_jerk_max*m_vel_start*2.0-m_acc_max*m_acc_max+m_acc_start*m_acc_start-m_acc_end*m_jerk_max*t7*2.0)/(m_acc_max*m_jerk_max*2.0)-((m_jerk_max*m_jerk_max)*(t7*t7)+m_acc_end*m_acc_end-m_acc_max*m_acc_max-m_acc_end*m_jerk_max*t7*2.0)/(m_acc_max*m_jerk_min*2.0);



	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		/*if (t[i] < 0)
			return calculate_trajectory_case8(t);*/
	}
	
	return true;
}
bool Trajectory_without_drag::calculate_trajectory_case8(double t[]) {
	
	double t1, t2, t3, t4, t5, t6, t7;

	t2 = 0;
	t4 = 0;
	t6 = 0;



	double a1  = (m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_max-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*6.0-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*8.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*3.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_max*3.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_max*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_max*6.0-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_p_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_p_start*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max)*m_vel_start*1.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*3.2E1-(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end*m_vel_end)*8.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*8.0+(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start*m_vel_start)*3.2E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_vel_start*m_vel_start)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*3.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*3.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_vel_end*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*m_vel_start*1.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*1.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*2.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*4.8E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*m_vel_end*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*m_vel_start*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_end*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_start*6.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_vel_end*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max)*m_vel_start*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*2.4E1-(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_end*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_start*2.4E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*7.2E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*2.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*4.8E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*m_jerk_max*m_vel_start*2.4E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*7.2E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*m_vel_start*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*m_vel_start*7.2E1+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*7.2E1-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*7.2E1;
	double a2   = (m_jerk_min-m_jerk_max)*((m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_max-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*4.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*2.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_max*4.0-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_max*4.0-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*1.2E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*1.2E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_vel_end*4.0-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.0+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max)*m_vel_start*8.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*2.0-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_max-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*m_vel_end*4.0-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*m_vel_start*4.0-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_end*4.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.2E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_vel_end*4.0+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*m_vel_start*4.0-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_p_end*1.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_p_start*1.2E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0+m_acc_end*(m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*4.0+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.2E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.2E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_end*8.0-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*m_vel_start*8.0+(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*m_jerk_max*m_vel_start*1.2E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.6E1+m_acc_end*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*8.0+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*1.2E1)*6.0;
	double a3   = (m_jerk_min-m_jerk_max)*((m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max)*5.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*4.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*8.0-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.0+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*4.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*4.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*2.0-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.2E1-m_acc_end*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*2.4E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*3.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*3.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.6E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.4E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*1.2E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*1.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*8.0+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*8.0-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*2.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*8.0+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*8.0-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*4.0-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*8.0+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*8.0+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*3.6E1-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.2E1)*-3.0;
	double a4   = m_jerk_max*pow(m_jerk_min-m_jerk_max,2.0)*(m_jerk_min-m_jerk_max*2.0)*(m_acc_end*(m_acc_start*m_acc_start)*3.0+(m_jerk_max*m_jerk_max)*m_p_end*3.0-(m_jerk_max*m_jerk_max)*m_p_start*3.0-(m_acc_end*m_acc_end*m_acc_end)*2.0-m_acc_start*m_acc_start*m_acc_start+m_acc_end*m_jerk_max*m_vel_end*3.0-m_acc_end*m_jerk_max*m_vel_start*6.0+m_acc_start*m_jerk_max*m_vel_start*3.0)*-4.0;
	double a5   = (m_jerk_max*m_jerk_max)*pow(m_jerk_min-m_jerk_max,2.0)*(m_jerk_min-m_jerk_max*2.0)*(m_jerk_max*m_vel_end*2.0-m_jerk_max*m_vel_start*2.0-m_acc_end*m_acc_end+m_acc_start*m_acc_start)*3.0;

	t7 = fifth_2deriv(a1, a2, a3, a4, a5);

	t1 = -(-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max)+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*2.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*2.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max)*2.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*8.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.0E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_end*4.0-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*1.6E1-m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*1.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*8.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_start*4.0-m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*8.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*3.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*2.0+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*8.0+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*1.6E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*4.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*4.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*3.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*5.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max)*4.0-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*8.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*2.4E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*2.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*m_jerk_max*4.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_end*8.0+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_start*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_vel_start*8.0+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*2.4E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*t7*1.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*t7*8.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.5E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*2.0E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.8E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*1.08E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*5.4E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*1.08E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*1.08E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*5.4E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*6.6E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*1.08E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*6.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*7.2E1+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*4.0E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*1.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*9.6E1+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*3.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*9.6E1+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*3.2E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start)*3.2E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.8E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*4.0E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*1.2E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*1.6E1-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start)*4.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*4.7E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(t7*t7)*1.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.6E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(t7*t7)*2.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.2E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*3.2E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*6.0E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*7.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*5.6E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*1.6E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*2.4E1-m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*4.0E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.92E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.2E2-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*2.4E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.92E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*4.8E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.2E2+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*2.4E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*4.8E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*1.08E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*3.6E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*1.08E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*1.2E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*8.0E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*8.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*9.6E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*3.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*8.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*m_jerk_max+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*8.0-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*2.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.0-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.12E2+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*6.8E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*5.3E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*8.4E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*2.4E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*8.4E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.1E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.0E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*m_jerk_max*1.9E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.2E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.8E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*t7*4.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*7.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*6.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*7.2E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*7.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*6.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.2E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*3.6E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*7.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*3.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*9.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*8.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*1.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*7.2E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*7.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*3.8E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*8.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*1.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*4.0+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*7.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*2.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*1.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7*2.0E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*m_jerk_max*t7*1.6E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*2.4E1+m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*4.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.2E2+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*2.4E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.2E2-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*7.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*7.2E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*4.8E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*4.8E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*4.8E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*2.4E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.56E2-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.08E2+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*1.92E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*6.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.02E2+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*6.2E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*9.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.68E2+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*9.6E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*2.52E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*7.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.32E2-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*6.8E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.0E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*4.0-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*3.6E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*6.0E1+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*1.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*5.0E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*1.0E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*6.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*6.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*4.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*3.0E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7*3.2E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*t7*2.0E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*2.16E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*2.16E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*3.2E1+m_acc_end*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*9.6E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*8.0E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*4.8E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.2E2+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*4.8E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.2E2-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*4.8E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*4.8E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.2E2+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*4.8E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.2E2-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*1.44E2+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*1.44E2+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*1.92E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*3.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*1.6E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*7.2E1-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.16E2-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*5.6E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*4.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*t7*3.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.58E2+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*t7*1.6E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*4.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*t7*3.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*6.0-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*2.16E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*7.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*1.16E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*4.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7*t7)*1.2E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*4.8E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*9.6E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*2.4E1+m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.8E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*1.68E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*1.68E2+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*1.68E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*1.68E2-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*9.6E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*9.6E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*9.6E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*1.68E2-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*1.2E2-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*1.68E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*7.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_end*3.6E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*7.2E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*7.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_p_start*3.6E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*7.2E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*8.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_end*1.6E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*4.4E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start)*1.2E1+(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*6.0E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(m_vel_start*m_vel_start)*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*5.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*m_vel_start*4.0E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*8.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7*4.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*t7*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*t7*4.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*t7*7.9E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*t7*6.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*6.3E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(t7*t7)*1.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.4E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*4.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*6.0E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*7.2E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*3.6E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*3.6E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.92E2-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.2E2+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.92E2-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.2E2+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*1.32E2-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*t7*3.6E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*1.32E2+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*t7*3.6E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*7.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*t7*2.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*4.8E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*7.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*t7*2.4E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*4.8E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*3.36E2-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.68E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*8.8E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.88E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.2E2+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*4.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*4.8E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*1.0E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*1.56E2+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*1.72E2-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.68E2-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*1.28E2+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*6.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*7.8E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*2.52E2-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.2E2+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*3.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*3.6E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*9.6E1-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*8.0E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*6.0E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*6.0E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*6.0E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*1.2E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*1.2E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*6.0E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*3.12E2+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*3.36E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*1.44E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*1.92E2+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*1.68E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*4.8E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*1.44E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*1.44E2-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*1.2E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*9.6E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*9.6E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*9.6E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.4E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*2.4E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*1.44E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*6.3E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.6E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(t7*t7)*2.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.1E2+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*(t7*t7)*3.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*4.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*9.6E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*3.0E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*6.0+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*3.0E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*6.0+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*3.0E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*6.0-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*3.0E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*6.0-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.2E2-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.2E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*1.12E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*9.6E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*t7*3.6E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.04E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.2E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*1.0E2-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.68E2+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*1.32E2+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.52E2+(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*7.2E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.44E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*1.32E2+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.12E2-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.44E2+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*6.0E1-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*1.2E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.44E2+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.44E2+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.44E2-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.44E2-m_acc_end*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*9.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*3.36E2-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.92E2-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*3.36E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.92E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.44E2+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.44E2-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*2.4E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.22E2-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*6.6E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*9.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.52E2+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*7.2E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*3.0E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*6.0+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*7.2E1-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*2.16E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*7.2E1+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*9.6E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*4.8E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*2.4E1-(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*m_vel_start*3.6E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*1.12E2-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*t7*3.2E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.92E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.6E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*t7*2.0E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.44E2-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*t7*3.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*3.36E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*9.6E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.2E2+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*3.36E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*2.4E1-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*9.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.2E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*2.4E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*2.16E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*2.16E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*7.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.92E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.08E2-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.68E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*9.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*6.0E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*t7*1.2E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*6.0E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*t7*1.2E1-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.36E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*9.6E1+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.36E2)/((m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*2.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*4.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*1.6E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*2.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*3.0+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*9.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*3.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*9.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*6.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*3.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*1.6E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*2.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*1.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*7.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*7.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*7.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*7.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.8E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*9.6E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*1.44E2+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.8E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*6.0E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*2.16E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*9.6E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.8E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*2.16E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*7.2E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*2.16E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.16E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.16E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*7.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*7.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.44E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.44E2+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.44E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.44E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*4.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*4.8E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.44E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*7.2E1);
	t3 = -((m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*4.0-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*2.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)*4.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*2.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*2.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.0E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*1.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*8.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*3.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*2.0-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*8.0-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*1.6E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*8.0+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*8.0-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*3.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*8.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*4.0+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*8.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*2.4E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*8.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*4.0-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*2.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*2.0E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*2.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*1.8E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*1.9E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*4.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*2.4E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.0E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.5E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*2.0E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.8E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*1.08E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*5.4E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*1.08E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*5.4E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*1.2E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*1.2E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*7.2E1-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*4.0E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*1.6E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*9.6E1-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*3.2E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*8.0-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*2.8E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*2.4E1+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*1.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*t7*1.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*4.7E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.2E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*3.2E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*6.0E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*7.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*5.6E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*1.6E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*2.4E1+m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*4.8E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*4.0E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.92E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.2E2+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.92E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*4.8E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.2E2-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*2.4E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*4.8E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*3.6E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*1.2E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*8.0E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*8.4E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*9.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*3.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*8.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*2.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*1.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.6E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*8.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*6.4E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.0E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*6.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*5.3E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*1.6E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*2.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*2.1E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*6.0+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*8.4E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7*t7)*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*8.4E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7*t7)*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*t7*1.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*7.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*6.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*7.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*6.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*6.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*6.0-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*9.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*8.0-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*7.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*3.8E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.6E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.8E1-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.4E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.8E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*2.0E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*2.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.0+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*4.0E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_end*m_vel_end)*1.2E2-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(m_vel_start*m_vel_start)*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_end*m_vel_end)*1.2E2+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(m_vel_start*m_vel_start)*2.4E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*t7*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*t7*7.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*t7*2.4E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*t7*7.2E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*4.8E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*2.4E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.56E2+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.08E2-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*1.92E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*6.0E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.02E2-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*6.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*9.6E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*1.68E2-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*9.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*2.52E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*7.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.32E2+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*6.8E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.0E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*4.0+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*3.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*6.0E1-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*1.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*5.0E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*1.0E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.8E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.6E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*3.2E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*t7*4.0+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*2.16E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*4.8E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*3.2E1-m_acc_end*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*8.0E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*4.8E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.2E2-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*4.8E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.2E2+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*4.8E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*4.8E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*4.8E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*1.44E2-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*1.44E2-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*4.8E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*1.92E2-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*3.2E1-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.2E2+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*8.0E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.6E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.58E2-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*4.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*6.0+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*2.16E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*7.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.0E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*2.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.16E2-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*3.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*4.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*1.2E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*4.8E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*9.6E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*2.4E1-m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.8E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*1.68E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*1.68E2-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*1.68E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*1.68E2+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*9.6E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*1.68E2+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*1.2E2+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*1.68E2-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*7.2E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*7.2E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*2.4E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*8.4E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.0E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.6E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*8.0+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1-(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.0E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*6.0E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*8.0+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*3.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*5.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*t7*2.0E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*4.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*t7*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*t7*7.9E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*t7*8.0+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*6.3E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.4E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*4.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*6.0E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*7.2E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*3.6E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*3.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*8.4E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.4E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*4.8E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*8.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*1.32E2+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*3.6E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*1.32E2-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*3.6E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*7.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*2.4E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*4.8E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*4.8E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*7.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*2.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*4.8E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*3.36E2+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.68E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*8.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.6E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*1.68E2-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.92E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*1.68E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*1.2E2-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.4E2+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.4E2+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*4.8E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*1.0E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*3.2E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*(t7*t7)*1.56E2-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*1.72E2+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.68E2+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.28E2-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*1.6E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*t7*6.0E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*7.8E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*3.6E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*2.52E2+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.2E2-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*3.2E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7*t7)*8.4E1-(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*3.6E1-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*(t7*t7)*9.6E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*8.0E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7*t7)*8.4E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*6.0E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*(t7*t7)*1.2E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*6.0E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*6.0E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*(t7*t7)*1.2E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*(t7*t7)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*6.0E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*(t7*t7)*1.2E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*3.12E2-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*3.36E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*t7*1.44E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*t7*1.92E2-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*1.68E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7*t7)*4.8E1+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*1.44E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*1.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*3.6E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*3.6E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*9.6E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*9.6E1-m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*9.6E1+m_acc_end*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*9.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*8.4E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*2.8E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*4.4E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*2.4E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*7.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*2.0E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*8.0E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*t7*1.44E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*t7*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*t7*4.4E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*t7*6.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*6.3E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(t7*t7)*1.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.6E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7)*1.1E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(t7*t7*t7)*4.2E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*m_vel_start*9.6E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*m_vel_start*9.6E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*3.0E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*6.0-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*3.0E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*6.0-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*3.0E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*(t7*t7)*6.0+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*3.0E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*(t7*t7)*6.0+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*3.36E2-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*8.4E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*3.6E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*1.12E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*3.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.04E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.2E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*2.4E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.0E2-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*2.0E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.68E2-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.32E2+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*t7*3.6E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.52E2-(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*4.8E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*7.2E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.44E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*1.32E2-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.12E2+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.44E2-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*6.0E1+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*1.2E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.44E2+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.44E2+m_acc_end*m_acc_start*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*9.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*3.36E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.92E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*3.36E2-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.92E2+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.44E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.44E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*t7*3.6E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*2.22E2+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(t7*t7)*6.6E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*2.52E2-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*7.2E1+(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*3.0E1-(m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*(t7*t7)*6.0+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.16E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*9.6E1-m_acc_end*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.92E2+(m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.2E2+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*t7*1.12E2+m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*t7*1.92E2-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.6E2+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*t7*1.44E2+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*3.36E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*t7*9.6E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*1.2E2-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*3.36E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*2.4E1+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*t7*9.6E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*1.2E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*2.4E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*2.16E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*t7*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*2.16E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*t7*7.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.92E2-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*4.08E2+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*1.68E2+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*(t7*t7)*9.6E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*6.0E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*t7*1.2E1-m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*6.0E1+m_acc_end*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*t7*1.2E1+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.36E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*9.6E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*t7*3.36E2)/(m_jerk_min*((m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max)*2.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*4.0+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*1.6E1-(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max)+(m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*2.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*6.0+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*7.2E1+m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*4.8E1+(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*2.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*8.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max)*3.0+(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*2.4E1-(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1+(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_end*m_p_end)*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*1.08E2+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_p_start*m_p_start)*3.6E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*9.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end*m_vel_end)*3.2E1-(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*9.6E1+(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start*m_vel_start)*6.4E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*3.0-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*1.6E1+(m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*2.4E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max)*1.2E1-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*m_jerk_max*1.2E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*1.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*7.2E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*7.2E1-(m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*7.2E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*7.2E1+(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_p_start*2.4E1-(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*1.2E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*2.4E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.8E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*9.6E1-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*9.6E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*1.44E2+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*3.6E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*1.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1+(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*4.8E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*4.8E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*3.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*(m_vel_start*m_vel_start)*6.0E1+(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1-(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*4.8E1+(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*4.8E1+(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0-(m_acc_end*m_acc_end*m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*6.0-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1-(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1+(m_acc_start*m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1+(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*2.16E2-(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_p_start*7.2E1-(m_acc_end*m_acc_end)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*(m_vel_start*m_vel_start)*2.4E1-m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*(m_vel_end*m_vel_end)*m_vel_start*9.6E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.8E1+m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_end*2.4E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*4.8E1+(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*2.4E1+m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*2.16E2-m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*7.2E1-m_acc_end*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*2.16E2+m_acc_end*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*7.2E1-m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*2.16E2+m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*7.2E1+m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*2.16E2-m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*7.2E1+(m_acc_end*m_acc_end)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1+(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.4E1-m_acc_end*(m_acc_start*m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*7.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*4.8E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_start*1.2E1-(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*m_jerk_max*m_vel_start*4.8E1-(m_acc_end*m_acc_end*m_acc_end)*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*7.2E1-(m_acc_end*m_acc_end)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1-(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1+(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*9.6E1-m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_end*1.44E2+m_acc_end*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_end*1.44E2+m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_end*m_vel_start*1.44E2-m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_p_start*m_vel_start*1.44E2+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*4.8E1+(m_acc_end*m_acc_end)*(m_acc_start*m_acc_start)*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_start*4.8E1-m_acc_end*m_acc_start*m_jerk_min*(m_jerk_max*m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*1.44E2+m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*2.16E2-m_acc_end*m_acc_start*(m_jerk_min*m_jerk_min*m_jerk_min)*(m_jerk_max*m_jerk_max)*m_vel_end*m_vel_start*7.2E1));
	t5 = (m_acc_end-m_jerk_max*t7)/m_jerk_min;

	t[0] = t1;
	t[1] = t2;
	t[2] = t3;
	t[3] = t4;
	t[4] = t5;
	t[5] = t6;
	t[6] = t7;

	for (unsigned int i = 0; i < num_time_segments; i++)
	{
		if (t[i] < 0)
			return false;
	}
	
	return true;
}


void Trajectory_without_drag::calculate_theta_for_two_points(char axis, std::vector<Point3D>& theta_refs)
{

	double sign = check_point();
	double t[7] = { 0,0,0,0,0,0,0 };

	if (calculate_trajectory_case7(t))
	{
		double time_in_milisec = 0.0;
		double acc_current = m_acc_start;


		std::cout << t[0] << std::endl << t[1] << std::endl << t[2] << std::endl << t[3] << std::endl << t[4] << std::endl << t[5] << std::endl << t[6] << std::endl;

		double acc_start = m_acc_start;
		int num_of_elem = static_cast<int>((t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6])*1000.0);
		int last_point_theta = theta_refs.size() - 1;
		theta_refs.reserve(num_of_elem + theta_refs.size());

		std::cout << num_of_elem << std::endl;


		for (int k = 0; k < num_time_segments; k++)
		{
			double jerk_current = 0.0;

			switch (k + 1)
			{
			case 1:
				jerk_current = m_jerk_max;
				break;

			case 2:
				jerk_current = 0.0;
				break;

			case 3:
				jerk_current = m_jerk_min;
				break;

			case 4:
				jerk_current = 0.0;
				break;

			case 5:
				jerk_current = m_jerk_min;
				break;

			case 6:
				jerk_current = 0.0;
				break;

			case 7:
				jerk_current = m_jerk_max;
				break;
			default:
				break;
			}

			time_in_milisec = t[k] * 1000.0;
			int num_of_time_slices = static_cast<int>(time_in_milisec);
			double time_in_seconds = 0.0;

			for (int k = 0; k < num_of_time_slices; k++)
			{
				Point3D theta_ref;



				acc_current = calculate_acceleration(time_in_seconds, acc_start, jerk_current);
				switch (static_cast<int>(axis)) {

				case 'x':
					theta_ref.x = sign * (acc_current / 9.78) + (theta_refs[last_point_theta]).x;
					theta_ref.y = 0.0;
					theta_ref.z = 0.0;
					break;
				case 'y':

					theta_ref.x = 0.0;
					theta_ref.y = sign * (acc_current / 9.78) + (theta_refs[last_point_theta]).y;
					theta_ref.z = 0.0;
					break;
				case 'z':

					theta_ref.x = 0.0;
					theta_ref.y = 0.0;
					theta_ref.z = sign * (acc_current / 9.78) + (theta_refs[last_point_theta]).y;
					break;
				}

				time_in_seconds += 0.001;
				theta_refs.push_back(theta_ref);

			}
			acc_start = acc_current;

		}




	}

}


std::vector<Point3D> Trajectory_without_drag::calculate_theta_ref(std::vector<Point3D> points)
{

	int number_of_spaces = points.size() - 1;
	double sign = 1.0;
	std::vector<Point3D> theta_refs;

	Point3D first;
	theta_refs.push_back(first);

	for (int i = 0; i < number_of_spaces; i++)
	{

		m_p_start = (points[i]).x;
		m_p_end = (points[i + 1]).x;

		if (abs(m_p_end - m_p_start) > 1e-5)
			calculate_theta_for_two_points('x', theta_refs);

		m_p_start = (points[i]).y;
		m_p_end = (points[i + 1]).y;

		if (abs(m_p_end - m_p_start) > 1e-5)
			calculate_theta_for_two_points('y', theta_refs);

		m_p_start = (points[i]).z;
		m_p_end = (points[i + 1]).z;

		if (abs(m_p_end - m_p_start) > 1e-5)
			calculate_theta_for_two_points('z', theta_refs);



	}

	return theta_refs;
}