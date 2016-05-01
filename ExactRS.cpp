#include "ExactRS.hpp"
#include "universal_error.hpp"
#include <cmath>
#include <algorithm>

namespace
{
	double CalcFrarefraction(Primitive const& cell, double p,double gamma)
	{
		double cs = sqrt(gamma*cell.pressure/ cell.density);
		return 2 * cs*(pow(p / cell.pressure, (gamma - 1) / (2 * gamma)) - 1) / (gamma - 1);
	}

	double CalcFshock(Primitive const& cell, double p,double gamma)
	{
		double A = 2 / ((gamma + 1) *cell.density);
		double B = (gamma - 1)*cell.pressure / (gamma + 1);
		return (p - cell.pressure)*sqrt(A / (p + B));
	}

	double dCalcFrarefraction(Primitive const& cell, double p, double gamma)
	{
		return (p + gamma*p + 3 * gamma*cell.pressure - cell.pressure)*pow(p + gamma*p - cell.pressure + 
			gamma*cell.pressure, -1.5) / sqrt(2 * cell.density);
	}

	double dCalcFshock(Primitive const& cell, double p, double gamma)
	{
		double cs = sqrt(gamma*cell.pressure/ cell.density);
		return pow(cell.pressure / p, (gamma + 1) / (2 * gamma)) / (cell.density*cs);
	}

	RSsolution GetFirstGuess(Primitive const & left, Primitive const & right,double gamma)
	{
		RSsolution res;
		double csl = sqrt(gamma*left.pressure / left.density);
		double csr = sqrt(gamma*right.pressure / right.density);
		res.pressure = pow((csl+csr-0.5*(gamma-1)*(right.velocity-left.velocity))/(csl*pow(
			left.pressure,(-gamma+1)/(2*gamma))+ csr*pow(
				right.pressure, (-gamma + 1) / (2 * gamma))),2*gamma/(gamma-1));
		return res;
	}

	double GetValue(Primitive const & left, Primitive const & right,double p,double gamma)
	{
		double res = right.velocity-left.velocity;
		if (p > left.pressure)
			res += CalcFshock(left, p, gamma);
		else
			res += CalcFrarefraction(left, p,gamma);
		if (p > right.pressure)
			res += CalcFshock(right, p, gamma);
		else
			res += CalcFrarefraction(right, p, gamma);
		return res;
	}

	double GetdValue(Primitive const & left, Primitive const & right, double p, double gamma)
	{
		double res = 0;
		if (p > left.pressure)
			res += dCalcFshock(left, p, gamma);
		else
			res += dCalcFrarefraction(left, p, gamma);
		if (p > right.pressure)
			res += dCalcFshock(right, p, gamma);
		else
			res += dCalcFrarefraction(right, p, gamma);
		return res;
	}
}

ExactRS::ExactRS(double gamma):gamma_(gamma)
{}

ExactRS::~ExactRS()
{}

RSsolution ExactRS::Solve(Primitive const & left, Primitive const & right)const
{
	const double eps = 1e-7;
	// Is there a vaccum?
	double dv = right.velocity - left.velocity;
	double soundspeeds = 2 * (sqrt(gamma_*left.pressure / left.density) + sqrt(gamma_*right.pressure / right.density)) / (gamma_ - 1);
	if (dv > soundspeeds)
	{
		RSsolution res;
		res.pressure = 0;
		res.velocity = 0;
		return res;
	}
	RSsolution res = GetFirstGuess(left, right,gamma_);
	res.velocity = 0;
	if (res.pressure < 0)
	{
		res.pressure = 0;
		return res;
	}
	double value = GetValue(left, right, res.pressure, gamma_);
	//	double Cs = std::max(sqrt(gamma_*left.pressure / left.density), sqrt(gamma_*right.pressure / right.density));
	double dp = 0;
	double p = res.pressure;
	size_t counter = 0;
	do
	{
		p = res.pressure;
		dp = value / GetdValue(left, right, res.pressure, gamma_);
		res.pressure -= std::max(std::min(0.5*dp,res.pressure*0.1),-0.1*res.pressure);
		value = GetValue(left, right, res.pressure, gamma_);
		++counter;
		if (counter > 200)
		{
			UniversalError eo("Too many iterations in RS");
			eo.AddEntry("Left density", left.density);
			eo.AddEntry("Left pressure", left.pressure);
			eo.AddEntry("Left velocity", left.velocity);
			eo.AddEntry("Right density", right.density);
			eo.AddEntry("Right pressure", right.pressure);
			eo.AddEntry("Right velocity", right.velocity);
			throw eo;
		}
	} while ((fabs(dp) > eps*(p+res.pressure))&&(dv*eps>value));
	double fr = (res.pressure > right.pressure) ? CalcFshock(right, res.pressure, gamma_) : CalcFrarefraction(
		right, res.pressure, gamma_);
	double fl = (res.pressure > left.pressure) ? CalcFshock(left, res.pressure, gamma_) : CalcFrarefraction(
		left, res.pressure, gamma_);
	res.velocity = 0.5*(left.velocity + right.velocity) + 0.5*(fr-fl);
	return res;
}
