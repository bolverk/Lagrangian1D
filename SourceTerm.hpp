#ifndef SOURCETERM_HPP
#define SOURCETERM_HPP 1

#include "Primitive.hpp"
#include "Extensive.hpp"
#include <vector>

using namespace std;

class SourceTerm
{
public:
	virtual void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double time,
		vector<Extensive> &extensives,double dt)const=0;

	virtual ~SourceTerm();
};


class ZeroForce : public SourceTerm
{
  void CalcForce(vector<double> const& /*edges*/, vector<Primitive> const& /*cells*/, double /*time*/,
		 vector<Extensive> & /*extensives*/,double /*dt*/)const
	{
		return;
	}
};
#endif //SOURCETERM_HPP
