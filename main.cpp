#define _USE_MATH_DEFINES
#include "hdsim.hpp"
#include "hdf_util.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <boost/math/tools/roots.hpp>
#include <float.h>
#ifdef _MSC_VER // Checks if code is compiled by Micro$oft visual studio
#include <Windows.h>

unsigned int oldd = 0;
unsigned int fp_control_state = _controlfp_s(&oldd,_EM_INEXACT, _MCW_EM);
#endif // _MSC_VER

namespace
{
	string int2str(int n)
	{
		stringstream ss;
		ss << n;
		return ss.str();
	}

	template <typename T>
	T LinearInterpolation(const vector<T> &x, const vector<T> &y, T xi)
	{
		typename vector<T>::const_iterator it = upper_bound(x.begin(), x.end(), xi);
		assert(it != x.end() && it != x.begin() &&
			"X out of range in Linear Interpolation");
		if (*it == xi)
			return y[static_cast<size_t>(it - x.begin())];

		// Are we near the edge?
		if (it == x.end() - 1)
			return y[static_cast<size_t>(it - x.begin())] + (xi - *it)*
			(y[static_cast<size_t>(it - 1 - x.begin())] -
				y[static_cast<size_t>(it - x.begin())]) / (*(it - 1) - *it);
		else
			return y[static_cast<size_t>(it - x.begin())] + (xi - *it)*
			(y[static_cast<size_t>(it + 1 - x.begin())] -
				y[static_cast<size_t>(it - x.begin())]) / (*(it + 1) - *it);
	}

	bool missing_file_data(string const& fname)
	{
		std::cout << "Could not find file " << fname << std::endl;
		return false;
	}

	vector<double> read_vector(string const& fname)
	{
		double buf = 0;
		vector<double> res;
		ifstream f(fname.c_str());
		assert(f || missing_file_data(fname));
		while (f >> buf)
			res.push_back(buf);
		f.close();
		return res;
	}

	vector<Primitive> calc_init(vector<double> const& edges, double M, double R, double Nemden,double stargamma)
	{
		// Get lane emden result, G=1 R=Rsun M=Msun
		vector<double> xsi, theta, dtheta;
		if (stargamma < 1.6)
		{
			xsi = read_vector("xsin3.txt");
			theta = read_vector("thetan3.txt");
			dtheta = read_vector("dthetan3.txt");
		}
		else
		{
			xsi = read_vector("xsi.txt");
			theta = read_vector("theta.txt");
			dtheta = read_vector("dtheta.txt");
		}
		size_t n = xsi.size();
		double rhoc = -M*xsi[n - 1] / (4 * M_PI*R*R*R*dtheta[n - 1]);
		double K = M*M*pow(4 * M_PI, 1.0 / Nemden)*pow(-M*xsi[n - 1] / (dtheta[n - 1] * R*R*R), -(Nemden + 1.0) / Nemden) /
			((dtheta[n - 1] * dtheta[n - 1])*(1 + Nemden)*pow(R, 4));
		vector<Primitive> res(edges.size() - 1);
		for (size_t i = 0; i < res.size(); ++i)
		{
			double y = 0.5 * xsi[n - 1] * (edges[i + 1] + edges[i]) / R;
			double Theta = 0;
			bool outside = false;
			if (y > xsi[n - 1])
			{
				Theta = theta[n - 2] * 0.1;
				outside = true;
			}
			else
				Theta = LinearInterpolation(xsi, theta, y);
			res[i].density = rhoc*pow(Theta, Nemden);
			res[i].pressure = K*pow(res[i].density, (Nemden + 1) / Nemden);
			if (outside)
				res[i].pressure = res[i].density*0.1;
			res[i].velocity = 0;
		}
		return res;
	}

	vector<double> getedges(size_t N, double L)
	{
	  double dx = L / static_cast<double>(N);
		vector<double> res(N + 1);
		for (size_t i = 0; i <= N; ++i)
		  res[i] = static_cast<double>(i)*dx;
		return res;
	}

	class ParaboleAnomaly
	{
	private:
		double t_, Rp_, Mbh_;
	public:
		ParaboleAnomaly(double t, double Rp, double Mbh) :t_(t), Rp_(Rp), Mbh_(Mbh) {}
		double operator()(double f)
		{
			return t_ - sqrt(2 * pow(Rp_, 3) / Mbh_)*tan(f / 2)*(3 + pow(tan(f / 2), 2)) / 3;
		}
	};

	struct TerminationCondition
	{
		bool operator() (double min, double max)
		{
			return abs(min - max) <= 0.000001;
		}
	};

	double GetTrueAnomaly(double t, double Rp, double Mbh)
	{

		ParaboleAnomaly Eanom(t, Rp, Mbh);
		pair<double, double> res = boost::math::tools::bisect(Eanom, -M_PI, M_PI,
			TerminationCondition());
		return res.first;
	}

	class Gravity : public SourceTerm
	{
	private:
		vector<double> acc_;
		double rhoc_;
		double alpha_;
		double Mbh_;
		double Rp_;
		mutable double f_;
		bool selfgravity_;
	public:
		Gravity(double M, double R, double Mbh, double Rp,bool selfgravity,double stargamma, vector<double> const& edges) :
			acc_(vector<double>()),rhoc_(0),alpha_(0), Mbh_(Mbh), Rp_(Rp), f_(0),selfgravity_(selfgravity)
		{
			vector<double> xi, dtheta;
			if (stargamma > 1.6)
			{
				xi = read_vector("xsi.txt");
				dtheta = read_vector("dtheta.txt");
			}
			else
			{
				xi = read_vector("xsin3.txt");
				dtheta = read_vector("dthetan3.txt");
			}
			size_t n = xi.size();
			rhoc_ = -M*xi[n - 1] / (4 * M_PI*R*R*R*dtheta[n - 1]);
			alpha_ = R / xi[n - 1];

			size_t N = edges.size()-1;
			acc_.resize(N);
			for (size_t i = 0; i < N; ++i)
			{
				double x = 0.5*(edges[i + 1] + edges[i]);
				double xsi = x / alpha_;
				if (xsi > xi.back())
					xsi = xi.back() - 0.000001;
				double Mtemp = -rhoc_*alpha_*alpha_*alpha_ * 4 * M_PI*LinearInterpolation(xi, dtheta, xsi)*xsi*xsi;
				acc_[i] = -Mtemp / (x*x);
			}
		}

		void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double time,
			vector<Extensive> & extensives, double dt)const
		{
			size_t N = cells.size();
			f_ = GetTrueAnomaly(time, Rp_, Mbh_);
			for (size_t i = 0; i < N; ++i)
			{
				double acc = acc_[i];
				if (!selfgravity_)
					acc = 0;
				// Do Mbh
				double R = 2 * Rp_ / (1 + cos(f_));
				double x = 0.5*(edges[i + 1] + edges[i]);
				acc -= Mbh_*x*pow(sqrt(R*R + x*x), -3);
				extensives[i].momentum += extensives[i].mass*acc*dt;
				extensives[i].energy += extensives[i].mass*acc*dt*cells[i].velocity;
			}
		}
	};

	double read_number(string const& fname)
	{
		double buf = 0;
		ifstream f(fname.c_str());
		assert(f || missing_file_data(fname));
		f >> buf;
		f.close();
		return buf;
	}

  class RawInputData
  {
  public:

    const double beta;
    const double star_gamma;
    const double gas_gamma;
    const bool self_gravity;
    const string output_path;

    RawInputData
    (const double& beta_i,
     const double& star_gamma_i,
     const double& gas_gamma_i,
     const bool& self_gravity_i,
     const string& output_path_i):
      beta(beta_i),
      star_gamma(star_gamma_i),
      gas_gamma(gas_gamma_i),
      self_gravity(self_gravity_i),
      output_path(output_path_i) {}
  };

  string read_string(const string& fname)
  {
    string res = ".";
    ifstream f(fname.c_str());
    assert(f);
    assert(f.is_open());
    f >> res;
    f.close();
    return res;
  }

  RawInputData read_input(const string& input_path)
  {
    const double beta = read_number(input_path+"/beta.txt");
    const double star_gamma = 
      read_number(input_path+"/star_gamma.txt");
    const double gas_gamma =
      read_number(input_path+"/gas_gamma.txt");
    const double sg = read_number(input_path+"/selfgravity.txt");
    const string output_path = read_string(input_path+"/output_dir.txt");
    return RawInputData
      (beta,
       star_gamma,
       gas_gamma,
       sg>0.5,
       output_path);
  }

  class SimData
  {
  public:
    SimData(void):
      rid_(read_input(".")),
      cfl_(0.2),
      rs_(rid_.gas_gamma),
      eos_(rid_.gas_gamma),
      Np_(512),
      R_(1),
      M_(1),
      edges_(getedges(Np_,R_*1.01)),
      bl_(),
      br_
      (Primitive
       (1e-25, 
	1e-26, 
	0, 
	eos_.dp2s(1e-25, 1e-26))),
      boundary_(bl_, br_),
      interp_(boundary_),
      Nemden_(1./rid_.star_gamma-1.0),
      cells_
      (calc_init
       (edges_,
	M_,
	R_,
	Nemden_,
	rid_.star_gamma)),
      Mbh_(1e6),
      Rt_(R_*pow(Mbh_/M_,1.0/3.0)),
      Rp_(Rt_/rid_.beta),
      source_
      (M_,
       R_,
       Mbh_,
       Rp_,
       rid_.self_gravity,
       rid_.star_gamma,
       edges_),
      sim_
      (cfl_,
       cells_,
       edges_,
       interp_,
       eos_,
       rs_,
       source_)
    {}

    double getCFL(void) const
    {
      return cfl_;
    }

    const ExactRS& getRS(void) const
    {
      return rs_;
    }

    const IdealGas& getEOS(void) const
    {
      return eos_;
    }

    const vector<double>& getEdges(void) const
    {
      return edges_;
    }

    const MinMod& getInterp(void) const
    {
      return interp_;
    }

    const vector<Primitive>& getCells(void) const
    {
      return cells_;
    }

    const Gravity& getSourceTerm(void) const
    {
      return source_;
    }

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const RawInputData rid_;
    const double cfl_;
    const ExactRS rs_;
    const IdealGas eos_;
    const size_t Np_;
    const double R_;
    const double M_;
    const vector<double> edges_;
    const RigidWall bl_;
    const ConstantPrimitive br_;
    const SeveralBoundary boundary_;
    const MinMod interp_;
    const double Nemden_;
    const vector<Primitive> cells_;
    const double Mbh_;
    const double Rt_;
    const double Rp_;
    const Gravity source_;
    hdsim sim_;
  };
}

int main(void)
{
	// Units G=1 M=solar R=solar t=1.592657944577715e+03
	RawInputData raw_input_data = read_input(".");
	SimData sim_data;
	double R = 1;
	double M = 1;
	double Mbh = 1e6;
	double Rt = R*pow(Mbh / M, 1.0 / 3.0);
	double fstart = -acos(2 / raw_input_data.beta - 1);
	double tstart = sqrt(2 * pow(Rt / raw_input_data.beta, 3) / 
			     Mbh)*tan(fstart / 2)*
	  (3 + pow(tan(fstart / 2), 2)) / 3;
	hdsim& sim = sim_data.getSim();
	sim.SetTime(tstart);

	double dt = 0.05;
	double initd = sim_data.getCells().front().density;
	double maxd = sim_data.getCells().front().density;
	double last = sim.GetTime();
	double mind = maxd;
	int counter = 0;

	while (sim.GetCells()[0].density> 
	       max(0.25*initd,0.1*maxd) && 
	       sim.GetTime()<0.6)
	{
		if (sim.GetCycle() % 500 == 0)
			write_snapshot_to_hdf5(sim, "temp.h5");
		if (sim.GetCycle() % 100 == 0)
			cout << "Time = " << sim.GetTime() << " Cycle = " << sim.GetCycle() << endl;
		sim.TimeAdvance2();
		if (sim.GetTime() - last > dt || sim.GetCycle() == 0 || sim.GetCells()[0].density>1.02*maxd || 
			sim.GetCells()[0].density*1.02<mind)
		{
		  write_snapshot_to_hdf5
		    (sim,
		     raw_input_data.output_path+"/tide_" + 
		     int2str(counter) + ".h5");
		  last = sim.GetTime();
		  ++counter;
		  maxd = max(maxd, sim.GetCells()[0].density);
		  mind = sim.GetCells()[0].density;
		}
	}
	return 0;
}
