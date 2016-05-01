#include "hdf_util.hpp"

using namespace H5;

void write_std_vector_to_hdf5
(const CommonFG& file,
	const vector<double>& data,
	const string& caption)
{
	FloatType datatype(PredType::NATIVE_DOUBLE);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
		(file,
			data,
			caption,
			datatype);
}

void write_std_vector_to_hdf5
(const CommonFG& file,
	const vector<int>& data,
	const string& caption)
{
	IntType datatype(PredType::NATIVE_INT);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
		(file,
			data,
			caption,
			datatype);
}


Snapshot::Snapshot(void) :
	edges(),
	cells(),
	time(),
	cycle() {}

Snapshot::Snapshot(const Snapshot& source) :
	edges(source.edges),
	cells(source.cells),
	time(source.time),
	cycle(source.cycle) {}

namespace 
{
	template<class T> vector<T> read_vector_from_hdf5
		(const CommonFG& file,
			const string& caption,
			const DataType& datatype)
	{
		DataSet dataset = file.openDataSet(caption);
		DataSpace filespace = dataset.getSpace();
		hsize_t dims_out[2];
		filespace.getSimpleExtentDims(dims_out, NULL);
		const size_t NX = static_cast<size_t>(dims_out[0]);
		vector<T> result(NX);
		dataset.read(&result[0], datatype);
		return result;
	}

	vector<double> read_double_vector_from_hdf5
		(CommonFG& file, string const& caption)
	{
		return read_vector_from_hdf5<double>
			(file,
				caption,
				PredType::NATIVE_DOUBLE);
	}

	vector<int> read_int_vector_from_hdf5
		(const CommonFG& file,
			const string& caption)
	{
		return read_vector_from_hdf5<int>
			(file,
				caption,
				PredType::NATIVE_INT);
	}
}


void write_snapshot_to_hdf5(hdsim const& sim, string const& fname)
{
	H5File file(H5std_string(fname), H5F_ACC_TRUNC);
	Group geometry = file.createGroup("/geometry");
	Group hydrodynamic = file.createGroup("/hydrodynamic");

	// General
	write_std_vector_to_hdf5
		(file,
			vector<double>(1, sim.GetTime()),
			"time");
	write_std_vector_to_hdf5
		(file,
		 vector<int>(1, static_cast<int>(sim.GetCycle())),
			"cycle");

	// Geometry  
	write_std_vector_to_hdf5
		(geometry,sim.GetEdges(),"edges");
	
	// Hydrodynamic
	vector<Primitive> const& cells = sim.GetCells();
	size_t N = cells.size();
	vector<double> density(N), pressure(N), velocity(N);
	for (size_t i = 0; i < N; ++i)
	{
		density[i] = cells[i].density;
		pressure[i] = cells[i].pressure;
		velocity[i] = cells[i].velocity;
	}

	write_std_vector_to_hdf5
		(hydrodynamic,density,
			"density");
	write_std_vector_to_hdf5
		(hydrodynamic,
			pressure,
			"pressure");
	write_std_vector_to_hdf5
		(hydrodynamic,
			velocity,
			"velocity");
}

Snapshot read_hdf5_snapshot
(const string &fname)
{
	Snapshot res;
	H5File file(fname, H5F_ACC_RDONLY);
	Group g_geometry = file.openGroup("geometry");
	Group g_hydrodynamic = file.openGroup("hydrodynamic");

	// Mesh points
	res.edges = read_double_vector_from_hdf5(g_geometry, "edges");

	// Hydrodynamic
	{
		const vector<double> density =
			read_double_vector_from_hdf5(g_hydrodynamic, "density");
		const vector<double> pressure =
			read_double_vector_from_hdf5(g_hydrodynamic, "pressure");
		const vector<double> velocity =
			read_double_vector_from_hdf5(g_hydrodynamic, "velocity");
		res.cells.resize(density.size());
		for (size_t i = 0; i<res.cells.size(); ++i) 
		{
			res.cells.at(i).density = density.at(i);
			res.cells.at(i).pressure = pressure.at(i);
			res.cells.at(i).velocity = velocity.at(i);
		}
	}

	// Misc
	{
		const vector<double> time =
			read_double_vector_from_hdf5(file, "time");
		res.time = time.at(0);
		const vector<int> cycle =
			read_int_vector_from_hdf5(file, "cycle");
		res.cycle = cycle.at(0);
	}
	return res;
}

