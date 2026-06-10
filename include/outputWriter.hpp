#pragma once

#include "ErrorHandler.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cstring>    // memcpy
#include <cstdint>

// ============================================================
//  VTUWriter — XML VTK UnstructuredGrid output (.vtu)
//              for ParaView (modern format, replaces .vtk)
//
//  Structured 2D/3D grid written as an unstructured grid:
//    2D → VTK_QUAD cells        (type=9)
//    3D → VTK_HEXAHEDRON cells  (type=12)
//
//  Supports two encodings:
//    ASCII  → human-readable, convenient for debugging
//    Base64 → inline binary, files ~3x smaller and faster to read
//
//  Minimal usage (2D, ASCII):
//    VTUWriter w("out_0000.vtu", nx, ny, 1, dx, dy, 1.0);
//    w.addScalar("phi", phi_vec);
//    w.write();
//
//  With Base64:
//    VTUWriter w("out.vtu", nx, ny, 1, dx, dy, 1.0,
//                0,0,0, true, VTUWriter::BASE64);
//    w.addScalar("phi", phi);
//    w.write();
// ============================================================

// ---- filename helper ----
inline std::string makeVTUFilename(const std::string& prefix, int step, int width = 6) {
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(width) << std::setfill('0') << step << ".vtu";
    return oss.str();
}

class VTUWriter {
public:
    enum Encoding { ASCII, BASE64 };

    // nx, ny, nz : number of CELLS per dimension (nz=1 for 2D)
    // cellCentered: true  → CELL_DATA
    //               false → POINT_DATA
    VTUWriter(const std::string& filename,
              int nx, int ny, int nz,
              double dx, double dy, double dz,
              double x0 = 0.0, double y0 = 0.0, double z0 = 0.0,
              bool cellCentered = false,
              Encoding enc = ASCII)
        : filename_(filename),
          nx_(nx), ny_(ny), nz_(nz),
          dx_(dx), dy_(dy), dz_(dz),
          x0_(x0), y0_(y0), z0_(z0),
          cellCentered_(cellCentered),
          enc_(enc)
    {
        is2D_   = (nz_ == 1);
        nCells_ = nx_ * ny_ * nz_;
        // nodes: (nx+1)*(ny+1)*(nz+1)
        nNodes_ = (nx_+1) * (ny_+1) * (is2D_ ? 1 : (nz_+1));
    }

    void addScalar(const std::string& name, const std::vector<double>& data) {
        int expected = cellCentered_ ? nCells_ : nNodes_;
        if ((int)data.size() != expected)
            throw std::runtime_error("VTUWriter::addScalar — wrong size for '" + name + "'");
        scalars_.push_back({name, data});
    }

    void addVector(const std::string& name,
                   const std::vector<double>& u,
                   const std::vector<double>& v,
                   const std::vector<double>& w = {})
    {
        int expected = cellCentered_ ? nCells_ : nNodes_;
        if ((int)u.size() != expected || (int)v.size() != expected)
            throw std::runtime_error("VTUWriter::addVector — wrong size for '" + name + "'");
        std::vector<double> ww = w.empty() ? std::vector<double>(expected, 0.0) : w;
        vectors_.push_back({name, u, v, ww});
    }

    void write() const {
        std::ofstream f(filename_);
        if (!f.is_open())
        	throw CannotOpenFile(filename_,"output");

        f << std::scientific << std::setprecision(8);

        // XML header
        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "  <UnstructuredGrid>\n";
        f << "    <Piece NumberOfPoints=\"" << nNodes_
          << "\" NumberOfCells=\"" << nCells_ << "\">\n";

        writePoints(f);
        writeCells(f);
        writeData(f);

        f << "    </Piece>\n";
        f << "  </UnstructuredGrid>\n";
        f << "</VTKFile>\n";
        f.close();
    }

private:
    struct ScalarField { std::string name; std::vector<double> data; };
    struct VectorField { std::string name; std::vector<double> u, v, w; };

    std::string filename_;
    int nx_, ny_, nz_, nCells_, nNodes_;
    double dx_, dy_, dz_, x0_, y0_, z0_;
    bool cellCentered_, is2D_;
    Encoding enc_;
    std::vector<ScalarField> scalars_;
    std::vector<VectorField> vectors_;

    // ---- node (i,j,k) → flat index ----
    //  k varies last (z), i first (x) → row-major x-fast
    int nodeIndex(int i, int j, int k) const {
        return k*(ny_+1)*(nx_+1) + j*(nx_+1) + i;
    }

    // ---- cells: same ordering as the phi array ----
    int cellIndex(int i, int j, int k) const {
        return k*ny_*nx_ + j*nx_ + i;
    }

    // ---- POINTS ----
    void writePoints(std::ofstream& f) const {
        f << "      <Points>\n";
        f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\""
          << fmtStr() << "\">\n";

        std::vector<double> coords;
        coords.reserve(nNodes_ * 3);
        int kmax = is2D_ ? 0 : nz_;
		for (int k = 0; k <= kmax; ++k)
        for (int i = 0; i <= nx_; ++i)
        for (int j = 0; j <= ny_; ++j) {
            coords.push_back(x0_ + i * dx_);
            coords.push_back(y0_ + j * dy_);
            coords.push_back(z0_ + k * dz_);
        }

        writeArray(f, coords);
        f << "        </DataArray>\n";
        f << "      </Points>\n";
    }

    // ---- CELLS (connectivity, offsets, types) ----
    void writeCells(std::ofstream& f) const {
        f << "      <Cells>\n";

        // --- connectivity ---
        f << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\""
          << fmtStr() << "\">\n";
        {
            std::vector<double> conn;
            conn.reserve(nCells_ * (is2D_ ? 4 : 8));
            for (int k = 0; k < nz_; ++k)
            for (int i = 0; i < nx_; ++i) 
            for (int j = 0; j < ny_; ++j) {
                if (is2D_) {
                    // VTK_QUAD, nodes in counter-clockwise order
                    conn.push_back(nodeIndex(i,   j,   0));
                    conn.push_back(nodeIndex(i+1, j,   0));
                    conn.push_back(nodeIndex(i+1, j+1, 0));
                    conn.push_back(nodeIndex(i,   j+1, 0));
                } else {
                    // VTK_HEXAHEDRON: bottom face first, then top (same order)
                    conn.push_back(nodeIndex(i,   j,   k));
                    conn.push_back(nodeIndex(i+1, j,   k));
                    conn.push_back(nodeIndex(i+1, j+1, k));
                    conn.push_back(nodeIndex(i,   j+1, k));
                    conn.push_back(nodeIndex(i,   j,   k+1));
                    conn.push_back(nodeIndex(i+1, j,   k+1));
                    conn.push_back(nodeIndex(i+1, j+1, k+1));
                    conn.push_back(nodeIndex(i,   j+1, k+1));
                }
            }
            writeArray(f, conn, true);
        }
        f << "        </DataArray>\n";

        // --- offsets ---
        f << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\""
          << fmtStr() << "\">\n";
        {
            int nodesPerCell = is2D_ ? 4 : 8;
            std::vector<double> off(nCells_);
            for (int c = 0; c < nCells_; ++c)
                off[c] = (c+1) * nodesPerCell;
            writeArray(f, off, true);
        }
        f << "        </DataArray>\n";

        // --- types ---
        f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\""
          << fmtStr() << "\">\n";
        {
            int cellType = is2D_ ? 9 : 12;   // QUAD=9, HEXAHEDRON=12
            std::vector<double> types(nCells_, cellType);
            writeArray(f, types, true);
        }
        f << "        </DataArray>\n";

        f << "      </Cells>\n";
    }

    // ---- DATA (scalars + vectors) ----
    void writeData(std::ofstream& f) const {
        if (scalars_.empty() && vectors_.empty()) return;

        // First scalar/vector name used as default colormap in ParaView
        std::string defaultScalar = scalars_.empty() ? "" : scalars_[0].name;
        std::string defaultVector = vectors_.empty() ? "" : vectors_[0].name;

        std::string tag = cellCentered_ ? "CellData" : "PointData";

        f << "      <" << tag;
        if (!defaultScalar.empty()) f << " Scalars=\"" << defaultScalar << "\"";
        if (!defaultVector.empty()) f << " Vectors=\"" << defaultVector << "\"";
        f << ">\n";

        for (const auto& s : scalars_) {
            f << "        <DataArray type=\"Float64\" Name=\"" << s.name
              << "\" NumberOfComponents=\"1\" format=\"" << fmtStr() << "\">\n";
            writeArray(f, s.data);
            f << "        </DataArray>\n";
        }

        for (const auto& vec : vectors_) {
            f << "        <DataArray type=\"Float64\" Name=\"" << vec.name
              << "\" NumberOfComponents=\"3\" format=\"" << fmtStr() << "\">\n";
            int n = (int)vec.u.size();
            std::vector<double> interleaved(n * 3);
            for (int i = 0; i < n; ++i) {
                interleaved[3*i+0] = vec.u[i];
                interleaved[3*i+1] = vec.v[i];
                interleaved[3*i+2] = vec.w[i];
            }
            writeArray(f, interleaved);
            f << "        </DataArray>\n";
        }

        f << "      </" << tag << ">\n";
    }

    // ---- write values: ASCII or Base64 ----
    void writeArray(std::ofstream& f, const std::vector<double>& v, bool asInt = false) const {
        if (enc_ == ASCII) {
            f << "          ";
            if (asInt){
	            for (size_t i = 0; i < v.size(); ++i) {
	                f << std::setw(8) << static_cast<int>(v[i]);
	                f << ((i % 8 == 7) ? "\n          " : " ");
	            }	
            }
            else{
	            for (size_t i = 0; i < v.size(); ++i) {
	                f << v[i];
	                f << ((i % 6 == 5) ? "\n          " : " ");
	            }
            }
            f << "\n";
        } else {
            // base64: prepend byte length (uint32) as per VTK spec
            uint32_t byteLen = (uint32_t)(v.size() * sizeof(double));
            std::vector<unsigned char> raw(4 + byteLen);
            std::memcpy(raw.data(), &byteLen, 4);
            std::memcpy(raw.data() + 4, v.data(), byteLen);
            f << "          " << base64Encode(raw) << "\n";
        }
    }

    const char* fmtStr() const { return enc_ == ASCII ? "ascii" : "binary"; }

    // ---- minimal Base64 encoder ----
    static std::string base64Encode(const std::vector<unsigned char>& in) {
        static const char* tab =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::string out;
        out.reserve(((in.size() + 2) / 3) * 4);
        for (size_t i = 0; i < in.size(); i += 3) {
            uint32_t b = (uint32_t)in[i] << 16;
            if (i+1 < in.size()) b |= (uint32_t)in[i+1] << 8;
            if (i+2 < in.size()) b |= (uint32_t)in[i+2];
            out += tab[(b >> 18) & 63];
            out += tab[(b >> 12) & 63];
            out += (i+1 < in.size()) ? tab[(b >> 6) & 63] : '=';
            out += (i+2 < in.size()) ? tab[b & 63]        : '=';
        }
        return out;
    }
};

// ============================================================
//  PVDWriter — .pvd file for time series in ParaView
//  Collects .vtu files and enables animation via the time slider.
//
//  Usage:
//    PVDWriter pvd("simulation.pvd");
//    for (int step = 0; ...) {
//        std::string vtu = makeVTUFilename("phi", step);
//        // write the vtu...
//        pvd.addStep(step * DT, vtu);
//    }
//    pvd.write();
// ============================================================
class PVDWriter {
public:
    PVDWriter(const std::string& filename) : filename_(filename) {}

    void addStep(double time, const std::string& vtuFile) {
        steps_.push_back({time, vtuFile});
    }

    void write() const {
        std::ofstream f(filename_);
        if (!f.is_open())
            throw std::runtime_error("PVDWriter: cannot open '" + filename_ + "'");
        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
        f << "  <Collection>\n";
        for (const auto& s : steps_)
            f << "    <DataSet timestep=\"" << s.time
              << "\" file=\"" << s.file << "\"/>\n";
        f << "  </Collection>\n";
        f << "</VTKFile>\n";
    }

private:
    struct Step { double time; std::string file; };
    std::string filename_;
    std::vector<Step> steps_;
};


class OutputWriter1D {
public:
    OutputWriter1D(const ConfigParser& cfg)
    {
    	filename_ = cfg.getString("OUTPUT_DIR") + cfg.getString("OUTPUT_FILE") + ".cfg";
    	nx = cfg.getInt("GRID_POINTS");
	    std::ofstream file(filename_);
	    if (!file){
	        throw CannotOpenFile(filename_,"output");
	    }
	    file << "t";
	    for (int i = 0; i < nx; i++){
	        file << ",u" << i;
	    }
	    file.close();
    }
    void saveCurrentTimeStep(int timestep, const std::vector<double> u) {
	    std::ofstream file(filename_, std::ios::app);
	    if (!file){
	        std::cerr << "Error: impossible to open file " << filename_ << std::endl;
	        return;
	    }
	    file << '\n';

	    file << timestep << ",";
	    file << u[0];
	    for (int i = 1; i < nx; i++){
	        file << "," << u[i];
	    }
	    file.close();
	    std::cout << "Step: " << timestep << " saved into: " << filename_ << std::endl;
	}
private:
	std::string filename_;

	int nx;
};
