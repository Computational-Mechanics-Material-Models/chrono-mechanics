// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// MPI multidomain explicit LDPM displacement-loading demo.
//
// This demo is intended for worker-domain scaling studies using the Chrono
// multidomain module and LDPM mesh/materials.
// =============================================================================

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChLink.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkMotorLinearPosition.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/functions/ChFunctionPoly.h"
#include "chrono/functions/ChFunctionRamp.h"
#include "chrono/functions/ChFunctionSequence.h"
#include "chrono/functions/ChFunctionPositionXYZFunctions.h"

#include "chrono_ldpm/ChBuilderLDPM.h"
#include "chrono_ldpm/ChElementLDPM.h"
#include "chrono_ldpm/ChMaterialVECT.h"
#include "chrono_ldpm/ChSectionLDPM.h"

#include "chrono_multidomain/ChDomainBuilder.h"
#include "chrono_multidomain/ChDomainManagerMPI.h"
#include "chrono_multidomain/ChSolverLumpedMultidomain.h"
#include "chrono_multidomain/ChSystemDescriptorMultidomain.h"

#include "demo_MDOM_TagUtils.h"

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::ldpm;
using namespace chrono::multidomain;

namespace {
enum class BoundaryMode { FIX, FREE };

BoundaryMode ParseBoundaryMode(std::string mode) {
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (mode == "free")
        return BoundaryMode::FREE;
    return BoundaryMode::FIX;
}

const char* BoundaryModeName(BoundaryMode mode) {
    return (mode == BoundaryMode::FREE) ? "free" : "fix";
}

struct RunSettings {
    double dt = 5.0e-7;
    double t_end = 0.4;
    double disp_rate = -5.0;
    double poly_a = -1250.0;
    double switch_time = 2.0e-3;
    double load_scale = 1.0;
    double target_top_displacement = -0.02;
    double stop_disp_frac = -1.0;
    double ref_x = 53.0;
    double ref_y = 50.0;
    double z_tol_ref = 0.2;
    double z_tol_surface = 1.0e-3;
    double constraint_penalty = 2.0e6;
    double constraint_damping = 0.0;
    int hist_stride = 1;
    int vtk_stride = 800;
    bool use_diagonal_lumping = false;
    bool do_partition_update_once = false;
    bool do_partition_update_each_step = false;
    bool disable_master_after_init = true;
    bool use_zslice_surface_owner = true;
    bool use_tag_reaction_selector = true;
    bool fix_top_rotations = false;
    BoundaryMode boundary_mode = BoundaryMode::FIX;
    bool verbose_partition = false;
    bool verbose_serialization = false;
    bool verbose_variable_updates = false;
    bool clean_output_dirs = true;
    DomainSerializerFormat serializer_format = DomainSerializerFormat::BINARY;
    ChAxis slice_axis = ChAxis::Z;
    bool use_grid_decomposition = false;
    bool use_complete_worker_graph = true;
    int grid_x = 0;
    int grid_y = 0;
    int grid_z = 0;
};

bool StartsWith(const std::string& text, const std::string& prefix) {
    return text.rfind(prefix, 0) == 0;
}

bool EndsWith(const std::string& text, const std::string& suffix) {
    if (suffix.size() > text.size())
        return false;
    return text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string TrimCopy(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

constexpr int kLDPMFacetTagBase = 1000000;
constexpr int kTagConstrBot = 90000001;
constexpr int kTagMotor = 90000002;
constexpr int kTagTopNodeMin = 91000000;
constexpr int kTagTopNodeMax = 91099999;
constexpr int kTagBottomNodeMin = 91100000;
constexpr int kTagBottomNodeMax = 91199999;
std::vector<std::array<ChVector3d, 3>> g_ldpm_facet_vertices;

bool LoadFacetVertexLookup(const std::filesystem::path& data_dir) {
    if (!g_ldpm_facet_vertices.empty())
        return true;

    const std::filesystem::path facet_file = data_dir / "LDPMgeo000-data-facets.dat";
    const std::filesystem::path facet_vertices_file = data_dir / "LDPMgeo000-data-facetsVertices.dat";

    std::ifstream fv_in(facet_vertices_file);
    if (!fv_in.is_open())
        return false;

    std::vector<ChVector3d> all_vertices;
    std::string line;
    while (std::getline(fv_in, line)) {
        const std::string t = TrimCopy(line);
        if (t.empty() || t.rfind("//", 0) == 0)
            continue;
        std::istringstream iss(t);
        double x, y, z;
        if (!(iss >> x >> y >> z))
            continue;
        all_vertices.emplace_back(x, y, z);
    }
    fv_in.close();

    if (all_vertices.empty())
        return false;

    std::ifstream f_in(facet_file);
    if (!f_in.is_open())
        return false;

    while (std::getline(f_in, line)) {
        const std::string t = TrimCopy(line);
        if (t.empty() || t.rfind("//", 0) == 0)
            continue;

        std::istringstream iss(t);
        int tet, id0, id1, id2;
        if (!(iss >> tet >> id0 >> id1 >> id2))
            continue;

        if (id0 >= 0 && id1 >= 0 && id2 >= 0 && static_cast<size_t>(id0) < all_vertices.size() &&
            static_cast<size_t>(id1) < all_vertices.size() && static_cast<size_t>(id2) < all_vertices.size()) {
            g_ldpm_facet_vertices.push_back({all_vertices[static_cast<size_t>(id0)], all_vertices[static_cast<size_t>(id1)],
                                             all_vertices[static_cast<size_t>(id2)]});
        } else {
            g_ldpm_facet_vertices.push_back({ChVector3d(0, 0, 0), ChVector3d(0, 0, 0), ChVector3d(0, 0, 0)});
        }
    }
    f_in.close();

    return !g_ldpm_facet_vertices.empty();
}

std::array<int, 3> ComputeGridDivisions(int domains) {
    int best_x = 1;
    int best_y = 1;
    int best_z = domains;
    int best_spread = std::numeric_limits<int>::max();
    for (int nx = 1; nx <= domains; ++nx) {
        if (domains % nx != 0)
            continue;
        const int rem = domains / nx;
        for (int ny = 1; ny <= rem; ++ny) {
            if (rem % ny != 0)
                continue;
            const int nz = rem / ny;
            const int spread = std::max({nx, ny, nz}) - std::min({nx, ny, nz});
            if (spread < best_spread ||
                (spread == best_spread && std::tie(nx, ny, nz) > std::tie(best_x, best_y, best_z))) {
                best_spread = spread;
                best_x = nx;
                best_y = ny;
                best_z = nz;
            }
        }
    }
    return {best_x, best_y, best_z};
}

RunSettings ParseArgs(int argc, char* argv[]) {
    RunSettings cfg;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (StartsWith(arg, "--dt=")) {
            cfg.dt = std::stod(arg.substr(std::string("--dt=").size()));
            continue;
        }
        if (StartsWith(arg, "--t_end=")) {
            cfg.t_end = std::stod(arg.substr(std::string("--t_end=").size()));
            continue;
        }
        if (arg == "--pressure") {
            continue;
        }
        if (StartsWith(arg, "--pressure=")) {
            continue;
        }
        if (StartsWith(arg, "--target_pressure=")) {
            continue;
        }
        if (StartsWith(arg, "--pressure_ramp=")) {
            continue;
        }
        if (StartsWith(arg, "--target_disp=")) {
            cfg.target_top_displacement = std::stod(arg.substr(std::string("--target_disp=").size()));
            continue;
        }
        if (StartsWith(arg, "--disp_rate=")) {
            cfg.disp_rate = std::stod(arg.substr(std::string("--disp_rate=").size()));
            continue;
        }
        if (StartsWith(arg, "--poly_a=")) {
            cfg.poly_a = std::stod(arg.substr(std::string("--poly_a=").size()));
            continue;
        }
        if (StartsWith(arg, "--switch_time=")) {
            cfg.switch_time = std::stod(arg.substr(std::string("--switch_time=").size()));
            continue;
        }
        if (StartsWith(arg, "--load_scale=")) {
            cfg.load_scale = std::stod(arg.substr(std::string("--load_scale=").size()));
            continue;
        }
        if (StartsWith(arg, "--disp=")) {
            cfg.target_top_displacement = std::stod(arg.substr(std::string("--disp=").size()));
            continue;
        }
        if (StartsWith(arg, "--stop_disp_frac=")) {
            cfg.stop_disp_frac = std::stod(arg.substr(std::string("--stop_disp_frac=").size()));
            continue;
        }
        if (StartsWith(arg, "--hist_stride=")) {
            cfg.hist_stride = std::max(0, std::stoi(arg.substr(std::string("--hist_stride=").size())));
            continue;
        }
        if (StartsWith(arg, "--vtk_stride=")) {
            cfg.vtk_stride = std::max(0, std::stoi(arg.substr(std::string("--vtk_stride=").size())));
            continue;
        }
        if (StartsWith(arg, "--constraint_penalty=")) {
            cfg.constraint_penalty = std::stod(arg.substr(std::string("--constraint_penalty=").size()));
            continue;
        }
        if (StartsWith(arg, "--constraint_damping=")) {
            cfg.constraint_damping = std::stod(arg.substr(std::string("--constraint_damping=").size()));
            continue;
        }
        if (StartsWith(arg, "--serializer=")) {
            const std::string mode = arg.substr(std::string("--serializer=").size());
            if (mode == "json") {
                cfg.serializer_format = DomainSerializerFormat::JSON;
            } else if (mode == "xml") {
                cfg.serializer_format = DomainSerializerFormat::XML;
            } else {
                cfg.serializer_format = DomainSerializerFormat::BINARY;
            }
            continue;
        }
        if (StartsWith(arg, "--slice_axis=")) {
            const std::string axis = arg.substr(std::string("--slice_axis=").size());
            if (axis == "x" || axis == "X") {
                cfg.slice_axis = ChAxis::X;
            } else if (axis == "y" || axis == "Y") {
                cfg.slice_axis = ChAxis::Y;
            } else {
                cfg.slice_axis = ChAxis::Z;
            }
            continue;
        }
        if (StartsWith(arg, "--grid=")) {
            std::string dimensions = arg.substr(std::string("--grid=").size());
            std::replace(dimensions.begin(), dimensions.end(), 'x', ' ');
            std::replace(dimensions.begin(), dimensions.end(), 'X', ' ');
            std::replace(dimensions.begin(), dimensions.end(), ',', ' ');
            std::istringstream stream(dimensions);
            if (!(stream >> cfg.grid_x >> cfg.grid_y >> cfg.grid_z) || cfg.grid_x < 1 || cfg.grid_y < 1 ||
                cfg.grid_z < 1) {
                throw std::invalid_argument("--grid requires three positive dimensions, for example --grid=2x2x2");
            }
            cfg.use_grid_decomposition = true;
            continue;
        }
        if (arg == "--neighbor_worker_graph") {
            cfg.use_complete_worker_graph = false;
            continue;
        }
        if (arg == "--no_diagonal_lumping") {
            cfg.use_diagonal_lumping = false;
            continue;
        }
        if (arg == "--diagonal_lumping") {
            cfg.use_diagonal_lumping = true;
            continue;
        }

        if (arg == "--verbose_partition") {
            cfg.verbose_partition = true;
            continue;
        }
        if (arg == "--verbose_serialization") {
            cfg.verbose_serialization = true;
            continue;
        }
        if (arg == "--verbose_variables") {
            cfg.verbose_variable_updates = true;
            continue;
        }
        if (arg == "--no_clean_output") {
            cfg.clean_output_dirs = false;
            continue;
        }
        if (arg == "--clean_output") {
            cfg.clean_output_dirs = true;
            continue;
        }
        if (arg == "--no_partition_update") {
            cfg.do_partition_update_once = false;
            continue;
        }
        if (arg == "--partition_update") {
            cfg.do_partition_update_once = true;
            continue;
        }
        if (arg == "--partition_update_each_step") {
            cfg.do_partition_update_each_step = true;
            continue;
        }
        if (arg == "--keep_master_enabled") {
            cfg.disable_master_after_init = false;
            continue;
        }
        if (arg == "--no_zslice_surface_owner") {
            cfg.use_zslice_surface_owner = false;
            continue;
        }
        if (arg == "--zslice_surface_owner") {
            cfg.use_zslice_surface_owner = true;
            continue;
        }
        if (arg == "--no_tag_reaction_selector") {
            cfg.use_tag_reaction_selector = false;
            continue;
        }
        if (arg == "--tag_reaction_selector") {
            cfg.use_tag_reaction_selector = true;
            continue;
        }
        if (arg == "--fix_top_rotations") {
            cfg.fix_top_rotations = true;
            continue;
        }
        if (StartsWith(arg, "--boundary_mode=")) {
            cfg.boundary_mode = ParseBoundaryMode(arg.substr(std::string("--boundary_mode=").size()));
            continue;
        }
        if (arg == "--free_boundary") {
            cfg.boundary_mode = BoundaryMode::FREE;
            continue;
        }
        if (arg == "--fix_boundary") {
            cfg.boundary_mode = BoundaryMode::FIX;
            continue;
        }
    }

    return cfg;
}

std::shared_ptr<ChNodeFEAxyzrot> FindClosestNodeAtPlane(std::shared_ptr<ChMesh> mesh,
                                                        double z_target,
                                                        double z_tol,
                                                        double x_ref,
                                                        double y_ref) {
    std::shared_ptr<ChNodeFEAxyzrot> best;
    double best_dist2 = std::numeric_limits<double>::max();

    for (auto& node_base : mesh->GetNodes()) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
        if (!node)
            continue;

        const ChVector3d p = node->GetPos();
        if (std::abs(p.z() - z_target) > z_tol)
            continue;

        const double dx = p.x() - x_ref;
        const double dy = p.y() - y_ref;
        const double dist2 = dx * dx + dy * dy;

        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best = node;
        }
    }

    return best;
}

std::vector<std::shared_ptr<ChNodeFEAxyzrot>> CollectSurfaceNodes(std::shared_ptr<ChMesh> mesh,
                                                                  double z_target,
                                                                  double z_tol,
                                                                  std::shared_ptr<ChNodeFEAxyzrot> exclude) {
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> out;

    for (auto& node_base : mesh->GetNodes()) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
        if (!node)
            continue;

        if (node == exclude)
            continue;

        const ChVector3d p = node->GetPos();
        if (std::abs(p.z() - z_target) <= z_tol)
            out.push_back(node);
    }

    return out;
}

std::shared_ptr<ChNodeFEAxyzrot> FindNodeByTag(ChSystem& sys, int tag) {
    for (auto& mesh : sys.GetMeshes()) {
        for (auto& node_base : mesh->GetNodes()) {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
            if (node && node->GetTag() == tag)
                return node;
        }
    }

    return nullptr;
}

std::shared_ptr<ChLinkMotorLinearPosition> FindMotorByTag(ChSystem& sys, int tag) {
    for (auto& link_item : sys.GetLinks()) {
        auto motor = std::dynamic_pointer_cast<ChLinkMotorLinearPosition>(link_item);
        if (motor && motor->GetTag() == tag)
            return motor;
    }
    return nullptr;
}

std::shared_ptr<ChLinkMateGeneric> FindMateByTag(ChSystem& sys, int tag) {
    for (auto& link_item : sys.GetLinks()) {
        auto mate = std::dynamic_pointer_cast<ChLinkMateGeneric>(link_item);
        if (mate && mate->GetTag() == tag)
            return mate;
    }
    return nullptr;
}

std::shared_ptr<ChNodeFEAxyzrot> FindClosestNodeAtPlaneInSystem(ChSystem& sys,
                                                                double z_target,
                                                                double z_tol,
                                                                double x_ref,
                                                                double y_ref) {
    std::shared_ptr<ChNodeFEAxyzrot> best;
    double best_dist2 = std::numeric_limits<double>::max();

    for (auto& mesh : sys.GetMeshes()) {
        auto candidate = FindClosestNodeAtPlane(mesh, z_target, z_tol, x_ref, y_ref);
        if (!candidate)
            continue;

        const ChVector3d p = candidate->GetPos();
        const double dx = p.x() - x_ref;
        const double dy = p.y() - y_ref;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best = candidate;
        }
    }

    return best;
}

ChNodeFEAxyzrot* FindClosestNodeAtPlaneInLinks(ChSystem& sys,
                                               double z_target,
                                               double z_tol,
                                               double x_ref,
                                               double y_ref) {
    ChNodeFEAxyzrot* best = nullptr;
    double best_dist2 = std::numeric_limits<double>::max();

    auto try_candidate = [&](ChNodeFEAxyzrot* node) {
        if (!node)
            return;
        const ChVector3d p = node->GetPos();
        if (std::abs(p.z() - z_target) > z_tol)
            return;
        const double dx = p.x() - x_ref;
        const double dy = p.y() - y_ref;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best = node;
        }
    };

    for (auto& link_item : sys.GetLinks()) {
        auto link = std::dynamic_pointer_cast<ChLink>(link_item);
        if (!link)
            continue;

        try_candidate(dynamic_cast<ChNodeFEAxyzrot*>(link->GetBody1()));
        try_candidate(dynamic_cast<ChNodeFEAxyzrot*>(link->GetBody2()));
    }

    return best;
}

void ApplyLDPMMaterialParameters(const std::shared_ptr<ChMaterialVECT>& mat) {
    mat->Set_density(2.338E-9);
    mat->Set_E0(60273);
    mat->Set_alpha(0.25);
    mat->Set_sigmat(3.44);
    mat->Set_sigmas(8.944);
    mat->Set_nt(0.4);
    mat->Set_lt(500);
    mat->Set_Ed(60273);
    mat->Set_sigmac0(150);
    mat->Set_Hc0(24109);
    mat->Set_Hc1(6027.3);
    mat->Set_beta(0);
    mat->Set_kc0(4);
    mat->Set_kc1(1);
    mat->Set_kc2(5);
    mat->Set_kc3(0.1);
    mat->Set_mu0(0.4);
    mat->Set_muinf(0);
    mat->Set_sigmaN0(600);
}

void TagLDPMSharedObjects(const std::shared_ptr<ChMesh>& mesh, const std::shared_ptr<ChMaterialVECT>& mat, int start_tag) {
    int next_tag = start_tag;
    mat->SetTag(next_tag++);

    for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie));
        if (!elem)
            continue;

        if (!elem->GetMaterial())
            elem->SetMaterial(mat);

        for (auto& facet : elem->GetSection()) {
            facet->SetTag(next_tag++);
            facet->Set_material(mat);
        }
    }
}

bool MeshHasLDPMElements(const std::shared_ptr<ChMesh>& mesh) {
    for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
        if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie)))
            return true;
    }
    return false;
}

bool WriteLDPMMeshVTK(const std::shared_ptr<ChMesh>& mesh, const std::filesystem::path& vtk_path) {
    std::unordered_map<ChNodeFEAxyzrot*, int> node_lut;
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodes;
    std::vector<std::array<int, 4>> tets;

    for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie));
        if (!elem)
            continue;

        std::array<int, 4> tet = {-1, -1, -1, -1};
        bool valid = true;
        for (int j = 0; j < 4; ++j) {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(elem->GetNode(j));
            if (!node) {
                valid = false;
                break;
            }

            auto it = node_lut.find(node.get());
            if (it == node_lut.end()) {
                const int new_id = static_cast<int>(nodes.size());
                node_lut.emplace(node.get(), new_id);
                nodes.push_back(node);
                tet[j] = new_id;
            } else {
                tet[j] = it->second;
            }
        }

        if (valid)
            tets.push_back(tet);
    }

    if (tets.empty())
        return false;

    std::ofstream out(vtk_path, std::ios::out);
    if (!out.is_open())
        return false;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << tets.size() << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        const ChVector3d p = node->GetPos();
        out << "          " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& tet : tets) {
        out << "          " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (size_t i = 0; i < tets.size(); ++i) {
        offset += 4;
        out << "          " << offset << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < tets.size(); ++i) {
        out << "          10\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
    out << "      <PointData Vectors=\"displacement\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        const ChVector3d d = node->GetPos() - node->GetX0().GetPos();
        out << "          " << d.x() << " " << d.y() << " " << d.z() << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        const ChVector3d v = node->GetPosDt();
        out << "          " << v.x() << " " << v.y() << " " << v.z() << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </PointData>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    return true;
}

bool WriteLDPMCrackVTK(const std::shared_ptr<ChMesh>& mesh, const std::filesystem::path& vtk_path) {
    auto map_ref_to_deformed = [](const std::array<ChVector3d, 4>& x_ref, const std::array<ChVector3d, 4>& x_cur,
                                  const ChVector3d& p_ref, ChVector3d& p_cur) {
        const auto det3 = [](const ChVector3d& a, const ChVector3d& b, const ChVector3d& c) {
            return a.x() * (b.y() * c.z() - b.z() * c.y()) - a.y() * (b.x() * c.z() - b.z() * c.x()) +
                   a.z() * (b.x() * c.y() - b.y() * c.x());
        };

        const ChVector3d v0 = x_ref[1] - x_ref[0];
        const ChVector3d v1 = x_ref[2] - x_ref[0];
        const ChVector3d v2 = x_ref[3] - x_ref[0];
        const ChVector3d rhs = p_ref - x_ref[0];

        const double det = det3(v0, v1, v2);
        if (std::abs(det) < 1e-18) {
            p_cur = p_ref;
            return false;
        }

        const double w1 = det3(rhs, v1, v2) / det;
        const double w2 = det3(v0, rhs, v2) / det;
        const double w3 = det3(v0, v1, rhs) / det;
        const double w0 = 1.0 - w1 - w2 - w3;

        // Validate that the point is inside (or very near) the reference tetrahedron.
        // If not, this facet->tet association is likely mismatched and we should fall back.
        const double eps = 1e-6;
        if (w0 < -eps || w1 < -eps || w2 < -eps || w3 < -eps ||
            w0 > 1.0 + eps || w1 > 1.0 + eps || w2 > 1.0 + eps || w3 > 1.0 + eps) {
            p_cur = p_ref;
            return false;
        }

        p_cur = w0 * x_cur[0] + w1 * x_cur[1] + w2 * x_cur[2] + w3 * x_cur[3];
        return true;
    };

    struct FacetTri {
        ChVector3d p0;
        ChVector3d p1;
        ChVector3d p2;
        double crack;
    };
    std::vector<FacetTri> facet_tris;

    for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie));
        if (!elem)
            continue;

        std::array<ChVector3d, 4> tet_x0;
        std::array<ChVector3d, 4> tet_x;
        bool have_tet_nodes = true;
        for (int k = 0; k < 4; ++k) {
            auto node = elem->GetTetrahedronNode(k);
            if (!node) {
                have_tet_nodes = false;
                break;
            }
            tet_x0[k] = node->GetX0().GetPos();
            tet_x[k] = node->GetPos();
        }

        const auto& sections = elem->GetSection();
        const int nfacets = std::min<int>(12, static_cast<int>(sections.size()));
        for (int j = 0; j < nfacets; ++j) {
            auto facet = sections[j];
            if (!facet)
                continue;

            const auto& statev = facet->Get_StateVarRef();
            FacetTri tri;
            bool used_exact = false;

            const int ftag = facet->GetTag();
            if (ftag > kLDPMFacetTagBase) {
                const size_t fidx = static_cast<size_t>(ftag - (kLDPMFacetTagBase + 1));
                if (have_tet_nodes && fidx < g_ldpm_facet_vertices.size()) {
                    const auto& pref = g_ldpm_facet_vertices[fidx];
                    ChVector3d p0, p1, p2;
                    const bool ok0 = map_ref_to_deformed(tet_x0, tet_x, pref[0], p0);
                    const bool ok1 = map_ref_to_deformed(tet_x0, tet_x, pref[1], p1);
                    const bool ok2 = map_ref_to_deformed(tet_x0, tet_x, pref[2], p2);
                    if (ok0 && ok1 && ok2) {
                        tri.p0 = p0;
                        tri.p1 = p1;
                        tri.p2 = p2;
                        used_exact = true;
                    }
                }
            }

            if (!used_exact) {
                // Fallback: In multidomain runs, facet-vertex auxiliary nodes are not serialized.
                // Reconstruct a small triangle from center/frame/area if exact facet vertices are unavailable.
                const ChVector3d c = facet->Get_center();
                const ChMatrix33<double> F = facet->Get_facetFrame();
                ChVector3d u(F(1, 0), F(1, 1), F(1, 2));
                ChVector3d v(F(2, 0), F(2, 1), F(2, 2));

                if (u.Length() < 1e-12 || v.Length() < 1e-12) {
                    u = ChVector3d(1, 0, 0);
                    v = ChVector3d(0, 1, 0);
                } else {
                    u = u.GetNormalized();
                    v = v.GetNormalized();
                }

                const double area = std::max(facet->Get_area(), 1e-20);
                const double r = std::sqrt(4.0 * area / (3.0 * std::sqrt(3.0)));
                tri.p0 = c + r * u;
                tri.p1 = c + r * (-0.5 * u + (0.8660254037844386) * v);
                tri.p2 = c + r * (-0.5 * u - (0.8660254037844386) * v);
            }

            tri.crack = statev(11);
            facet_tris.push_back(tri);
        }
    }

    if (facet_tris.empty())
        return false;

    std::ofstream out(vtk_path, std::ios::out);
    if (!out.is_open())
        return false;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << (facet_tris.size() * 3) << "\" NumberOfCells=\"" << facet_tris.size()
        << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& tri : facet_tris) {
        out << "          " << tri.p0.x() << " " << tri.p0.y() << " " << tri.p0.z() << "\n";
        out << "          " << tri.p1.x() << " " << tri.p1.y() << " " << tri.p1.z() << "\n";
        out << "          " << tri.p2.x() << " " << tri.p2.y() << " " << tri.p2.z() << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < facet_tris.size(); ++i) {
        const int i0 = static_cast<int>(3 * i + 0);
        const int i1 = static_cast<int>(3 * i + 1);
        const int i2 = static_cast<int>(3 * i + 2);
        out << "          " << i0 << " " << i1 << " " << i2 << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (size_t i = 0; i < facet_tris.size(); ++i) {
        offset += 3;
        out << "          " << offset << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < facet_tris.size(); ++i) {
        out << "          5\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
    out << "      <CellData Scalars=\"Crack\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"Crack\" format=\"ascii\">\n";
    for (const auto& tri : facet_tris) {
        out << "          " << tri.crack << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </CellData>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    return true;
}

int WriteConcreteSnapshotVTK(ChSystem& sys,
                             const std::filesystem::path& out_dir,
                             int snapshot_id,
                             int step,
                             double time_now) {
    int written_files = 0;
    int mesh_id = 0;
    const long long time_us = static_cast<long long>(std::llround(time_now * 1.0e6));
    for (auto& mesh : sys.GetMeshes()) {
        if (!MeshHasLDPMElements(mesh))
            continue;

        std::ostringstream fname;
        fname << "concrete_snapshot" << snapshot_id << "_step" << step << "_t" << time_us << "us_mesh" << mesh_id
              << ".vtu";
        if (WriteLDPMMeshVTK(mesh, out_dir / fname.str()))
            ++written_files;
        ++mesh_id;
    }
    return written_files;
}

int WriteCrackSnapshotVTK(ChSystem& sys,
                          const std::filesystem::path& out_dir,
                          int snapshot_id,
                          int step,
                          double time_now) {
    int written_files = 0;
    int mesh_id = 0;
    const long long time_us = static_cast<long long>(std::llround(time_now * 1.0e6));
    for (auto& mesh : sys.GetMeshes()) {
        if (!MeshHasLDPMElements(mesh))
            continue;

        std::ostringstream fname;
        fname << "crack_snapshot" << snapshot_id << "_step" << step << "_t" << time_us << "us_mesh" << mesh_id
              << ".vtu";
        if (WriteLDPMCrackVTK(mesh, out_dir / fname.str()))
            ++written_files;
        ++mesh_id;
    }
    return written_files;
}

void ComputeLocalEnergies(ChSystem& sys, double& internal_work, double& kinetic_energy) {
    internal_work = 0.0;
    kinetic_energy = 0.0;

    for (auto& mesh : sys.GetMeshes()) {
        for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
            auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie));
            if (!elem)
                continue;

            for (auto& facet : elem->GetSection()) {
                const auto& statev = facet->Get_StateVarRef();
                internal_work += statev(10);
            }

            ChVectorN<double, 24> V;
            V.setZero();
            for (int j = 0; j < 4; ++j) {
                auto node = elem->GetTetrahedronNode(j);
                V((j + 1) * 6 - 6) = node->GetPosDt().x();
                V((j + 1) * 6 - 5) = node->GetPosDt().y();
                V((j + 1) * 6 - 4) = node->GetPosDt().z();
                V((j + 1) * 6 - 3) = node->GetAngVelLocal().x();
                V((j + 1) * 6 - 2) = node->GetAngVelLocal().y();
                V((j + 1) * 6 - 1) = node->GetAngVelLocal().z();
            }

            ChMatrixNM<double, 24, 24> M;
            M.setZero();
            elem->ComputeMmatrixGlobal(M);
            kinetic_energy += 0.5 * V.transpose() * M * V;
        }
    }
}

std::vector<std::string> CollectSnapshotFiles(const std::filesystem::path& rank_dir,
                                              const std::string& file_prefix,
                                              int snapshot_id) {
    auto parse_snapshot_file_meta = [](const std::string& file_name, int& step, long long& time_us, int& mesh_id) {
        const size_t step_pos = file_name.find("_step");
        const size_t t_pos = file_name.find("_t", step_pos);
        const size_t us_mesh_pos = file_name.find("us_mesh", t_pos);
        const size_t dot_pos = file_name.rfind(".vtu");
        if (step_pos == std::string::npos || t_pos == std::string::npos || us_mesh_pos == std::string::npos ||
            dot_pos == std::string::npos) {
            return false;
        }

        const size_t step_begin = step_pos + std::string("_step").size();
        const size_t t_begin = t_pos + std::string("_t").size();
        const size_t mesh_begin = us_mesh_pos + std::string("us_mesh").size();
        if (step_begin >= t_pos || t_begin >= us_mesh_pos || mesh_begin >= dot_pos) {
            return false;
        }

        try {
            step = std::stoi(file_name.substr(step_begin, t_pos - step_begin));
            time_us = std::stoll(file_name.substr(t_begin, us_mesh_pos - t_begin));
            mesh_id = std::stoi(file_name.substr(mesh_begin, dot_pos - mesh_begin));
        } catch (...) {
            return false;
        }
        return true;
    };

    struct SnapshotPick {
        int step = -1;
        long long time_us = -1;
        std::string file;
    };

    std::vector<std::string> files;
    if (!std::filesystem::exists(rank_dir))
        return files;

    const std::string wanted_prefix = file_prefix + std::to_string(snapshot_id) + "_";
    std::unordered_map<int, SnapshotPick> best_by_mesh;
    for (const auto& entry : std::filesystem::directory_iterator(rank_dir)) {
        if (!entry.is_regular_file())
            continue;
        const std::string candidate = entry.path().filename().string();
        if (candidate.rfind(wanted_prefix, 0) != 0)
            continue;
        if (entry.path().extension() != ".vtu")
            continue;

        int step = -1;
        long long time_us = -1;
        int mesh_id = -1;
        if (!parse_snapshot_file_meta(candidate, step, time_us, mesh_id))
            continue;

        auto it = best_by_mesh.find(mesh_id);
        const bool replace =
            (it == best_by_mesh.end()) || (step > it->second.step) ||
            (step == it->second.step && time_us > it->second.time_us) ||
            (step == it->second.step && time_us == it->second.time_us && candidate > it->second.file);
        if (replace) {
            best_by_mesh[mesh_id] = SnapshotPick{step, time_us, candidate};
        }
    }

    std::vector<std::pair<int, std::string>> sorted_by_mesh;
    sorted_by_mesh.reserve(best_by_mesh.size());
    for (const auto& kv : best_by_mesh) {
        sorted_by_mesh.emplace_back(kv.first, kv.second.file);
    }
    std::sort(sorted_by_mesh.begin(), sorted_by_mesh.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    files.reserve(sorted_by_mesh.size());
    for (const auto& kv : sorted_by_mesh) {
        files.push_back(kv.second);
    }

    return files;
}

void WriteRankedSnapshotPVD(const std::filesystem::path& binary_dir,
                            const std::string& rank_dir_prefix,
                            const std::string& file_prefix,
                            const std::vector<double>& snapshot_times,
                            const std::filesystem::path& pvd_path,
                            int mpi_total_ranks) {
    struct Entry {
        double time;
        int part;
        std::string file;
    };

    std::vector<Entry> entries;
    for (size_t isnap = 0; isnap < snapshot_times.size(); ++isnap) {
        for (int irank = 0; irank < mpi_total_ranks; ++irank) {
            std::ostringstream rank_tag;
            rank_tag << std::setw(3) << std::setfill('0') << irank;
            const std::string rank_dir_name = rank_dir_prefix + rank_tag.str();
            const std::filesystem::path rank_dir = binary_dir / rank_dir_name;
            if (!std::filesystem::exists(rank_dir))
                continue;

            auto files = CollectSnapshotFiles(rank_dir, file_prefix, static_cast<int>(isnap + 1));
            for (size_t imesh = 0; imesh < files.size(); ++imesh) {
                entries.push_back({snapshot_times[isnap], static_cast<int>(1000 * irank + imesh),
                                   rank_dir_name + "/" + files[imesh]});
            }
        }
    }

    std::ofstream out(pvd_path, std::ios::out);
    if (!out.is_open())
        return;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (const auto& entry : entries) {
        out << "    <DataSet timestep=\"" << entry.time << "\" group=\"\" part=\"" << entry.part << "\" file=\""
            << entry.file << "\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    const RunSettings cfg = ParseArgs(argc, argv);
    const std::filesystem::path example_dir(LDPM_EXAMPLE_DIR);
    const std::filesystem::path data_dir = example_dir / "LDPMgeo000Box000";

    ChDomainManagerMPI domain_manager(&argc, &argv);
    const int mpi_rank = domain_manager.GetMPIrank();
    const int mpi_total_ranks = domain_manager.GetMPItotranks();

    if (mpi_total_ranks < 2) {
        if (mpi_rank == 0) {
            std::cerr << "Error: this demo requires at least 2 MPI ranks (1 worker + 1 master)." << std::endl;
        }
        return 1;
    }

    const int worker_domains = mpi_total_ranks - 1;
    const char* no_master_mig_env = std::getenv("CH_MDOM_NO_MASTER_MIGRATION");
    const bool no_master_migration =
        no_master_mig_env && std::string(no_master_mig_env) == "1";

    if (cfg.vtk_stride > 0) {
        const bool facet_lookup_ok = LoadFacetVertexLookup(data_dir);
        if (!facet_lookup_ok && mpi_rank == 0) {
            std::cerr << "Warning: could not load LDPM facet-vertex lookup; crack output will use fallback triangles."
                      << std::endl;
        }
    }

    domain_manager.verbose_partition = cfg.verbose_partition;
    domain_manager.verbose_serialization = cfg.verbose_serialization;
    domain_manager.verbose_variable_updates = cfg.verbose_variable_updates;
    domain_manager.serializer_type = cfg.serializer_format;

    constexpr double specimen_min_x = 0.0;
    constexpr double specimen_max_x = 100.0;
    constexpr double specimen_min_y = 0.0;
    constexpr double specimen_max_y = 100.0;
    constexpr double specimen_min_z = 0.0;
    constexpr double specimen_max_z = 200.0;

    std::array<int, 3> grid_divisions = {1, 1, 1};
    if (cfg.use_grid_decomposition) {
        grid_divisions = {cfg.grid_x, cfg.grid_y, cfg.grid_z};
        if (grid_divisions[0] * grid_divisions[1] * grid_divisions[2] != worker_domains) {
            if (mpi_rank == 0) {
                std::cerr << "Error: grid dimensions " << grid_divisions[0] << "x" << grid_divisions[1] << "x"
                          << grid_divisions[2] << " do not equal the " << worker_domains << " worker ranks."
                          << std::endl;
            }
            return 1;
        }
    } else if (cfg.slice_axis == ChAxis::X) {
        grid_divisions[0] = worker_domains;
    } else if (cfg.slice_axis == ChAxis::Y) {
        grid_divisions[1] = worker_domains;
    } else {
        grid_divisions[2] = worker_domains;
    }

    double slice_min = specimen_min_z;
    double slice_max = specimen_max_z;
    if (cfg.slice_axis == ChAxis::X || cfg.slice_axis == ChAxis::Y) {
        slice_max = specimen_max_x;
    }

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    const int master_rank = worker_domains;
    const bool load_local_mesh =
        no_master_migration ? (mpi_rank != master_rank) : (mpi_rank == master_rank);
    std::shared_ptr<chrono::multidomain::ChDomain> local_domain;
    if (cfg.use_grid_decomposition) {
        ChDomainBuilderGrid domain_builder(grid_divisions[0], specimen_min_x, specimen_max_x, grid_divisions[1],
                                           specimen_min_y, specimen_max_y, grid_divisions[2], specimen_min_z,
                                           specimen_max_z, true);
        if (mpi_rank == master_rank) {
            local_domain = domain_builder.BuildMasterDomain(&sys);
        } else {
            local_domain = domain_builder.BuildDomain(&sys, mpi_rank);
            if (cfg.use_complete_worker_graph) {
                for (int peer = 0; peer < worker_domains; ++peer) {
                    if (peer == mpi_rank || local_domain->GetInterfaces().count(peer))
                        continue;
                    int i = 0;
                    int j = 0;
                    int k = 0;
                    domain_builder.FromRankToIJKindexes(peer, i, j, k);
                    auto grid_bound = [](int index, int divisions, double minimum, double maximum) {
                        const double width = (maximum - minimum) / static_cast<double>(divisions);
                        const double lower = (index == 0) ? -1.0e34 : minimum + index * width;
                        const double upper = (index == divisions - 1) ? 1.0e34 : minimum + (index + 1) * width;
                        return std::pair<double, double>(lower, upper);
                    };
                    const auto xb = grid_bound(i, grid_divisions[0], specimen_min_x, specimen_max_x);
                    const auto yb = grid_bound(j, grid_divisions[1], specimen_min_y, specimen_max_y);
                    const auto zb = grid_bound(k, grid_divisions[2], specimen_min_z, specimen_max_z);
                    auto peer_domain = chrono_types::make_shared<ChDomainBox>(
                        nullptr, peer,
                        ChAABB(ChVector3d(xb.first, yb.first, zb.first), ChVector3d(xb.second, yb.second, zb.second)));
                    auto& interface = local_domain->GetInterfaces()[peer];
                    interface.side_IN = local_domain.get();
                    interface.side_OUT = peer_domain;
                }
            }
        }
    } else {
        ChDomainBuilderSlices domain_builder(worker_domains, slice_min, slice_max, cfg.slice_axis, true);
        if (mpi_rank == master_rank) {
            local_domain = domain_builder.BuildMasterDomain(&sys);
        } else {
            local_domain = domain_builder.BuildDomain(&sys, mpi_rank);
            if (cfg.use_complete_worker_graph) {
                const double width = (slice_max - slice_min) / static_cast<double>(worker_domains);
                for (int peer = 0; peer < worker_domains; ++peer) {
                    if (peer == mpi_rank || local_domain->GetInterfaces().count(peer))
                        continue;
                    const double lower = (peer == 0) ? -1.0e34 : slice_min + peer * width;
                    const double upper =
                        (peer == worker_domains - 1) ? 1.0e34 : slice_min + (peer + 1) * width;
                    auto peer_domain =
                        chrono_types::make_shared<ChDomainSlice>(nullptr, peer, lower, upper, cfg.slice_axis);
                    auto& interface = local_domain->GetInterfaces()[peer];
                    interface.side_IN = local_domain.get();
                    interface.side_OUT = peer_domain;
                }
            }
        }
    }
    domain_manager.SetDomain(local_domain);

    auto explicit_timestepper = chrono_types::make_shared<ChTimestepperLeapfrog>(&sys);
    if (cfg.use_diagonal_lumping) {
        explicit_timestepper->SetDiagonalLumpingON(cfg.constraint_penalty, cfg.constraint_damping);
    } else {
        explicit_timestepper->SetDiagonalLumpingOFF();
    }
    sys.SetTimestepper(explicit_timestepper);

    auto lumped_solver = chrono_types::make_shared<ChSolverLumpedMultidomain>();
    sys.SetSolver(lumped_solver);

    double top_surface_area = 0.0;
    double specimen_height = 0.0;

    if (load_local_mesh) {
        std::string ldpm_data_path = data_dir.string() + "/";
        std::string ldpm_geo_name = "LDPMgeo000";

        auto mesh = chrono_types::make_shared<ChMesh>();
        auto mat = chrono_types::make_shared<ChMaterialVECT>();
        ApplyLDPMMaterialParameters(mat);

        ChBuilderLDPM builder;
        builder.read_LDPM_info(mesh, mat, ldpm_data_path, ldpm_geo_name);
        mesh->SetAutomaticGravity(false);
        sys.Add(mesh);

        auto rp_top = FindClosestNodeAtPlane(mesh, 200.0, cfg.z_tol_ref, cfg.ref_x, cfg.ref_y);
        auto rp_bottom = FindClosestNodeAtPlane(mesh, 0.0, cfg.z_tol_ref, cfg.ref_x, cfg.ref_y);

        if (!rp_top || !rp_bottom) {
            std::cerr << "Failed to identify top/bottom reference nodes in LDPM mesh." << std::endl;
            return 2;
        }

        auto top_nodes = CollectSurfaceNodes(mesh, 200.0, cfg.z_tol_surface, rp_top);
        auto bottom_nodes = CollectSurfaceNodes(mesh, 0.0, cfg.z_tol_surface, rp_bottom);

        if (top_nodes.empty() || bottom_nodes.empty()) {
            std::cerr << "Failed to identify top/bottom surface nodes for boundary conditions." << std::endl;
            return 3;
        }

        double top_min_x = std::numeric_limits<double>::max();
        double top_max_x = std::numeric_limits<double>::lowest();
        double top_min_y = std::numeric_limits<double>::max();
        double top_max_y = std::numeric_limits<double>::lowest();
        for (const auto& node : top_nodes) {
            const ChVector3d p = node->GetPos();
            top_min_x = std::min(top_min_x, p.x());
            top_max_x = std::max(top_max_x, p.x());
            top_min_y = std::min(top_min_y, p.y());
            top_max_y = std::max(top_max_y, p.y());
        }

        top_surface_area = std::max(0.0, top_max_x - top_min_x) * std::max(0.0, top_max_y - top_min_y);
        specimen_height = std::abs(rp_top->GetX0().GetPos().z() - rp_bottom->GetX0().GetPos().z());

        if (top_surface_area <= 0.0 || specimen_height <= 0.0) {
            std::cerr << "Invalid geometry metrics (top area or specimen height)." << std::endl;
            return 4;
        }

        // First assign generic unique tags to all entities, then overwrite boundary-node tags
        // with reserved ranges used for deterministic boundary selection/reaction reduction.
        SetUniqueTags(sys);
        TagLDPMSharedObjects(mesh, mat, 1000000);

        for (size_t i = 0; i < top_nodes.size(); ++i) {
            top_nodes[i]->SetTag(kTagTopNodeMin + static_cast<int>(i));
        }
        for (size_t i = 0; i < bottom_nodes.size(); ++i) {
            bottom_nodes[i]->SetTag(kTagBottomNodeMin + static_cast<int>(i));
        }
        // Keep reference nodes in dedicated trailing boundary tags so ownership and BC
        // reduction include them (historically this yields 134 owned top nodes).
        rp_top->SetTag(kTagTopNodeMin + static_cast<int>(top_nodes.size()));
        rp_bottom->SetTag(kTagBottomNodeMin + static_cast<int>(bottom_nodes.size()));

        // NOTE: explicit top/bottom tie links and motor links are disabled for MPI multidomain runs
        // because ChLinkMateGeneric ownership during domain migration can segfault in initialization.

        std::cout << "LDPM multidomain setup complete on master rank."
                  << " worker_domains=" << worker_domains
                  << " top_nodes=" << top_nodes.size()
                  << " bottom_nodes=" << bottom_nodes.size()
                  << " top_surface_area=" << top_surface_area
                  << " specimen_height=" << specimen_height
                  << " dt=" << cfg.dt
                  << " t_end=" << cfg.t_end
                  << " disp_rate=" << cfg.disp_rate
                  << " poly_a=" << cfg.poly_a
                  << " switch_time=" << cfg.switch_time
                  << " load_scale=" << cfg.load_scale
                  << " decomposition=" << (cfg.use_grid_decomposition ? "grid" : "slices")
                  << " grid=" << grid_divisions[0] << "x" << grid_divisions[1] << "x" << grid_divisions[2]
                  << " worker_graph=" << (cfg.use_complete_worker_graph ? "complete" : "neighbor")
                  << " partition_update_each_step=" << (cfg.do_partition_update_each_step ? "on" : "off")
                  << " slice_axis="
                  << (cfg.slice_axis == ChAxis::X ? "x" : (cfg.slice_axis == ChAxis::Y ? "y" : "z"))
                  << " diagonal_lumping=" << (cfg.use_diagonal_lumping ? "on" : "off")
                  << " Ck=" << cfg.constraint_penalty
                  << " Cr=" << cfg.constraint_damping
                  << " zslice_surface_owner=" << (cfg.use_zslice_surface_owner ? "on" : "off")
                  << " tag_reaction_selector=" << (cfg.use_tag_reaction_selector ? "on" : "off")
                  << " fix_top_rotations=" << (cfg.fix_top_rotations ? "on" : "off")
                  << " boundary_mode=" << BoundaryModeName(cfg.boundary_mode)
                  << " clean_output=" << (cfg.clean_output_dirs ? "on" : "off")
                  << " master_after_init=" << (cfg.disable_master_after_init ? "disabled" : "enabled")
                  << std::endl;
    }

    double geom_bcast[2] = {top_surface_area, specimen_height};
    const int metadata_root_rank = no_master_migration ? 0 : master_rank;
    MPI_Bcast(geom_bcast, 2, MPI_DOUBLE, metadata_root_rank, MPI_COMM_WORLD);
    top_surface_area = geom_bcast[0];
    specimen_height = geom_bcast[1];

    if (top_surface_area <= 0.0 || specimen_height <= 0.0) {
        if (mpi_rank == 0) {
            std::cerr << "Invalid metadata received after master setup; aborting run." << std::endl;
        }
        return 5;
    }

    if (no_master_migration) {
        domain_manager.master_domain_enabled = false;
    }
    domain_manager.DoDomainInitialize(mpi_rank);
    if (cfg.disable_master_after_init) {
        domain_manager.master_domain_enabled = false;
    }
    if (cfg.do_partition_update_once && worker_domains > 1) {
        // Run one worker-only repartition exchange to finalize interfaces.
        // Use cleanup so outside elements/nodes from the initial master spill are pruned.
        domain_manager.DoDomainPartitionUpdate(mpi_rank, true);
    }

    // Ensure variable offsets are rebuilt after any repartition exchange before
    // we form owned-node sets and interface diagnostics.
    const char* skip_setup_env = std::getenv("CH_MDOM_SKIP_SETUP");
    const bool skip_setup = skip_setup_env && std::string(skip_setup_env) == "1";
    if (!skip_setup) {
        sys.Setup();
    }

    auto local_motor = FindMotorByTag(sys, kTagMotor);
    auto local_constr_bot = FindMateByTag(sys, kTagConstrBot);
    const int no_owner_rank = mpi_total_ranks + 1;
    int motor_owner_rank_local = local_motor ? mpi_rank : no_owner_rank;
    int constr_owner_rank_local = local_constr_bot ? mpi_rank : no_owner_rank;
    int motor_owner_rank = no_owner_rank;
    int constr_owner_rank = no_owner_rank;
    MPI_Allreduce(&motor_owner_rank_local, &motor_owner_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&constr_owner_rank_local, &constr_owner_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    const bool use_motor_drive = (motor_owner_rank < no_owner_rank);

    if (mpi_rank == 0) {
        std::cout << "Calibration boundary ownership: motor_owner_rank=" << motor_owner_rank
                  << " constr_bot_owner_rank=" << constr_owner_rank
                  << " motor_drive=" << (use_motor_drive ? "on" : "off") << std::endl;
    }

    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> local_top_all;
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> local_bottom_all;
    const double axis_min = slice_min;
    const double axis_max = slice_max;
    const double slice_w = (axis_max - axis_min) / static_cast<double>(worker_domains);
    auto coord_on_axis = [&](const ChVector3d& p) {
        if (cfg.slice_axis == ChAxis::X)
            return p.x();
        if (cfg.slice_axis == ChAxis::Y)
            return p.y();
        return p.z();
    };
    auto owner_from_axis_coord = [&](double coord) {
        const double c_shift = std::clamp(coord - axis_min - 1.0e-12, 0.0, (axis_max - axis_min) - 1.0e-12);
        const int owner = static_cast<int>(std::floor(c_shift / slice_w));
        return std::clamp(owner, 0, worker_domains - 1);
    };
    auto grid_index = [](double coord, double minimum, double maximum, int divisions) {
        const double shifted = std::clamp(coord - minimum - 1.0e-12, 0.0, (maximum - minimum) - 1.0e-12);
        const double width = (maximum - minimum) / static_cast<double>(divisions);
        return std::clamp(static_cast<int>(std::floor(shifted / width)), 0, divisions - 1);
    };
    auto owner_from_position = [&](const ChVector3d& position) {
        if (!cfg.use_grid_decomposition)
            return owner_from_axis_coord(coord_on_axis(position));
        const int i = grid_index(position.x(), specimen_min_x, specimen_max_x, grid_divisions[0]);
        const int j = grid_index(position.y(), specimen_min_y, specimen_max_y, grid_divisions[1]);
        const int k = grid_index(position.z(), specimen_min_z, specimen_max_z, grid_divisions[2]);
        return i + j * grid_divisions[0] + k * grid_divisions[0] * grid_divisions[1];
    };

    for (auto& mesh : sys.GetMeshes()) {
        for (auto& node_base : mesh->GetNodes()) {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
            if (!node)
                continue;
            const int tag = node->GetTag();
            if (tag >= kTagTopNodeMin && tag <= kTagTopNodeMax)
                local_top_all.push_back(node);
            if (tag >= kTagBottomNodeMin && tag <= kTagBottomNodeMax)
                local_bottom_all.push_back(node);
        }
    }

    auto get_node_weight = [&domain_manager](ChNodeFEAxyzrot* node) {
        if (!node || !domain_manager.domain)
            return 0.0;
        const auto& wv = domain_manager.domain->CoordWeightsWv();
        const int off = node->Variables().GetOffset();
        if (off < 0 || off >= static_cast<int>(wv.size()))
            return 0.0;
        return wv(off);
    };

    std::vector<ChNodeFEAxyzrot*> local_top_bc_nodes;
    std::vector<ChNodeFEAxyzrot*> local_bottom_bc_nodes;
    std::vector<ChNodeFEAxyzrot*> local_top_nodes;
    std::unordered_set<ChNodeFEAxyzrot*> local_top_node_ptrs;
    std::unordered_set<ChNodeFEAxyzrot*> local_bottom_node_ptrs;
    std::unordered_set<int> local_top_node_tags;
    std::unordered_set<int> local_bottom_node_tags;
    std::unordered_set<ChNodeFEAxyzrot*> seen_top_bc_ptrs;
    std::unordered_set<ChNodeFEAxyzrot*> seen_bottom_bc_ptrs;
    seen_top_bc_ptrs.reserve(local_top_all.size() * 2 + 1);
    seen_bottom_bc_ptrs.reserve(local_bottom_all.size() * 2 + 1);
    for (const auto& node : local_top_all) {
        if (!no_master_migration && !node->Variables().IsActive())
            continue;
        if (seen_top_bc_ptrs.insert(node.get()).second) {
            // Apply prescribed top displacement to every local copy of tagged top nodes.
            local_top_bc_nodes.push_back(node.get());
        }
        bool keep_owned = true;
        if (cfg.use_zslice_surface_owner) {
            const int owner_rank = owner_from_position(node->GetX0().GetPos());
            keep_owned = (mpi_rank == owner_rank);
        }
        if (!cfg.use_zslice_surface_owner) {
            const double w = get_node_weight(node.get());
            keep_owned = (w > 0.0);
        }
        if (!keep_owned) {
            continue;
        }
        local_top_nodes.push_back(node.get());
        local_top_node_ptrs.insert(node.get());
        local_top_node_tags.insert(node->GetTag());
    }
    for (const auto& node : local_bottom_all) {
        if (!no_master_migration && !node->Variables().IsActive())
            continue;
        node->SetFixed((cfg.boundary_mode == BoundaryMode::FIX) && !use_motor_drive);
        if (seen_bottom_bc_ptrs.insert(node.get()).second) {
            // Apply bottom boundary prescription to every local copy of tagged bottom nodes.
            local_bottom_bc_nodes.push_back(node.get());
        }

        bool keep_owned = true;
        if (cfg.use_zslice_surface_owner) {
            const int owner_rank = owner_from_position(node->GetX0().GetPos());
            keep_owned = (mpi_rank == owner_rank);
        }
        if (!cfg.use_zslice_surface_owner) {
            const double w = get_node_weight(node.get());
            keep_owned = (w > 0.0);
        }
        if (!keep_owned) {
            continue;
        }
        local_bottom_node_ptrs.insert(node.get());
        local_bottom_node_tags.insert(node->GetTag());
    }

    auto enforce_unique_owner_by_tag = [&](std::unordered_set<int>& local_tags) {
        const int local_count = static_cast<int>(local_tags.size());
        std::vector<int> recv_counts(mpi_total_ranks, 0);
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> displs(mpi_total_ranks, 0);
        int total_entries = 0;
        for (int r = 0; r < mpi_total_ranks; ++r) {
            displs[r] = total_entries;
            total_entries += recv_counts[r];
        }

        std::vector<int> send_tags(local_tags.begin(), local_tags.end());
        std::vector<int> all_tags(total_entries, -1);
        MPI_Allgatherv(send_tags.data(), local_count, MPI_INT, all_tags.data(), recv_counts.data(), displs.data(),
                       MPI_INT, MPI_COMM_WORLD);

        std::unordered_map<int, int> owner_rank_by_tag;
        owner_rank_by_tag.reserve(static_cast<size_t>(total_entries));
        for (int r = 0; r < mpi_total_ranks; ++r) {
            for (int i = 0; i < recv_counts[r]; ++i) {
                const int tag = all_tags[displs[r] + i];
                auto it = owner_rank_by_tag.find(tag);
                if (it == owner_rank_by_tag.end() || r < it->second) {
                    owner_rank_by_tag[tag] = r;
                }
            }
        }

        std::vector<int> tags_to_remove;
        tags_to_remove.reserve(local_tags.size());
        for (int tag : local_tags) {
            auto it = owner_rank_by_tag.find(tag);
            if (it != owner_rank_by_tag.end() && it->second != mpi_rank) {
                tags_to_remove.push_back(tag);
            }
        }
        for (int tag : tags_to_remove) {
            local_tags.erase(tag);
        }

        int unique_tag_count = static_cast<int>(owner_rank_by_tag.size());
        int global_occurrences = 0;
        MPI_Allreduce(&local_count, &global_occurrences, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        int duplicate_owners = global_occurrences - unique_tag_count;
        return duplicate_owners;
    };

    int duplicate_top_owners = 0;
    int duplicate_bottom_owners = 0;
    if (!no_master_migration) {
        duplicate_top_owners = enforce_unique_owner_by_tag(local_top_node_tags);
        duplicate_bottom_owners = enforce_unique_owner_by_tag(local_bottom_node_tags);
    }

    local_top_nodes.erase(std::remove_if(local_top_nodes.begin(), local_top_nodes.end(),
                                         [&](ChNodeFEAxyzrot* n) { return local_top_node_tags.count(n->GetTag()) == 0; }),
                          local_top_nodes.end());
    local_top_node_ptrs.clear();
    for (auto* n : local_top_nodes) {
        local_top_node_ptrs.insert(n);
    }
    for (auto it = local_bottom_node_ptrs.begin(); it != local_bottom_node_ptrs.end();) {
        if (local_bottom_node_tags.count((*it)->GetTag()) == 0) {
            it = local_bottom_node_ptrs.erase(it);
        } else {
            ++it;
        }
    }

    const int local_top_count = static_cast<int>(local_top_nodes.size());
    int global_top_count = 0;
    MPI_Allreduce(&local_top_count, &global_top_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_top_count <= 0 && cfg.t_end > 0) {
        if (mpi_rank == 0) {
            std::cerr << "Could not determine owned top-surface nodes after multidomain initialization." << std::endl;
        }
        return 6;
    }

    std::vector<int> gathered_top_counts;
    if (mpi_rank == 0)
        gathered_top_counts.resize(mpi_total_ranks, 0);
    MPI_Gather(&local_top_count, 1, MPI_INT, mpi_rank == 0 ? gathered_top_counts.data() : nullptr, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::cout << "Owned top-node counts per rank:";
        for (int r = 0; r < mpi_total_ranks; ++r)
            std::cout << " r" << r << "=" << gathered_top_counts[r];
        std::cout << " global=" << global_top_count << " dedup_top=" << duplicate_top_owners
                  << " dedup_bottom=" << duplicate_bottom_owners << std::endl;
    }

    auto report_global_tag_multiplicity = [&](const std::unordered_set<int>& local_tags, const char* label) {
        const int local_count = static_cast<int>(local_tags.size());
        std::vector<int> recv_counts(mpi_total_ranks, 0);
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> displs(mpi_total_ranks, 0);
        int total_entries = 0;
        for (int r = 0; r < mpi_total_ranks; ++r) {
            displs[r] = total_entries;
            total_entries += recv_counts[r];
        }

        std::vector<int> send_tags(local_tags.begin(), local_tags.end());
        std::vector<int> all_tags(total_entries, -1);
        MPI_Allgatherv(send_tags.data(), local_count, MPI_INT, all_tags.data(), recv_counts.data(), displs.data(),
                       MPI_INT, MPI_COMM_WORLD);

        std::unordered_map<int, int> multiplicity;
        multiplicity.reserve(static_cast<size_t>(total_entries) * 2 + 1);
        for (int tag : all_tags) {
            if (tag >= 0)
                multiplicity[tag] += 1;
        }

        int duplicated_tags = 0;
        int duplicate_owners = 0;
        for (const auto& kv : multiplicity) {
            if (kv.second > 1) {
                ++duplicated_tags;
                duplicate_owners += (kv.second - 1);
            }
        }

        if (mpi_rank == 0) {
            std::cout << "[mdom-tag-audit] " << label << " local_total_entries=" << total_entries
                      << " unique_tags=" << multiplicity.size() << " duplicated_tags=" << duplicated_tags
                      << " duplicate_owners=" << duplicate_owners << std::endl;
        }
    };

    auto report_shared_tag_coverage = [&](const std::unordered_set<int>& local_all_tags,
                                          const std::unordered_set<int>& local_shared_tags,
                                          const char* label) {
        const int local_all_count = static_cast<int>(local_all_tags.size());
        const int local_shared_count = static_cast<int>(local_shared_tags.size());

        std::vector<int> recv_all_counts(mpi_total_ranks, 0);
        std::vector<int> recv_shared_counts(mpi_total_ranks, 0);
        MPI_Allgather(&local_all_count, 1, MPI_INT, recv_all_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&local_shared_count, 1, MPI_INT, recv_shared_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> all_displs(mpi_total_ranks, 0);
        std::vector<int> shared_displs(mpi_total_ranks, 0);
        int total_all_entries = 0;
        int total_shared_entries = 0;
        for (int r = 0; r < mpi_total_ranks; ++r) {
            all_displs[r] = total_all_entries;
            shared_displs[r] = total_shared_entries;
            total_all_entries += recv_all_counts[r];
            total_shared_entries += recv_shared_counts[r];
        }

        std::vector<int> send_all_tags(local_all_tags.begin(), local_all_tags.end());
        std::vector<int> send_shared_tags(local_shared_tags.begin(), local_shared_tags.end());
        std::vector<int> all_tags(total_all_entries, -1);
        std::vector<int> shared_tags(total_shared_entries, -1);

        MPI_Allgatherv(send_all_tags.data(), local_all_count, MPI_INT, all_tags.data(), recv_all_counts.data(),
                       all_displs.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(send_shared_tags.data(), local_shared_count, MPI_INT, shared_tags.data(), recv_shared_counts.data(),
                       shared_displs.data(), MPI_INT, MPI_COMM_WORLD);

        std::unordered_map<int, int> owner_count;
        std::unordered_map<int, int> shared_owner_count;
        owner_count.reserve(static_cast<size_t>(total_all_entries) * 2 + 1);
        shared_owner_count.reserve(static_cast<size_t>(total_shared_entries) * 2 + 1);

        for (int tag : all_tags) {
            if (tag >= 0)
                owner_count[tag] += 1;
        }
        for (int tag : shared_tags) {
            if (tag >= 0)
                shared_owner_count[tag] += 1;
        }

        if (mpi_rank == 0) {
            int duplicated_tags = 0;
            int duplicated_not_fully_shared = 0;
            int duplicated_unshared_owners_total = 0;
            int duplicated_shared_on_all_owners = 0;
            for (const auto& kv : owner_count) {
                const int tag = kv.first;
                const int owners = kv.second;
                if (owners <= 1)
                    continue;
                duplicated_tags += 1;
                const int shared_owners = (shared_owner_count.count(tag) > 0) ? shared_owner_count[tag] : 0;
                if (shared_owners >= owners) {
                    duplicated_shared_on_all_owners += 1;
                } else {
                    duplicated_not_fully_shared += 1;
                    duplicated_unshared_owners_total += (owners - shared_owners);
                }
            }

            std::cout << "[mdom-shared-coverage] " << label
                      << " duplicated_tags=" << duplicated_tags
                      << " duplicated_shared_on_all_owners=" << duplicated_shared_on_all_owners
                      << " duplicated_not_fully_shared=" << duplicated_not_fully_shared
                      << " duplicated_unshared_owners_total=" << duplicated_unshared_owners_total
                      << std::endl;
        }
    };

    std::unordered_set<int> local_node_tags_all;
    std::unordered_set<int> local_element_tags_all;
    std::unordered_set<int> local_shared_node_tags_all;
    std::unordered_map<int, std::shared_ptr<ChNodeBase>> local_nodes_by_tag;
    for (const auto& mesh : sys.GetMeshes()) {
        for (const auto& node : mesh->GetNodes()) {
            if (node) {
                local_node_tags_all.insert(node->GetTag());
                local_nodes_by_tag[node->GetTag()] = node;
            }
        }
        for (const auto& elem : mesh->GetElements()) {
            if (elem)
                local_element_tags_all.insert(elem->GetTag());
        }
    }

    auto enforce_interface_sharing_from_global_duplicates = [&]() {
        if (!domain_manager.domain)
            return;

        const int local_count = static_cast<int>(local_node_tags_all.size());
        std::vector<int> recv_counts(mpi_total_ranks, 0);
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> displs(mpi_total_ranks, 0);
        int total_entries = 0;
        for (int r = 0; r < mpi_total_ranks; ++r) {
            displs[r] = total_entries;
            total_entries += recv_counts[r];
        }

        std::vector<int> send_tags(local_node_tags_all.begin(), local_node_tags_all.end());
        std::vector<int> all_tags(total_entries, -1);
        MPI_Allgatherv(send_tags.data(), local_count, MPI_INT, all_tags.data(), recv_counts.data(), displs.data(),
                       MPI_INT, MPI_COMM_WORLD);

        std::unordered_map<int, std::vector<int>> owners_by_tag;
        owners_by_tag.reserve(static_cast<size_t>(total_entries) * 2 + 1);
        for (int r = 0; r < mpi_total_ranks; ++r) {
            for (int i = 0; i < recv_counts[r]; ++i) {
                const int tag = all_tags[displs[r] + i];
                if (tag >= 0)
                    owners_by_tag[tag].push_back(r);
            }
        }

        auto& interfaces = domain_manager.domain->GetInterfaces();
        long long local_expected_peer_links = 0;
        long long local_missing_peer_links = 0;
        for (const auto& kv : owners_by_tag) {
            const int tag = kv.first;
            const auto& owners = kv.second;
            if (owners.size() < 2 || local_nodes_by_tag.find(tag) == local_nodes_by_tag.end())
                continue;
            for (int peer : owners) {
                if (peer == mpi_rank || peer == master_rank)
                    continue;
                ++local_expected_peer_links;
                auto iface_it = interfaces.find(peer);
                if (iface_it == interfaces.end()) {
                    ++local_missing_peer_links;
                    continue;
                }
                iface_it->second.shared_nodes[tag] = local_nodes_by_tag[tag];
            }
        }

        long long global_expected_peer_links = 0;
        long long global_missing_peer_links = 0;
        MPI_Allreduce(&local_expected_peer_links, &global_expected_peer_links, 1, MPI_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&local_missing_peer_links, &global_missing_peer_links, 1, MPI_LONG_LONG, MPI_SUM,
                      MPI_COMM_WORLD);

        auto make_shared_var_key = [](bool is_node, int tag, int local_var_index, int dof) {
            const std::uint64_t type_bit = is_node ? (1ull << 63) : 0ull;
            const std::uint64_t tag_bits =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(tag)) & 0xFFFFFFFFull) << 31;
            const std::uint64_t local_bits =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(local_var_index)) & 0xFFFFull) << 15;
            const std::uint64_t dof_bits =
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(dof)) & 0x7FFFull;
            return type_bit | tag_bits | local_bits | dof_bits;
        };

        long long local_key_mismatches = 0;
        for (auto& iface_kv : interfaces) {
            const int neigh_rank = iface_kv.first;
            auto& iface = iface_kv.second;
            iface.shared_vars.clear();
            iface.shared_var_keys.clear();
            ChSystemDescriptor temp_descr;

            for (auto& mshitem : iface.shared_items) {
                if (std::dynamic_pointer_cast<fea::ChMesh>(mshitem.second))
                    continue;
                const size_t before = temp_descr.GetVariables().size();
                mshitem.second->InjectVariables(temp_descr);
                const auto variables = temp_descr.GetVariables();
                for (size_t i = before; i < variables.size(); ++i) {
                    iface.shared_var_keys.push_back(
                        make_shared_var_key(false, mshitem.first, static_cast<int>(i - before), variables[i]->GetDOF()));
                }
            }
            for (auto& mshnode : iface.shared_nodes) {
                const size_t before = temp_descr.GetVariables().size();
                mshnode.second->InjectVariables(temp_descr);
                const auto variables = temp_descr.GetVariables();
                for (size_t i = before; i < variables.size(); ++i) {
                    iface.shared_var_keys.push_back(
                        make_shared_var_key(true, mshnode.first, static_cast<int>(i - before), variables[i]->GetDOF()));
                }
            }

            const auto raw_variables = temp_descr.GetVariables();
            std::unordered_set<std::uint64_t> seen_keys;
            std::vector<ChVariables*> unique_variables;
            std::vector<std::uint64_t> unique_keys;
            unique_variables.reserve(raw_variables.size());
            unique_keys.reserve(raw_variables.size());
            for (size_t i = 0; i < raw_variables.size() && i < iface.shared_var_keys.size(); ++i) {
                auto* variable = raw_variables[i];
                const auto key = iface.shared_var_keys[i];
                if (!variable || !seen_keys.insert(key).second)
                    continue;
                unique_variables.push_back(variable);
                unique_keys.push_back(key);
            }
            iface.shared_vars.swap(unique_variables);
            iface.shared_var_keys.swap(unique_keys);
            if (iface.shared_vars.size() != iface.shared_var_keys.size())
                ++local_key_mismatches;
        }

        long long global_key_mismatches = 0;
        MPI_Allreduce(&local_key_mismatches, &global_key_mismatches, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            std::cout << "[mdom-peer-coverage] expected_directed_links=" << global_expected_peer_links
                      << " missing_directed_links=" << global_missing_peer_links
                      << " key_size_mismatched_interfaces=" << global_key_mismatches << std::endl;
        }
        if (global_missing_peer_links != 0 || global_key_mismatches != 0) {
            throw std::runtime_error("incomplete shared-node peer graph or shared-variable key mismatch");
        }

        sys.Setup();
        domain_manager.domain->ComputeSharedCoordsWeights(domain_manager.domain->CoordWeightsWv());
    };
    enforce_interface_sharing_from_global_duplicates();

    if (domain_manager.domain) {
        for (const auto& iface : domain_manager.domain->GetInterfaces()) {
            for (const auto& sn : iface.second.shared_nodes) {
                local_shared_node_tags_all.insert(sn.first);
            }
        }
    }
    report_global_tag_multiplicity(local_node_tags_all, "nodes");
    report_global_tag_multiplicity(local_element_tags_all, "elements");
    report_shared_tag_coverage(local_node_tags_all, local_shared_node_tags_all, "nodes");

    // Elements have unique ownership after partitioning, so summing rho*volume
    // gives a duplicate-safe specimen mass audit without modifying the model.
    double local_specimen_mass = 0.0;
    for (const auto& mesh : sys.GetMeshes()) {
        for (const auto& element : mesh->GetElements()) {
            auto ldpm_element = std::dynamic_pointer_cast<ChElementLDPM>(element);
            if (ldpm_element) {
                local_specimen_mass += ldpm_element->GetVolume() * ldpm_element->GetDensity();
            }
        }
    }
    double global_specimen_mass = 0.0;
    MPI_Allreduce(&local_specimen_mass, &global_specimen_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::cout << "GLOBAL_SPECIMEN_MASS=" << std::setprecision(17) << global_specimen_mass << std::endl;
    }

    const std::filesystem::path binary_dir = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream rank_tag;
    rank_tag << std::setw(3) << std::setfill('0') << mpi_rank;
    const std::filesystem::path out_dir =
        binary_dir / ("out_mdom_ldpm_displacement_w" + std::to_string(worker_domains) + "_rank" + rank_tag.str());
    const std::string out_dir_prefix = "out_mdom_ldpm_displacement_w" + std::to_string(worker_domains) + "_rank";

    std::error_code out_ec;
    if (cfg.clean_output_dirs && std::filesystem::exists(out_dir)) {
        std::filesystem::remove_all(out_dir, out_ec);
        if (out_ec) {
            std::cerr << "Rank " << mpi_rank << " failed cleaning output directory: " << out_ec.message()
                      << std::endl;
            return 7;
        }
        out_ec.clear();
    }
    std::filesystem::create_directories(out_dir, out_ec);
    if (out_ec) {
        std::cerr << "Rank " << mpi_rank << " failed creating output directory: " << out_ec.message() << std::endl;
        return 7;
    }

    std::ofstream histfile;
    if (mpi_rank == 0 && cfg.hist_stride > 0) {
        histfile.open(binary_dir / "hist.dat", std::ios::out);
    }

    std::vector<double> vtk_targets;
    if (cfg.vtk_stride > 0) {
        vtk_targets.push_back(0.0);
        if (cfg.t_end > 0.0) {
            vtk_targets.push_back(0.5 * cfg.t_end);
            vtk_targets.push_back(cfg.t_end);
        }
    }
    std::vector<bool> vtk_written(vtk_targets.size(), false);

    for (size_t i = 0; i < vtk_targets.size(); ++i) {
        if (vtk_written[i])
            continue;
        if (sys.GetChTime() + 1e-15 < vtk_targets[i])
            continue;
        WriteConcreteSnapshotVTK(sys, out_dir, static_cast<int>(i + 1), 0, sys.GetChTime());
        if (cfg.vtk_stride > 0) {
            WriteCrackSnapshotVTK(sys, out_dir, static_cast<int>(i + 1), 0, sys.GetChTime());
        }
        vtk_written[i] = true;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const double t_wall_start = MPI_Wtime();

    int step = 0;
    int performed_steps = 0;
    double local_top_absdisp_peak = 0.0;
    double u = 0.0;
    double F = 0.0;
    double Wext = 0.0;
    int loop_partition_updates = 0;
    auto refresh_boundary_caches_after_partition = [&]() {
        local_top_all.clear();
        local_bottom_all.clear();
        local_top_bc_nodes.clear();
        local_bottom_bc_nodes.clear();
        local_top_nodes.clear();

        std::unordered_set<ChNodeFEAxyzrot*> refreshed_top_bc;
        std::unordered_set<ChNodeFEAxyzrot*> refreshed_bottom_bc;
        std::unordered_set<int> refreshed_owned_top_tags;
        for (auto& mesh : sys.GetMeshes()) {
            for (auto& node_base : mesh->GetNodes()) {
                auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
                if (!node || (!no_master_migration && !node->Variables().IsActive()))
                    continue;
                const int tag = node->GetTag();
                if (tag >= kTagTopNodeMin && tag <= kTagTopNodeMax) {
                    local_top_all.push_back(node);
                    if (refreshed_top_bc.insert(node.get()).second)
                        local_top_bc_nodes.push_back(node.get());
                    bool keep_owned = mpi_rank == owner_from_position(node->GetX0().GetPos());
                    if (!cfg.use_zslice_surface_owner)
                        keep_owned = get_node_weight(node.get()) > 0.0;
                    if (keep_owned && refreshed_owned_top_tags.insert(tag).second)
                        local_top_nodes.push_back(node.get());
                }
                if (tag >= kTagBottomNodeMin && tag <= kTagBottomNodeMax) {
                    local_bottom_all.push_back(node);
                    node->SetFixed((cfg.boundary_mode == BoundaryMode::FIX) && !use_motor_drive);
                    if (refreshed_bottom_bc.insert(node.get()).second)
                        local_bottom_bc_nodes.push_back(node.get());
                }
            }
        }

        const int local_owned_top = static_cast<int>(local_top_nodes.size());
        int global_owned_top = 0;
        MPI_Allreduce(&local_owned_top, &global_owned_top, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        return global_owned_top;
    };
    while (sys.GetChTime() < cfg.t_end) {
        if (cfg.do_partition_update_each_step) {
            domain_manager.DoDomainPartitionUpdate(mpi_rank, true);
            sys.Setup();
            ++loop_partition_updates;
            const int refreshed_top_count = refresh_boundary_caches_after_partition();
            if (refreshed_top_count != 134) {
                if (mpi_rank == 0) {
                    std::cerr << "Per-step partition update changed the owned top-node count: expected 134, got "
                              << refreshed_top_count << std::endl;
                }
                return 13;
            }
        }
        const double time_step_end = std::min(sys.GetChTime() + cfg.dt, cfg.t_end);
        const double poly_a_eff = cfg.poly_a * cfg.load_scale;
        const double disp_rate_eff = cfg.disp_rate * cfg.load_scale;
        const double t_switch = cfg.switch_time;
        const double u_switch = poly_a_eff * t_switch * t_switch;
        const double disp_cmd = (time_step_end <= t_switch) ? (poly_a_eff * time_step_end * time_step_end)
                                                             : (u_switch + disp_rate_eff * (time_step_end - t_switch));
        const double disp_rate_cmd = (time_step_end <= t_switch) ? (2.0 * poly_a_eff * time_step_end) : disp_rate_eff;
        if (!use_motor_drive) {
            for (auto* node : local_top_bc_nodes) {
                if (!node)
                    continue;
                const ChVector3d x0 = node->GetX0().GetPos();
                const ChVector3d p = node->GetPos();
                const ChVector3d v = node->GetPosDt();
                if (cfg.boundary_mode == BoundaryMode::FIX) {
                    // Constrained-top configuration (x,y fixed, imposed z).
                    node->SetPos(ChVector3d(x0.x(), x0.y(), x0.z() + disp_cmd));
                    node->SetPosDt(ChVector3d(0, 0, disp_rate_cmd));
                } else {
                    // Free-slip-top configuration (x,y free, imposed z).
                    node->SetPos(ChVector3d(p.x(), p.y(), x0.z() + disp_cmd));
                    node->SetPosDt(ChVector3d(v.x(), v.y(), disp_rate_cmd));
                }
                node->SetForce(ChVector3d(0, 0, 0));
                if (cfg.fix_top_rotations) {
                    node->SetRot(node->GetX0().GetRot());
                    node->SetAngVelLocal(ChVector3d(0, 0, 0));
                    node->SetTorque(ChVector3d(0, 0, 0));
                }
            }
            for (auto* node : local_bottom_bc_nodes) {
                if (!node)
                    continue;
                const ChVector3d x0 = node->GetX0().GetPos();
                const ChVector3d p = node->GetPos();
                const ChVector3d v = node->GetPosDt();
                if (cfg.boundary_mode == BoundaryMode::FIX) {
                    // Match the fixed lower reference plate used by the serial motor/tie setup.
                    node->SetPos(x0);
                    node->SetPosDt(ChVector3d(0, 0, 0));
                    node->SetRot(node->GetX0().GetRot());
                    node->SetAngVelLocal(ChVector3d(0, 0, 0));
                    node->SetTorque(ChVector3d(0, 0, 0));
                } else {
                    // Free-slip-bottom configuration (x,y free, fixed z).
                    node->SetPos(ChVector3d(p.x(), p.y(), x0.z()));
                    node->SetPosDt(ChVector3d(v.x(), v.y(), 0));
                }
                node->SetForce(ChVector3d(0, 0, 0));
            }
        }
        const double step_dt = std::max(0.0, time_step_end - sys.GetChTime());
        if (step_dt <= 0.0)
            break;
        sys.DoStepDynamics(step_dt);
        if (auto mdesc = std::dynamic_pointer_cast<ChSystemDescriptorMultidomain>(sys.GetSystemDescriptor())) {
            mdesc->SyncSharedStates(false);
        }
        ++step;
        ++performed_steps;

        for (size_t i = 0; i < vtk_targets.size(); ++i) {
            if (vtk_written[i])
                continue;
            if (sys.GetChTime() + 0.5 * cfg.dt < vtk_targets[i])
                continue;
            WriteConcreteSnapshotVTK(sys, out_dir, static_cast<int>(i + 1), step, sys.GetChTime());
            if (cfg.vtk_stride > 0) {
                WriteCrackSnapshotVTK(sys, out_dir, static_cast<int>(i + 1), step, sys.GetChTime());
            }
            vtk_written[i] = true;
        }

        double local_disp_weighted_sum = 0.0;
        double local_disp_weight_sum = 0.0;
        for (auto* node : local_top_nodes) {
            if (!node)
                continue;
            const double dz = node->GetPos().z() - node->GetX0().GetPos().z();
            if (!std::isfinite(dz))
                continue;
            local_disp_weighted_sum += dz;
            local_disp_weight_sum += 1.0;
            local_top_absdisp_peak = std::max(local_top_absdisp_peak, std::abs(dz));
        }

        double local_force = 0.0;
        ChVector3d local_support_from_internal(0, 0, 0);
        if (use_motor_drive) {
            if (local_motor && mpi_rank == motor_owner_rank) {
                local_force = local_motor->GetMotorForce();
            }
            if (local_constr_bot && mpi_rank == constr_owner_rank) {
                local_support_from_internal = local_constr_bot->GetReaction2().force;
            }
        } else {
            for (auto& mesh : sys.GetMeshes()) {
                for (int ie = 0; ie < mesh->GetNumElements(); ++ie) {
                    auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(ie));
                    if (!elem)
                        continue;

                    ChVectorDynamic<> Fi(24);
                    Fi.setZero();
                    elem->ComputeInternalForces(Fi);

                    for (int j = 0; j < 4; ++j) {
                        auto node = elem->GetTetrahedronNode(j);
                        if (!node)
                            continue;
                        const int ntag = node->GetTag();
                        const int ix = (j + 1) * 6 - 6;
                        const int iy = (j + 1) * 6 - 5;
                        const int iz = (j + 1) * 6 - 4;

                        // For reaction/support recovery, include all local element contributions
                        // connected to tagged top/bottom surface nodes. Elements are unique to ranks,
                        // so this avoids undercounting shared-node forces while remaining duplicate-safe.
                        const bool is_top = (ntag >= kTagTopNodeMin && ntag <= kTagTopNodeMax);
                        const bool is_bottom = (ntag >= kTagBottomNodeMin && ntag <= kTagBottomNodeMax);

                        if (is_top) {
                            local_force += (-Fi(iz));
                        }
                        if (is_bottom) {
                            local_support_from_internal += -ChVector3d(Fi(ix), Fi(iy), Fi(iz));
                        }
                    }
                }
            }
        }

        double global_disp_weighted_sum = 0.0;
        double global_disp_weight_sum = 0.0;
        double global_force = 0.0;
        MPI_Allreduce(&local_disp_weighted_sum, &global_disp_weighted_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_disp_weight_sum, &global_disp_weight_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_force, &global_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double global_disp = 0.0;
        if (global_disp_weight_sum > 0) {
            global_disp = global_disp_weighted_sum / global_disp_weight_sum;
        }
        if (!use_motor_drive) {
            global_disp = disp_cmd;
        }
        if (use_motor_drive) {
            double local_motor_pos = 0.0;
            if (local_motor && mpi_rank == motor_owner_rank) {
                local_motor_pos = local_motor->GetMotorPos();
            }
            double global_motor_pos = 0.0;
            MPI_Allreduce(&local_motor_pos, &global_motor_pos, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            global_disp = global_motor_pos;
        }
        if (mpi_rank == 0) {
            const double du = global_disp - u;
            Wext += std::abs(du * (global_force + F) * 0.5);
            u = global_disp;
            F = global_force;
        }

        const double disp_frac = std::abs(global_disp) / specimen_height;

        const bool do_hist_write = (cfg.hist_stride > 0) && ((step % cfg.hist_stride) == 0);
        double global_wint = 0.0;
        double global_ek = 0.0;
        ChVector3d global_support(0, 0, 0);
        ChVector3d global_motor_reaction(0, 0, 0);
        if (do_hist_write) {
            double local_wint = 0.0;
            double local_ek = 0.0;
            ChVector3d local_support(0, 0, 0);
            local_support = local_support_from_internal;
            ComputeLocalEnergies(sys, local_wint, local_ek);
            MPI_Allreduce(&local_wint, &global_wint, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_ek, &global_ek, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double local_support_buf[3] = {local_support.x(), local_support.y(), local_support.z()};
            double global_support_buf[3] = {0, 0, 0};
            MPI_Allreduce(local_support_buf, global_support_buf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            global_support = ChVector3d(global_support_buf[0], global_support_buf[1], global_support_buf[2]);
            global_motor_reaction = ChVector3d(0, 0, global_force);
        }

        if (mpi_rank == 0 && do_hist_write) {
            std::cout << " t=\t" << sys.GetChTime() << "\ttop_plate disp_z=\t" << global_disp << "\tforce=\t"
                      << global_force << "\tpressure=\t0"
                      << "\tSupport=\t" << global_support.z() << "\tinternal work\t" << global_wint
                      << "\texternal_work\t" << Wext << "\tkinetic energy\t" << global_ek
                      << "\tmotor_reaction=\t" << global_motor_reaction.z() << "\t\n";
            histfile << " t=\t" << sys.GetChTime() << "\ttop_plate disp_z=\t" << global_disp << "\tforce=\t"
                     << global_force << "\tpressure=\t0"
                     << "\tSupport=\t" << global_support << "\tinternal work\t" << global_wint
                     << "\texternal_work\t" << Wext << "\tkinetic energy\t" << global_ek
                     << "\tmotor_reaction=\t" << global_motor_reaction << "\t\n";
            histfile.flush();
        }

        if (cfg.stop_disp_frac > 0.0 && disp_frac >= cfg.stop_disp_frac) {
            break;
        }
    }

    auto report_node_state_spread = [&](const char* label, bool top_only) {
        std::vector<int> local_tags;
        std::vector<double> local_z;
        for (const auto& mesh : sys.GetMeshes()) {
            for (const auto& node_base : mesh->GetNodes()) {
                auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node_base);
                if (!node)
                    continue;
                const int tag = node->GetTag();
                if (top_only && !(tag >= kTagTopNodeMin && tag <= kTagTopNodeMax))
                    continue;
                local_tags.push_back(tag);
                local_z.push_back(node->GetPos().z());
            }
        }

        const int local_count = static_cast<int>(local_tags.size());
        std::vector<int> recv_counts(mpi_total_ranks, 0);
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> displs(mpi_total_ranks, 0);
        int total_entries = 0;
        for (int r = 0; r < mpi_total_ranks; ++r) {
            displs[r] = total_entries;
            total_entries += recv_counts[r];
        }

        std::vector<int> all_tags(total_entries, -1);
        std::vector<double> all_z(total_entries, 0.0);
        MPI_Allgatherv(local_tags.data(), local_count, MPI_INT, all_tags.data(), recv_counts.data(), displs.data(),
                       MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(local_z.data(), local_count, MPI_DOUBLE, all_z.data(), recv_counts.data(), displs.data(),
                       MPI_DOUBLE, MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            struct SpreadStats {
                int count = 0;
                double min_z = std::numeric_limits<double>::max();
                double max_z = std::numeric_limits<double>::lowest();
            };

            std::unordered_map<int, SpreadStats> stats;
            stats.reserve(static_cast<size_t>(total_entries) * 2 + 1);
            for (int i = 0; i < total_entries; ++i) {
                const int tag = all_tags[i];
                if (tag < 0)
                    continue;
                auto& s = stats[tag];
                s.count += 1;
                s.min_z = std::min(s.min_z, all_z[i]);
                s.max_z = std::max(s.max_z, all_z[i]);
            }

            std::vector<double> spreads;
            spreads.reserve(stats.size());
            for (const auto& kv : stats) {
                if (kv.second.count > 1) {
                    spreads.push_back(std::abs(kv.second.max_z - kv.second.min_z));
                }
            }

            double max_spread = 0.0;
            double mean_spread = 0.0;
            double p95_spread = 0.0;
            if (!spreads.empty()) {
                max_spread = *std::max_element(spreads.begin(), spreads.end());
                for (double v : spreads)
                    mean_spread += v;
                mean_spread /= static_cast<double>(spreads.size());
                std::sort(spreads.begin(), spreads.end());
                const size_t p95_idx = static_cast<size_t>(0.95 * static_cast<double>(spreads.size() - 1));
                p95_spread = spreads[p95_idx];
            }

            std::cout << "[mdom-state-spread] " << label
                      << " duplicated_tags=" << spreads.size()
                      << " max_abs_dz=" << max_spread
                      << " mean_abs_dz=" << mean_spread
                      << " p95_abs_dz=" << p95_spread << std::endl;
        }
    };

    report_node_state_spread("all_nodes", false);
    report_node_state_spread("top_nodes", true);

    MPI_Barrier(MPI_COMM_WORLD);
    const double local_wall_seconds = MPI_Wtime() - t_wall_start;

    double wall_seconds = 0.0;
    MPI_Reduce(&local_wall_seconds, &wall_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    std::vector<double> gathered_top_absdisp_peak;
    if (mpi_rank == 0) {
        gathered_top_absdisp_peak.resize(mpi_total_ranks, 0.0);
    }
    MPI_Gather(&local_top_absdisp_peak, 1, MPI_DOUBLE,
               mpi_rank == 0 ? gathered_top_absdisp_peak.data() : nullptr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        const int total_ranks = worker_domains + 1;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "RUN_SUMMARY worker_domains=" << worker_domains << " total_ranks=" << total_ranks
                  << " decomposition=" << (cfg.use_grid_decomposition ? "grid" : "slices") << " grid="
                  << grid_divisions[0] << "x" << grid_divisions[1] << "x" << grid_divisions[2]
                  << " worker_graph=" << (cfg.use_complete_worker_graph ? "complete" : "neighbor")
                  << " loop_partition_updates=" << loop_partition_updates
                  << " dt=" << cfg.dt << " t_end=" << cfg.t_end << " steps=" << performed_steps << std::endl;
        std::cout << "Top-node peak |dz| per rank:";
        for (int r = 0; r < mpi_total_ranks; ++r)
            std::cout << " r" << r << "=" << gathered_top_absdisp_peak[r];
        std::cout << std::endl;
        std::cout << "WALL_SECONDS=" << wall_seconds << std::endl;

        std::ofstream summary(out_dir / "run_summary.txt", std::ios::out);
        summary << "worker_domains=" << worker_domains << "\n";
        summary << "total_ranks=" << total_ranks << "\n";
        summary << "decomposition=" << (cfg.use_grid_decomposition ? "grid" : "slices") << "\n";
        summary << "grid=" << grid_divisions[0] << "x" << grid_divisions[1] << "x" << grid_divisions[2] << "\n";
        summary << "worker_graph=" << (cfg.use_complete_worker_graph ? "complete" : "neighbor") << "\n";
        summary << "loop_partition_updates=" << loop_partition_updates << "\n";
        summary << "dt=" << cfg.dt << "\n";
        summary << "t_end=" << cfg.t_end << "\n";
        summary << "steps=" << performed_steps << "\n";
        summary << "loading_mode=displacement_control\n";
        summary << "target_top_displacement=" << cfg.target_top_displacement << "\n";
        summary << "top_surface_area=" << top_surface_area << "\n";
        summary << "specimen_height=" << specimen_height << "\n";
        summary << "wall_seconds=" << wall_seconds << "\n";

        if (cfg.vtk_stride > 0) {
            WriteRankedSnapshotPVD(binary_dir, out_dir_prefix, "concrete_snapshot", vtk_targets,
                                   binary_dir / ("concrete_snapshots_w" + std::to_string(worker_domains) + ".pvd"),
                                   mpi_total_ranks);
            WriteRankedSnapshotPVD(binary_dir, out_dir_prefix, "crack_snapshot", vtk_targets,
                                   binary_dir / ("crack_snapshots_w" + std::to_string(worker_domains) + ".pvd"),
                                   mpi_total_ranks);
        }
    }

    if (histfile.is_open()) {
        histfile.close();
    }

    return 0;
}
