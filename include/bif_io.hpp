#pragma once

// Bifiltration reader adapted from hera_aida_matching/bif_to_scc.hpp.
// Expected .bif (hera/matching/README.md), comments starting with # ignored:
//
//   bifiltration_phat_like
//   <N>
//   dim x y i j k          # N lines, 0-based boundary indices
//
// Shape-dataset files may omit the keyword and start with N.

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <matching_distance.h>

namespace scc_hera {

struct BifSimplex {
    int dim = 0;
    double x = 0;
    double y = 0;
    std::vector<int> boundary;
};

inline std::string suffix_of(const std::string& path)
{
    auto dot = path.find_last_of('.');
    if (dot == std::string::npos)
        return {};
    std::string s = path.substr(dot);
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

inline bool is_bif_path(const std::string& path)
{
    return suffix_of(path) == ".bif";
}

// Kerber min-pres datasets name files foo.scc_2 / foo.scc.minpres_0.
// The trailing integer is the chain-complex position (reflected homology:
// .scc_2 is already an H_0 presentation). It is not --dim.
inline bool scc_suffix_with_index(const std::string& s, const char* prefix)
{
    const std::string p = prefix;
    if (s.size() < p.size() + 2)
        return false;
    if (s.compare(0, p.size(), p) != 0)
        return false;
    if (s[p.size()] != '_')
        return false;
    for (std::size_t i = p.size() + 1; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i])))
            return false;
    }
    return true;
}

inline bool is_scc_path(const std::string& path)
{
    auto s = suffix_of(path);
    return s == ".scc" || s == ".sccsum"
        || scc_suffix_with_index(s, ".scc")
        || scc_suffix_with_index(s, ".minpres");
}

inline std::vector<BifSimplex> read_bif(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open " + path);

    std::string line;
    auto next = [&]() -> bool {
        while (std::getline(in, line)) {
            auto hash = line.find('#');
            if (hash != std::string::npos)
                line = line.substr(0, hash);
            auto not_space = [](unsigned char c) { return !std::isspace(c); };
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
            if (!line.empty())
                return true;
        }
        return false;
    };

    if (!next())
        throw std::runtime_error(path + ": empty bifiltration");
    int n = 0;
    if (line == "bifiltration_phat_like") {
        if (!next())
            throw std::runtime_error(path + ": missing simplex count");
        std::istringstream iss(line);
        iss >> n;
        if (!iss || n < 0)
            throw std::runtime_error(path + ": bad simplex count");
    } else {
        std::istringstream iss(line);
        iss >> n;
        if (!iss || n < 0)
            throw std::runtime_error(path + ": expected bifiltration_phat_like or a simplex count, got '"
                + line + "'");
    }

    std::vector<BifSimplex> simplices;
    simplices.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!next())
            throw std::runtime_error(path + ": truncated at simplex " + std::to_string(i));
        std::istringstream iss(line);
        BifSimplex s;
        iss >> s.dim >> s.x >> s.y;
        if (!iss)
            throw std::runtime_error(path + ": bad simplex line: " + line);
        int idx = 0;
        while (iss >> idx) {
            if (idx < 0 || idx >= n)
                throw std::runtime_error(path + ": boundary index out of range");
            s.boundary.push_back(idx);
        }
        simplices.push_back(std::move(s));
    }
    return simplices;
}

// Presentation of H_d by d-simplices (generators) and (d+1)-simplices
// (relations). For d = 0 this is the usual presentation of H_0.
inline matching::Presentation presentation_from_bif(const std::vector<BifSimplex>& simplices,
                                                    int homology_dim)
{
    if (homology_dim < 0)
        throw std::runtime_error("homology dimension must be non-negative");

    const int n = static_cast<int>(simplices.size());
    std::vector<int> rank(n, -1);
    int n_gen = 0;
    for (int i = 0; i < n; ++i) {
        if (simplices[i].dim == homology_dim)
            rank[i] = n_gen++;
    }

    matching::Presentation P;
    P.generators.resize(n_gen);
    int gi = 0;
    for (const BifSimplex& s : simplices) {
        if (s.dim != homology_dim)
            continue;
        P.generators[gi].x = matching::FT(s.x);
        P.generators[gi].y = matching::FT(s.y);
        ++gi;
    }

    for (const BifSimplex& s : simplices) {
        if (s.dim != homology_dim + 1)
            continue;
        matching::Grade2 g;
        g.x = matching::FT(s.x);
        g.y = matching::FT(s.y);
        P.relations.push_back(g);
        std::vector<int> col;
        for (int b : s.boundary) {
            if (b < 0 || b >= n)
                continue;
            if (simplices[b].dim != homology_dim)
                continue;
            col.push_back(rank[b]);
        }
        matching::reduce_mod2(col);
        P.columns.push_back(std::move(col));
    }
    return P;
}

inline matching::Presentation load_bif_presentation(const std::string& path, int homology_dim)
{
    return presentation_from_bif(read_bif(path), homology_dim);
}

// Write C_{d+1} --∂_{d+1}--> C_d --∂_d--> C_{d-1} as scc2020 so that
// mpfree's default --dim=1 presents H_d = ker(∂_d)/im(∂_{d+1}).
//
// The size line is n_{d+1} n_d n_{d-1}. For d=0 there are no (-1)-simplices,
// so the third count is 0 and the file is already a presentation of the
// unminimized H_0. The previous writer used n_d as the third count and, for
// d=0, wrote the vertices twice; at d=2 that put edge indices into an
// n_triangles-row matrix and mpfree segfaulted.
inline void write_chain_complex_scc(const std::string& path,
                                    const std::vector<BifSimplex>& simplices,
                                    int homology_dim)
{
    int max_dim = 0;
    for (const auto& s : simplices)
        max_dim = std::max(max_dim, s.dim);
    if (homology_dim < 0 || homology_dim > max_dim)
        throw std::runtime_error("homology dimension out of range for this bifiltration");

    const int high = homology_dim + 1;
    const int mid = homology_dim;
    const int low = homology_dim - 1;
    const int n = static_cast<int>(simplices.size());
    std::vector<int> rank(n, -1);
    std::vector<int> by_dim(max_dim + 2, 0);
    for (int i = 0; i < n; ++i) {
        rank[i] = by_dim[simplices[i].dim];
        ++by_dim[simplices[i].dim];
    }
    const int n_high = high <= max_dim ? by_dim[high] : 0;
    const int n_mid = by_dim[mid];
    const int n_low = low >= 0 ? by_dim[low] : 0;

    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("cannot write " + path);
    out.setf(std::ios::fixed);
    out.precision(8);
    out << "scc2020\n2\n" << n_high << " " << n_mid << " " << n_low << "\n";

    auto write_dim = [&](int dim, int face_dim) {
        for (int i = 0; i < n; ++i) {
            if (simplices[i].dim != dim)
                continue;
            out << simplices[i].x << " " << simplices[i].y << " ;";
            if (face_dim >= 0) {
                for (int b : simplices[i].boundary) {
                    if (b < 0 || b >= n)
                        continue;
                    if (simplices[b].dim != face_dim)
                        continue;
                    out << " " << rank[b];
                }
            }
            out << "\n";
        }
    };
    if (n_high > 0)
        write_dim(high, mid);
    write_dim(mid, low);
    if (n_low > 0) {
        for (const auto& s : simplices)
            if (s.dim == low)
                out << s.x << " " << s.y << " ;\n";
    }
}

}  // namespace scc_hera
