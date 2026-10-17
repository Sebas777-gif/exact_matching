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

inline bool is_scc_path(const std::string& path)
{
    auto s = suffix_of(path);
    return s == ".scc" || s == ".sccsum";
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
        std::sort(col.begin(), col.end());
        col.erase(std::unique(col.begin(), col.end()), col.end());
        P.columns.push_back(std::move(col));
    }
    return P;
}

inline matching::Presentation load_bif_presentation(const std::string& path, int homology_dim)
{
    return presentation_from_bif(read_bif(path), homology_dim);
}

}  // namespace scc_hera
