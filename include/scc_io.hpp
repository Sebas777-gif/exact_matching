#pragma once

// Read an scc2020 presentation, or several concatenated ones (as in an
// AIDA .sccsum file), into matching::Presentation objects.
//
// Expected block:
//
//   scc2020
//   2
//   <n_relations> <n_generators> 0
//   x y ; i j k          # n_relations lines
//   x y ;                # n_generators lines
//
// Empty lines, comments, and AIDA type tags
// (free/cyclic/interval/non-interval, sccsum, scc2020sum) are skipped
// between blocks.

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <matching_distance.h>

namespace scc_hera {

inline std::string trim_scc(std::string s)
{
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

inline bool is_aida_tag(const std::string& s)
{
    return s == "free" || s == "cyclic" || s == "interval" || s == "non-interval"
        || s == "sccsum" || s == "scc2020sum";
}

inline bool is_count_line(const std::string& s)
{
    if (s.empty())
        return false;
    for (char c : s)
        if (!std::isdigit(static_cast<unsigned char>(c)))
            return false;
    return true;
}

// `in` is positioned just after a line that was exactly "scc2020".
inline matching::Presentation read_scc_block_body(std::istream& in, int& line_no)
{
    std::string line;
    auto next = [&]() -> bool {
        while (std::getline(in, line)) {
            ++line_no;
            line = trim_scc(line);
            if (line.empty() || line[0] == '#')
                continue;
            return true;
        }
        return false;
    };

    if (!next())
        throw std::runtime_error("truncated scc2020 block (parameter count)");
    int n_parameters = 0;
    {
        std::istringstream iss(line);
        iss >> n_parameters;
        if (!iss || n_parameters != 2)
            throw std::runtime_error("line " + std::to_string(line_no)
                + ": matching distance needs a 2-parameter module, got '" + line + "'");
    }

    if (!next())
        throw std::runtime_error("truncated scc2020 block (shape line)");
    int n_rel = 0, n_gen = 0, n_higher = 0;
    {
        std::istringstream iss(line);
        iss >> n_rel >> n_gen >> n_higher;
        if (!iss)
            throw std::runtime_error("line " + std::to_string(line_no)
                + ": expected '<n_relations> <n_generators> 0'");
        if (n_higher != 0)
            throw std::runtime_error("line " + std::to_string(line_no)
                + ": this is a chain complex, not a presentation (third count is "
                + std::to_string(n_higher) + "). Run mpfree first.");
    }

    matching::Presentation P;
    P.relations.reserve(n_rel);
    P.columns.reserve(n_rel);
    for (int i = 0; i < n_rel; ++i) {
        if (!next())
            throw std::runtime_error("truncated relations at line " + std::to_string(line_no));
        auto semi = line.find(';');
        std::string grade = semi == std::string::npos ? line : line.substr(0, semi);
        std::string entries = semi == std::string::npos ? std::string() : line.substr(semi + 1);
        std::istringstream gs(grade), es(entries);
        double x = 0, y = 0;
        gs >> x >> y;
        if (!gs)
            throw std::runtime_error("line " + std::to_string(line_no) + ": bad relation grade");
        matching::Grade2 g;
        g.x = matching::FT(x);
        g.y = matching::FT(y);
        P.relations.push_back(g);
        std::vector<int> col;
        int idx = 0;
        while (es >> idx) {
            if (idx < 0 || idx >= n_gen)
                throw std::runtime_error("line " + std::to_string(line_no)
                    + ": generator index " + std::to_string(idx) + " out of range");
            col.push_back(idx);
        }
        matching::reduce_mod2(col);
        P.columns.push_back(std::move(col));
    }

    P.generators.reserve(n_gen);
    for (int i = 0; i < n_gen; ++i) {
        if (!next())
            throw std::runtime_error("truncated generators at line " + std::to_string(line_no));
        auto semi = line.find(';');
        std::string grade = semi == std::string::npos ? line : line.substr(0, semi);
        std::istringstream gs(grade);
        double x = 0, y = 0;
        gs >> x >> y;
        if (!gs)
            throw std::runtime_error("line " + std::to_string(line_no) + ": bad generator grade");
        matching::Grade2 g;
        g.x = matching::FT(x);
        g.y = matching::FT(y);
        P.generators.push_back(g);
    }
    return P;
}

inline std::vector<matching::Presentation> read_scc_presentations(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open " + path);
    std::vector<matching::Presentation> out;
    int line_no = 0;
    std::string line;
    while (std::getline(in, line)) {
        ++line_no;
        line = trim_scc(line);
        if (line.empty() || line[0] == '#' || is_aida_tag(line) || is_count_line(line))
            continue;
        if (line != "scc2020")
            throw std::runtime_error("line " + std::to_string(line_no)
                + ": expected scc2020, got '" + line + "'");
        out.push_back(read_scc_block_body(in, line_no));
    }
    if (out.empty())
        throw std::runtime_error(path + ": no scc2020 block found");
    return out;
}

}  // namespace scc_hera
