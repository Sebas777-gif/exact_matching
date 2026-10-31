#pragma once

// mpfree + AIDA preprocessing for one input file.
//
//   .bif                         → chain → mpfree → AIDA
//   .scc / .scc_k / .minpres_k   → AIDA (already a presentation)
//   .sccsum                      → already decomposed; read as-is
//
// Each file is one module. The matching distance uses the union of summand
// grades for A0 and the concatenated slice barcodes for F.

#include "bif_io.hpp"
#include "scc_io.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace scc_hera {

inline std::string shell_quote(const std::string& s)
{
    std::string out = "'";
    for (char c : s) {
        if (c == '\'')
            out += "'\\''";
        else
            out += c;
    }
    out += "'";
    return out;
}

inline bool file_is_executable(const std::string& path)
{
    return !path.empty() && access(path.c_str(), X_OK) == 0;
}

inline std::string find_binary(const char* env_name,
                               const std::vector<std::string>& candidates,
                               const char* path_name)
{
    if (const char* env = std::getenv(env_name)) {
        if (*env)
            return env;
    }
    if (const char* home = std::getenv("HOME")) {
        for (const std::string& rel : candidates) {
            const std::string p = std::string(home) + rel;
            if (file_is_executable(p))
                return p;
        }
    }
    return path_name;
}

inline std::string find_mpfree()
{
    return find_binary("MPFREE",
        {
            "/hera_aida_matching/extern/mpfree/build/mpfree",
            "/mpfree/build/mpfree",
        },
        "mpfree");
}

inline std::string find_aida()
{
    if (const char* env = std::getenv("AIDA_WITHOUT_STATS")) {
        if (*env)
            return env;
    }
    return find_binary("AIDA",
        {
            "/AIDA/build/aida",
        },
        "aida");
}

inline void run_cmd(const std::string& cmd, const std::string& what)
{
    std::cerr << what << ": " << cmd << std::endl;
    const int status = std::system(cmd.c_str());
    if (status != 0)
        throw std::runtime_error(what + " failed (status " + std::to_string(status)
            + "). Set $MPFREE / $AIDA, or pass --no-aida.");
}

inline std::string temp_stem(const std::string& input)
{
    static std::string dir;
    if (dir.empty()) {
        char buf[] = "/tmp/exact_matching_XXXXXX";
        if (!mkdtemp(buf))
            throw std::runtime_error("mkdtemp failed");
        dir = buf;
    }
    auto slash = input.find_last_of('/');
    std::string base = slash == std::string::npos ? input : input.substr(slash + 1);
    return dir + "/" + base;
}

inline matching::Module presentations_to_module(std::vector<matching::Presentation> summands)
{
    matching::Module M;
    M.summands = std::move(summands);
    return M;
}

inline matching::Module run_aida_on_presentation(const std::string& scc_path)
{
    const std::string aida = find_aida();
    const std::string out = temp_stem(scc_path) + ".sccsum";
    const std::string cmd = shell_quote(aida) + " -l -o " + shell_quote(out)
        + " " + shell_quote(scc_path);
    run_cmd(cmd, "aida");
    std::ifstream check(out);
    if (!check)
        throw std::runtime_error("aida wrote no summands for " + scc_path);
    return presentations_to_module(read_scc_presentations(out));
}

inline std::string mpfree_presentation_from_bif(const std::string& bif_path, int homology_dim)
{
    const std::string stem = temp_stem(bif_path);
    const std::string chain = stem + ".chain.scc";
    const std::string presentation = stem + ".scc";
    write_chain_complex_scc(chain, read_bif(bif_path), homology_dim);

    const std::string mpfree = find_mpfree();
    const std::string cmd = shell_quote(mpfree) + " " + shell_quote(chain)
        + " " + shell_quote(presentation);
    run_cmd(cmd, "mpfree");
    std::ifstream check(presentation);
    if (!check)
        throw std::runtime_error("mpfree wrote no presentation for " + bif_path);
    return presentation;
}

// Load one module from a file. With preprocess=true (the default), .bif
// goes through mpfree then AIDA, .scc through AIDA, and .sccsum is read
// as a list of summands. With preprocess=false, .bif uses the in-process
// d-simplex presentation and .scc/.sccsum blocks are summands as written.
inline matching::Module load_module(const std::string& path, int homology_dim,
                                    bool preprocess)
{
    if (suffix_of(path) == ".sccsum")
        return presentations_to_module(read_scc_presentations(path));

    if (!preprocess) {
        if (is_bif_path(path))
            return matching::as_module(load_bif_presentation(path, homology_dim));
        return presentations_to_module(read_scc_presentations(path));
    }

    if (is_bif_path(path)) {
        const std::string scc = mpfree_presentation_from_bif(path, homology_dim);
        return run_aida_on_presentation(scc);
    }
    if (is_scc_path(path))
        return run_aida_on_presentation(path);
    throw std::runtime_error("unrecognized input " + path
        + " (expected .bif, .scc, .scc_k, .minpres_k, or .sccsum)");
}

}  // namespace scc_hera
