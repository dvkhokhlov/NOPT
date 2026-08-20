#pragma once
//
// dmrg_log — one block2 sweep log per solve, written under dmrg/.
//
// block2 prints to std::cout/std::cerr with no sink indirection, so solves are separated by
// swapping both rdbufs for the duration of one call. Restoring on every path is mandatory:
// block2 throws through its solves.

#include <fstream>
#include <string>

// block2 verbosity inside the log: 1 = per-sweep summary, 2 = + per-site lines and sweep timings.
constexpr int DMRG_LOG_IPRINT = 2;

// Which solve is running. Each tag names its own file series and carries its own counter.
enum class dmrg_log_tag { primary, canonical, scf, cdas };

// $PAR NAME, the log-file stem. Resets the counters when it changes.
void dmrg_log_set_job(const char* name);

// Announce the next solve. Sticky until the next call.
void dmrg_log_set_tag(dmrg_log_tag tag);

// RAII redirect of cout+cerr into dmrg/<job>_dmrg[_tag]<n>.log. fresh opens the next file for
// the current tag; !fresh appends to the last file opened, so the RDM sweeps land with the solve
// whose wavefunction they read. A failed open leaves both streams untouched.
class dmrg_log_guard {
public:
    explicit dmrg_log_guard(bool fresh);
    ~dmrg_log_guard();
    dmrg_log_guard(const dmrg_log_guard&) = delete;
    dmrg_log_guard& operator=(const dmrg_log_guard&) = delete;

private:
    std::ofstream   file_;
    std::streambuf* old_out_ = nullptr;
    std::streambuf* old_err_ = nullptr;
};
