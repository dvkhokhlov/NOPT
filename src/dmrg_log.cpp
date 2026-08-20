#include "dmrg_log.h"

#include <filesystem>
#include <iostream>
#include <system_error>

namespace {

const char* suffix_of(dmrg_log_tag t) {
    switch (t) {
    case dmrg_log_tag::canonical: return "_can";
    case dmrg_log_tag::scf:       return "_scf";
    case dmrg_log_tag::cdas:      return "_cdas";
    default:                      return "";
    }
}

std::string  g_job = "nopt";
dmrg_log_tag g_tag = dmrg_log_tag::primary;
int          g_count[4] = {0, 0, 0, 0};
std::string  g_last;   // file the last fresh guard opened; !fresh appends to it

} // namespace

void dmrg_log_set_job(const char* name) {
    const std::string n = (name != nullptr && *name != '\0') ? name : "nopt";
    if (n == g_job)
        return;
    g_job = n;
    for (int& c : g_count)
        c = 0;
    g_last.clear();
}

void dmrg_log_set_tag(dmrg_log_tag tag) { g_tag = tag; }

dmrg_log_guard::dmrg_log_guard(bool fresh) {
    std::string path;
    if (fresh) {
        std::error_code ec;
        std::filesystem::create_directory("dmrg", ec);
        const int i = static_cast<int>(g_tag);
        path = "dmrg/" + g_job + "_dmrg" + suffix_of(g_tag) + std::to_string(g_count[i]++) + ".log";
        file_.open(path, std::ios::out | std::ios::trunc);
    } else if (!g_last.empty()) {
        path = g_last;
        file_.open(path, std::ios::out | std::ios::app);
    }
    if (!file_.is_open())
        return;
    if (fresh)
        g_last = path;
    old_out_ = std::cout.rdbuf(file_.rdbuf());
    old_err_ = std::cerr.rdbuf(file_.rdbuf());
}

dmrg_log_guard::~dmrg_log_guard() {
    if (old_out_ != nullptr)
        std::cout.rdbuf(old_out_);
    if (old_err_ != nullptr)
        std::cerr.rdbuf(old_err_);
}
