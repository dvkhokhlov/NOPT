#pragma once
//
// block2_gpu_guard — the block2 GPU backend engaged for one two-site DMRG solve. Include
// ONLY from a block2 backend TU, after block2_core.hpp / block2_dmrg.hpp; without
// NOPT_BLOCK2_GPU every member compiles away, so callers carry no #ifdef of their own.

#ifdef NOPT_BLOCK2_GPU
#include "gpu/gpu_consumer.hpp"   // engage/disengage/engaged/stats trio exported by libblock2
#endif

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include "common_vars.h"    // out_stream
#include "inp_par_read.h"   // DMRG_WARM_ON

// The backend declines every step whose center arity is not two, so the one-site tail must run
// disengaged: finish() ends the engagement there, the destructor covers the other exits. The
// moving environment must outlive the guard — the restore goes back through it.
class block2_gpu_guard {
public:
    block2_gpu_guard(
        int gpu_opt,
        const std::shared_ptr<block2::MovingEnvironment<block2::SU2, double, double>> &me) {
#ifdef NOPT_BLOCK2_GPU
        if (gpu_opt != DMRG_WARM_ON)
            return;
        // No overwrite: an explicit BLOCK2_GPU_MODE (capture, crt15, strict) wins over the knob.
        setenv("BLOCK2_GPU_MODE", "on", 0);
        bool ok = false;
        try {
            ok = block2::gpu_consumer_engage(me);
        } catch (const std::exception &ex) {
            fprintf(out_stream, "ERROR: $DMRG gpu=on: the GPU backend refused to engage: %s\n",
                    ex.what());
            exit(EXIT_FAILURE);
        }
        if (!ok) {
            fprintf(out_stream, "ERROR: $DMRG gpu=on but the backend stayed disengaged"
                                " -- BLOCK2_GPU_MODE in the environment?\n");
            exit(EXIT_FAILURE);
        }
        engaged_ = true;
#else
        (void)gpu_opt;
        (void)me;
#endif
    }

    ~block2_gpu_guard() {
        try { finish(); } catch (...) {} // never throw from a destructor
    }

    // Disengage and report; idempotent, and a no-op when gpu=off never engaged.
    void finish() {
#ifdef NOPT_BLOCK2_GPU
        if (!engaged_)
            return;
        engaged_ = false;
        block2::gpu_consumer_disengage();
        const std::string stats = block2::gpu_consumer_stats_line();
        std::cout << stats << std::endl;             // per-solve dmrg/*.log
        fprintf(out_stream, "%s\n", stats.c_str());  // main output
#endif
    }

    block2_gpu_guard(const block2_gpu_guard &) = delete;
    block2_gpu_guard &operator=(const block2_gpu_guard &) = delete;

private:
#ifdef NOPT_BLOCK2_GPU
    bool engaged_ = false;
#endif
};
