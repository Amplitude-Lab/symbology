# compiler
# On macOS the default g++ is Apple Clang, which cannot compile SparseRREF
# (it uses std::execution::par, missing from Apple's libc++). Prefer a real
# GCC from Homebrew when available; override anytime with `make CXX=...`.
ifeq ($(shell uname -s),Darwin)
HOMEBREW_GXX := $(firstword $(sort $(wildcard /opt/homebrew/bin/g++-[0-9]* /usr/local/bin/g++-[0-9]*)))
CXX := $(if $(HOMEBREW_GXX),$(HOMEBREW_GXX),g++)
else
CXX := g++
endif
# -mcpu=native is ARM-speak; x86_64 GCC wants -march=native to actually
# dispatch to the host's vector ISA (Linux HPC boxes are x86_64).
ifeq ($(shell uname -m),x86_64)
ARCH_FLAGS := -march=native -mtune=native
else
ARCH_FLAGS := -mcpu=native -mtune=native
endif
CXXFLAGS := -O3 -flto $(ARCH_FLAGS) -std=c++20 -I.

# mimalloc is a performance add-on (it overrides the GMP/FLINT allocators
# when USE_MIMALLOC=1); HPC boxes often lack it, so link it only when a copy
# exists in the usual library locations. The math is identical without it.
# Override explicitly with MIMALLOC=1 (force) or MIMALLOC=0 (forbid).
MIMALLOC_LIB := $(firstword $(wildcard /opt/homebrew/lib/libmimalloc.* /usr/local/lib/libmimalloc.* /usr/lib/*/libmimalloc.* /usr/lib/libmimalloc.*))
ifeq ($(MIMALLOC),0)
LDLIBS := -lflint -lgmp
else ifneq ($(MIMALLOC_LIB)$(MIMALLOC),)
LDLIBS := -lflint -lgmp -lmimalloc
else
LDLIBS := -lflint -lgmp
endif

ifeq ($(OS),Windows_NT)
EXE := bootstrap.exe
LDLIBS += -ltbb12
else
EXE := bootstrap
LDLIBS += -ltbb
endif

# macOS + Homebrew: neither Apple Clang nor Homebrew GCC searches the Homebrew
# prefix by default, so flint/gmp/tbb/mimalloc installed by brew are invisible
# without -I/-L flags (symptoms: "flint/nmod.h: No such file or directory" at
# compile, "ld: library 'flint' not found" at link). Add the paths ourselves
# so a plain `make` works in ANY shell — CPATH/LIBRARY_PATH exports in a login
# profile (added by some Homebrew mirror setup scripts) are not guaranteed in
# CI runners, IDE build tasks, or non-login shells. No-op when brew or its
# flint headers are absent; override anytime with `make CXXFLAGS=... LDLIBS=...`.
ifeq ($(shell uname -s),Darwin)
BREW_BIN := $(shell command -v brew 2>/dev/null)
ifeq ($(BREW_BIN),)
BREW_BIN := $(firstword $(wildcard /opt/homebrew/bin/brew /usr/local/bin/brew))
endif
HOMEBREW_PREFIX := $(shell $(BREW_BIN) --prefix 2>/dev/null)
ifneq ($(wildcard $(HOMEBREW_PREFIX)/include/flint),)
CXXFLAGS += -I$(HOMEBREW_PREFIX)/include
LDLIBS := -L$(HOMEBREW_PREFIX)/lib -Wl,-rpath,$(HOMEBREW_PREFIX)/lib $(LDLIBS)
endif
endif

# SparseRREF sanity: the bundled race fixes (patches/) are required.
# rref publishes rewritten rows via plain int flags upstream; without the
# release/acquire fix the build crashes intermittently (SIGSEGV, ~10% of
# runs at weight >= 3). Warn loudly instead of failing silently.
ifeq ($(wildcard SparseRREF/sparse_mat.h),)
$(warning SparseRREF/ not found — clone it and apply patches/ (see README); see https://github.com/munuxi/SparseRREF)
else
ifneq ($(shell grep -c "std::atomic<int>" SparseRREF/sparse_mat.h),2)
$(warning SparseRREF lacks the required race fixes — apply patches/sparserref-5bbee55-race-and-init-fix.patch (or the race-only patch on newer upstream); intermittent SIGSEGV otherwise. See README.)
endif
endif

# target
all: $(EXE)

# executable file
$(EXE): bootstrap.cpp bootstrap.hpp projection.hpp solve_symmetry.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp tensor_shuffle.h
	$(CXX) bootstrap.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# compute_rhs executable (standalone RHS computation module)
compute_rhs: compute_rhs.cpp compute_rhs.hpp tensor_shuffle.h bootstrap.hpp projection.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp
	$(CXX) compute_rhs.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# tensor_add executable (weighted sum of sparse tensors)
tensor_add: tensor_add.cpp tensor_shuffle.h
	$(CXX) tensor_add.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# tensor_ops executable (ternary contraction, matrix power, tensor join, squeeze)
tensor_ops: tensor_ops.cpp
	$(CXX) tensor_ops.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# inspect_tensors diagnostic tool
inspect_tensors: inspect_tensors.cpp bootstrap.hpp projection.hpp tensor_shuffle.h
	$(CXX) inspect_tensors.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# letter_filter_bench: semantics verification + benchmark for the
# divergent/finite letter support filters (solve_collinear.hpp)
letter_filter_bench: letter_filter_bench.cpp bootstrap.hpp projection.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp
	$(CXX) letter_filter_bench.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# regression gate: re-run every recorded baseline (see
# audits/baseline-*/manifest*.json) in hermetic sandboxes and compare CRC32s
# of every output. Run this after ANY change to the calculation core —
# outputs must stay bit-identical.
regression:
	@for m in audits/baseline-*/manifest*.json; do \
		echo "== regression: $$m"; \
		python3 scripts/regression_check.py $$m || exit 1; \
	done

# aggregate local gate: regression baselines + front-end registry sync.
# (CI runs the machine-independent subset of this — see .github/workflows.)
check: regression
	python3 front-end/web/scripts/rhs-modes-sync.py

# bench: micro-benchmarks (crc_bench: per-byte istreambuf CRC vs chunked
# buffer CRC; letter_filter_bench: letter support-filter semantics + timing).
bench: bench/crc_bench
	./bench/crc_bench $(or $(BENCH_FILE),output/FEC_6.wxf)

bench/crc_bench: bench/crc_bench.cpp
	$(CXX) bench/crc_bench.cpp -O3 -std=c++20 -o $@

# clean target
clean:
	rm -f bootstrap bootstrap.exe compute_rhs inspect_tensors tensor_add tensor_ops letter_filter_bench bench/crc_bench
