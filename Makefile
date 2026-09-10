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
else ifeq ($(shell uname -m),aarch64)
ARCH_FLAGS := -mcpu=native -mtune=native
endif
ifeq ($(PORTABLE),1)
ARCH_FLAGS :=
endif
CXXFLAGS := -O3 -flto $(ARCH_FLAGS) -std=c++20 -I.

# mimalloc is a performance add-on (it overrides the GMP/FLINT allocators
# when USE_MIMALLOC=1); HPC boxes often lack it, so link it only when a copy
# exists in the usual library locations. The math is identical without it.
# Override explicitly with MIMALLOC=1 (force) or MIMALLOC=0 (forbid).
MIMALLOC_LIB := $(firstword $(wildcard /opt/homebrew/lib/libmimalloc.* /usr/local/lib/libmimalloc.* /usr/lib/*/libmimalloc.* /usr/lib/libmimalloc.*))
MIMALLOC_HEADER := $(firstword $(wildcard /opt/homebrew/include/mimalloc.h /usr/local/include/mimalloc.h /usr/include/mimalloc.h))
USE_ALLOCATOR := $(if $(filter 1,$(MIMALLOC)),1,$(if $(filter 0,$(MIMALLOC)),,$(if $(MIMALLOC_LIB),$(if $(MIMALLOC_HEADER),1))))
ALLOCATOR_FLAGS := $(if $(USE_ALLOCATOR),-DUSE_MIMALLOC,-UUSE_MIMALLOC)
LDLIBS := -lflint -lgmp $(if $(USE_ALLOCATOR),-lmimalloc)

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
all: $(EXE) compute_rhs tensor_add tensor_ops

CORE_HEADERS := $(wildcard *.hpp *.h SparseRREF/*.hpp SparseRREF/*.h)
$(EXE) compute_rhs tensor_add tensor_ops inspect_tensors letter_filter_bench tests/numerical_probe tests/wxf_probe: $(CORE_HEADERS)

# executable file
$(EXE): bootstrap.cpp bootstrap.hpp projection.hpp solve_symmetry.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp tensor_shuffle.h
	$(CXX) bootstrap.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# compute_rhs executable (standalone RHS computation module)
compute_rhs: compute_rhs.cpp compute_rhs.hpp tensor_shuffle.h bootstrap.hpp projection.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp
	$(CXX) compute_rhs.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# tensor_add executable (weighted sum of sparse tensors)
tensor_add: tensor_add.cpp tensor_shuffle.h numeric_parse.hpp
	$(CXX) tensor_add.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# tensor_ops executable (ternary contraction, matrix power, tensor join, squeeze)
tensor_ops: tensor_ops.cpp tensor_shuffle.h numeric_parse.hpp
	$(CXX) tensor_ops.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# inspect_tensors diagnostic tool
inspect_tensors: inspect_tensors.cpp bootstrap.hpp projection.hpp tensor_shuffle.h
	$(CXX) inspect_tensors.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# letter_filter_bench: semantics verification + benchmark for the
# divergent/finite letter support filters (solve_collinear.hpp)
letter_filter_bench: letter_filter_bench.cpp bootstrap.hpp projection.hpp solve_collinear.hpp linear_solve.hpp incremental_solve.hpp tensor_expand.hpp
	$(CXX) letter_filter_bench.cpp -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

# regression gate: re-run every recorded baseline (see
# tests/baselines/*/manifest*.json) in hermetic sandboxes and compare CRC32s
# of every output. Run this after ANY change to the calculation core —
# outputs must stay bit-identical.
regression:
	@for m in tests/baselines/*/manifest*.json; do \
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
	rm -f bootstrap bootstrap.exe compute_rhs inspect_tensors tensor_add tensor_ops letter_filter_bench bench/crc_bench tests/numerical_probe tests/wxf_probe tests/wxf_probe_sanitized

# Public, self-contained exact regression fixtures (no private project data).
tests/numerical_probe: tests/numerical_probe.cpp bootstrap.hpp linear_solve.hpp incremental_solve.hpp projection.hpp
	$(CXX) $< -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

check-public: bootstrap compute_rhs tensor_add tensor_ops tests/numerical_probe check-wxf check-cache check-templates
	python3 tests/native_regression.py
	python3 front-end/web/scripts/rhs-modes-sync.py

# Checked WXF decoding, including a bounded malformed-file corpus.
tests/wxf_probe: tests/wxf_probe.cpp
	$(CXX) $< -o $@ $(CPPFLAGS) $(CXXFLAGS) $(ALLOCATOR_FLAGS) $(LDLIBS)

check-wxf: tests/wxf_probe
	python3 tests/wxf_regression.py

# Keep sanitizer instrumentation on the reader, without an allocator override.
check-wxf-sanitized: tests/wxf_probe.cpp $(CORE_HEADERS)
	$(CXX) $< -o tests/wxf_probe_sanitized $(CPPFLAGS) $(filter-out -flto -O%,$(CXXFLAGS)) -O1 -g -UUSE_MIMALLOC -fsanitize=address,undefined -fno-omit-frame-pointer $(LDLIBS)
	python3 tests/wxf_regression.py --probe tests/wxf_probe_sanitized

.PHONY: all clean regression check check-public check-wxf check-wxf-sanitized bench

tests/cache_probe: tests/cache_probe.cpp native_cache.hpp
	$(CXX) $< -o $@ $(CPPFLAGS) $(CXXFLAGS)

check-cache: all tests/cache_probe
	python3 tests/cache_regression.py

.PHONY: check-cache

check-templates: all tests/numerical_probe
	python3 tests/template_regression.py

.PHONY: check-templates
