.DEFAULT_GOAL := all
-include config.mk
CXX ?= g++
ARCH_FLAGS ?=
CXXFLAGS ?= -O3 -std=c++20 $(ARCH_FLAGS)
CPPFLAGS += -Isrc -Isrc/reduced_d3 -Isrc/SparseRREF
$(foreach lib,FLINT MPFR GMP MIMALLOC TBB,$(if $($(lib)_DIR),$(eval CPPFLAGS += -I$($(lib)_DIR)/include)$(eval LDFLAGS += -L$($(lib)_DIR)/lib -L$($(lib)_DIR)/lib64)))
LDLIBS += -lflint -lmpfr -lgmp -lmimalloc
ifeq ($(OS),Windows_NT)
EXE := .exe
LDLIBS += -ltbb12 -lpsapi
else
LDLIBS += -ltbb -pthread
endif
PROGRAMS := irr_fec irr_lec irr_sew irr_convert irr_conditions irr_pipeline irr_cache
TARGETS := $(addprefix bin/,$(addsuffix $(EXE),$(PROGRAMS)))
.PHONY: all test test-transpose test-parallel test-ope-seed
all: $(TARGETS)
bin/irr_conditions_v2$(EXE): bin/irr_conditions$(EXE)
	cp $< $@
bin/%$(EXE): src/reduced_d3/%.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
-include $(wildcard bin/*.d)
bin/test_reconstruction_reuse$(EXE): tests/reconstruction_reuse.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/pipeline_compare$(EXE): tests/pipeline_compare.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/test_stage_options$(EXE): tests/stage_options.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/test_ope_rank$(EXE): tests/ope_rank.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/test_recurrence_store$(EXE): tests/recurrence_store.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
.PHONY: test-ope-library
test-ope-library:
	python3 tests/ope_library.py
.PHONY: test-reconstruction-reuse
test-reconstruction-reuse: bin/test_reconstruction_reuse$(EXE)
	@root="$${TEST_REUSE_ROOT:-output/reuse_test_$$(date +%Y%m%d_%H%M%S)_$$$$}"; ./bin/test_reconstruction_reuse$(EXE) "$$root" && ./bin/test_reconstruction_reuse$(EXE) "$$root" --resume
bin/irr_tensor_compare$(EXE): tests/irr_tensor_compare.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/reencode_wxf$(EXE): tools/reencode_wxf.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
bin/test_shift_cache$(EXE): tests/shift_cache.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
test: all bin/irr_tensor_compare$(EXE) test-transpose test-parallel test-ope-seed test-stream-hash test-block-cache
	bash tests/smoke.sh
	bash tests/physics.sh

bin/test_stream_hash$(EXE): tests/stream_hash.cpp src/reduced_d3/hex_hash.hpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 $< -o $@

.PHONY: test-stream-hash test-block-cache
test-stream-hash: bin/test_stream_hash$(EXE)
	python3 tests/stream_hash.py ./bin/test_stream_hash$(EXE) "$${TEST_HASH_ROOT:-output/hash_test_$$(date +%Y%m%d_%H%M%S)_$$$$}"

test-block-cache: bin/irr_fec$(EXE) bin/irr_cache$(EXE)
	python3 tests/block_cache.py "$${TEST_BLOCK_ROOT:-output/block_test_$$(date +%Y%m%d_%H%M%S)_$$$$}" --exe ./bin/irr_fec$(EXE)

bin/test_sparse_mat_transpose$(EXE): tests/sparse_mat_transpose.cpp src/SparseRREF/sparse_type.h
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
test-transpose: bin/test_sparse_mat_transpose$(EXE)
	./bin/test_sparse_mat_transpose$(EXE)

bin/irr_physics_compare$(EXE): tests/irr_physics_compare.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

bin/test_parallel$(EXE): tests/parallel.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)

test-parallel: bin/test_parallel$(EXE)
	./bin/test_parallel$(EXE) "$${TEST_PARALLEL_ROOT:-output/parallel_$$(date +%Y%m%d_%H%M%S)_$$$$}"

bin/test_ope_seed_parallel$(EXE): tests/ope_seed_parallel.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)

test-ope-seed: bin/test_ope_seed_parallel$(EXE)
	./bin/test_ope_seed_parallel$(EXE)

bin/test_nmhv_ope_germ$(EXE): tests/nmhv_ope_germ.cpp src/nmhv_ope_germ.hpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 $< -o $@ $(LDFLAGS) $(LDLIBS)
.PHONY: test-nmhv-ope
test-nmhv-ope: bin/test_nmhv_ope_germ$(EXE)
	./bin/test_nmhv_ope_germ$(EXE) data

bin/test_nmhv_ope_rows$(EXE): tests/nmhv_ope_rows.cpp src/nmhv_ope_germ.hpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) -O1 -std=c++20 $< -o $@ $(LDFLAGS) $(LDLIBS)
.PHONY: test-nmhv-ope-rows
test-nmhv-ope-rows: bin/test_nmhv_ope_rows$(EXE)
	./bin/test_nmhv_ope_rows$(EXE) data

bin/test_parallel_projection$(EXE): tests/parallel_projection.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)
