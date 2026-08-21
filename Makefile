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
CXXFLAGS := -O3 -flto -mcpu=native -mtune=native -std=c++20 -I.
LDLIBS := -lflint -lgmp -lmimalloc

ifeq ($(OS),Windows_NT)
EXE := bootstrap.exe
LDLIBS += -ltbb12
else
EXE := bootstrap
LDLIBS += -ltbb
endif

# target
all: $(EXE)

# executable file
$(EXE): bootstrap.cpp bootstrap.hpp projection.hpp solve_symmetry.hpp solve_collinear.hpp linear_solve.hpp tensor_expand.hpp tensor_shuffle.h
	$(CXX) bootstrap.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# compute_rhs executable (standalone RHS computation module)
compute_rhs: compute_rhs.cpp compute_rhs.hpp tensor_shuffle.h bootstrap.hpp projection.hpp solve_collinear.hpp linear_solve.hpp tensor_expand.hpp
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

# clean target
clean:
	rm -f bootstrap bootstrap.exe compute_rhs inspect_tensors tensor_add tensor_ops
