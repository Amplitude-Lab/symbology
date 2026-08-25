#include "bootstrap.hpp"
#include "projection.hpp"
#include "solve_symmetry.hpp"
#include "solve_collinear.hpp"

#include <cstdlib>

using index_t = int32_t;
using scalar_t = rat_t;

enum class bootstrap_mode_t {
	none,
	extend,
	sew,
	induce,
	project,
	solve_symmetry,
	solve_collinear
};

struct args_t {
	bootstrap_mode_t mode = bootstrap_mode_t::none;
	std::filesystem::path condition;
	std::filesystem::path first;
	std::filesystem::path last;
	std::filesystem::path transform;
	std::filesystem::path output;
	std::string symmetry;   // --project: collinear | cyclic | flip | parity
	std::string target;     // --project/--solve-*: e.g. SEW_5p1
	// --solve-collinear:
	std::filesystem::path rhs;
	std::string projection_type;  // "finite" or "divergent"
	std::vector<std::filesystem::path> basis_paths;  // expansion bases (highest weight first)
	std::string letter_projection;  // "identity" or a file path (resolved against exec dir)
	std::string solver = "incremental"; // --solve-collinear: "incremental" (default) or "sampled"
	// --data-dir / --output-dir (for --project, --solve-symmetry, --solve-collinear).
	// Empty means "use default" (resolved against the executable directory after parsing).
	std::filesystem::path data_dir;
	std::filesystem::path output_dir;
};

// bootstrap.cpp owns only the command-line contract and WXF I/O.
// The algebraic steps are kept in bootstrap.hpp so this file stays close to a dispatcher.
void print_usage(const char* program) {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "  " << program << " --extend -c <condition.wxf> -f <FEC_in.wxf> -o <FEC_out.wxf>" << std::endl;
	std::cerr << "  " << program << " --extend -c <condition.wxf> -l <LEC_in.wxf> -o <LEC_out.wxf>" << std::endl;
	std::cerr << "  " << program << " --sew -c <condition.wxf> -f <FEC.wxf> -l <LEC.wxf> -o <SEW.wxf>" << std::endl;
	std::cerr << "  " << program << " --project --symmetry <collinear|cyclic|flip|parity> --target <SEW_FpL|FEC_W|LEC_W> [--data-dir <dir>] [--output-dir <dir>]" << std::endl;
	std::cerr << "  " << program << " --solve-symmetry --symmetry <collinear|cyclic|flip|parity> --target <SEW_FpL|FEC_W|LEC_W> [--data-dir <dir>] [--output-dir <dir>] (note: collinear projections are usually non-square and will be rejected by the solver; prefer --solve-collinear)" << std::endl;
	std::cerr << "  " << program << " --solve-collinear --target <SEW_FpL|FEC_W> --rhs <rhs.wxf> --projection <finite|divergent> --letter-projection <file|identity> [--basis <basis.wxf> ...] [--solver <incremental|sampled>] [--data-dir <dir>] [--output-dir <dir>]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --data-dir defaults to <exec_dir>/data; --output-dir defaults to <exec_dir>/output." << std::endl;
	std::cerr << "  --extend / --sew ignore --data-dir / --output-dir (use explicit -c/-f/-l/-o paths)." << std::endl;
}

std::string take_value(int& i, int argc, char* argv[], const std::string& flag) {
	if (i + 1 >= argc) {
		throw std::runtime_error("Missing value after " + flag);
	}
	i++;
	return argv[i];
}

void set_mode(args_t& args, bootstrap_mode_t mode) {
	if (args.mode != bootstrap_mode_t::none) {
		throw std::runtime_error("Only one mode can be specified.");
	}
	args.mode = mode;
}

args_t parse_args(int argc, char* argv[]) {
	args_t args;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--extend") {
			set_mode(args, bootstrap_mode_t::extend);
		}
		else if (arg == "--sew") {
			set_mode(args, bootstrap_mode_t::sew);
		}
		else if (arg == "--induce") {
			set_mode(args, bootstrap_mode_t::induce);
		}
		else if (arg == "--project") {
			set_mode(args, bootstrap_mode_t::project);
		}
		else if (arg == "--solve-symmetry") {
			set_mode(args, bootstrap_mode_t::solve_symmetry);
		}
		else if (arg == "--solve-collinear") {
			set_mode(args, bootstrap_mode_t::solve_collinear);
		}
		else if (arg == "--symmetry") {
			args.symmetry = take_value(i, argc, argv, arg);
		}
		else if (arg == "--target") {
			args.target = take_value(i, argc, argv, arg);
		}
		else if (arg == "--rhs") {
			args.rhs = take_value(i, argc, argv, arg);
		}
		else if (arg == "--projection") {
			args.projection_type = take_value(i, argc, argv, arg);
		}
		else if (arg == "--letter-projection") {
			args.letter_projection = take_value(i, argc, argv, arg);
		}
		else if (arg == "--solver") {
			args.solver = take_value(i, argc, argv, arg);
		}
		else if (arg == "--basis") {
			args.basis_paths.push_back(take_value(i, argc, argv, arg));
		}
		else if (arg == "-c" || arg == "--condition") {
			args.condition = take_value(i, argc, argv, arg);
		}
		else if (arg == "-f" || arg == "--first") {
			args.first = take_value(i, argc, argv, arg);
		}
		else if (arg == "-l" || arg == "--last") {
			args.last = take_value(i, argc, argv, arg);
		}
		else if (arg == "-t" || arg == "--transform") {
			args.transform = take_value(i, argc, argv, arg);
		}
		else if (arg == "-o" || arg == "--output") {
			args.output = take_value(i, argc, argv, arg);
		}
		else if (arg == "--data-dir") {
			args.data_dir = take_value(i, argc, argv, arg);
		}
		else if (arg == "--output-dir") {
			args.output_dir = take_value(i, argc, argv, arg);
		}
		else if (arg == "-h" || arg == "--help") {
			print_usage(argv[0]);
			std::exit(0);
		}
		else {
			throw std::runtime_error("Unknown argument: " + arg);
		}
	}
	return args;
}

void validate_args(const args_t& args) {
	if (args.mode == bootstrap_mode_t::none) {
		throw std::runtime_error("Missing mode: use --extend, --sew, --induce, or --project.");
	}
	if (args.mode == bootstrap_mode_t::induce) {
		throw std::runtime_error("--induce is reserved for a future workflow stage.");
	}
	if (args.mode == bootstrap_mode_t::project) {
		if (args.symmetry.empty()) {
			throw std::runtime_error("--project requires --symmetry <collinear|cyclic|flip|parity>.");
		}
		if (args.target.empty()) {
			throw std::runtime_error("--project requires --target <SEW_FpL> (e.g. SEW_5p1).");
		}
		// Validate symmetry name early via get_symmetry_info (throws on unknown).
		get_symmetry_info(args.symmetry);
		// Validate target name early via parse_target (throws on bad format).
		parse_target(args.target);
		return;
	}
	if (args.mode == bootstrap_mode_t::solve_symmetry) {
		if (args.symmetry.empty()) {
			throw std::runtime_error("--solve-symmetry requires --symmetry <cyclic|flip|parity>.");
		}
		if (args.target.empty()) {
			throw std::runtime_error("--solve-symmetry requires --target <SEW_FpL> (e.g. SEW_5p1).");
		}
		get_symmetry_info(args.symmetry);
		parse_target(args.target);
		return;
	}
	if (args.mode == bootstrap_mode_t::solve_collinear) {
		if (args.target.empty()) {
			throw std::runtime_error("--solve-collinear requires --target <SEW_FpL|FEC_W> (e.g. SEW_5p1).");
		}
		if (args.rhs.empty()) {
			// Q7: exit cleanly (do not throw) when --rhs is missing.
			std::cerr << "Error: --solve-collinear requires --rhs <rhs.wxf>." << std::endl;
			std::cerr << "   The RHS is the collinear boundary expression (rank k, dims (11,...,11))." << std::endl;
			std::cerr << "   Use '--rhs 0' for an empty RHS (all-zero boundary)." << std::endl;
			std::cerr << "   The RHS can be computed by the compute_rhs module, or provided directly." << std::endl;
			std::exit(1);
		}
		if (args.projection_type.empty()) {
			throw std::runtime_error("--solve-collinear requires --projection <finite|divergent>.");
		}
		if (args.projection_type != "finite" && args.projection_type != "divergent") {
			throw std::runtime_error("--projection must be 'finite' or 'divergent', got: " + args.projection_type);
		}
		if (args.letter_projection.empty()) {
			std::cerr << "Error: --solve-collinear requires --letter-projection <file|identity>." << std::endl;
			std::cerr << "   Use '--letter-projection identity' to skip projection (solve in full letter space)." << std::endl;
			std::cerr << "   Use '--letter-projection <file>' to project each letter slot, e.g." << std::endl;
			std::cerr << "       --letter-projection output/collinear/colprojdiv_w1.wxf" << std::endl;
			std::exit(1);
		}
		if (args.solver != "sampled" && args.solver != "incremental") {
			throw std::runtime_error("--solver must be 'sampled' or 'incremental', got: " + args.solver);
		}
		// Validate target name early via parse_target (throws on bad format).
		parse_target(args.target);
		return;
	}
	if (args.condition.empty()) {
		throw std::runtime_error("Missing condition file: use -c/--condition.");
	}
	if (args.output.empty()) {
		throw std::runtime_error("Missing output file: use -o/--output.");
	}
	if (args.mode == bootstrap_mode_t::extend) {
		const bool has_first = !args.first.empty();
		const bool has_last = !args.last.empty();
		if (has_first == has_last) {
			throw std::runtime_error("--extend requires exactly one of -f/--first or -l/--last.");
		}
	}
	if (args.mode == bootstrap_mode_t::sew) {
		if (args.first.empty() || args.last.empty()) {
			throw std::runtime_error("--sew requires both -f/--first and -l/--last.");
		}
	}
	if (!args.transform.empty()) {
		throw std::runtime_error("-t/--transform is only valid for the future --induce mode.");
	}
}

// Relative paths are resolved against the executable directory. This keeps
// "symbology/bootstrap.exe -c data/..." working from the repository root.
std::filesystem::path resolve_path(const std::filesystem::path& base, const std::filesystem::path& path) {
	if (path.is_absolute()) {
		return path;
	}
	return base / path;
}

// Files are read as CSR tensors and moved into bootstrap.hpp, where each
// computation converts to COO only for the operations that need it.
sparse_tensor<scalar_t, index_t, SPARSE_CSR> read_tensor(
	const std::filesystem::path& path,
	const field_t& F,
	thread_pool* pool) {
	std::cout << "Reading file " << path.string() << " ..." << std::endl;
	auto tensor = sparse_tensor_read_wxf<scalar_t, index_t>(path, F, pool);
	print_crc32(path.string(), path);
	return tensor;
}

// SparseRREF emits its native WXF representation here. Use wxf_roundtrip.wls
// when byte-for-byte Mathematica-exported WXF is needed for archive comparison.
void write_tensor(
	const std::filesystem::path& path,
	sparse_tensor<scalar_t, index_t, SPARSE_CSR>&& tensor) {
	if (!path.parent_path().empty()) {
		std::filesystem::create_directories(path.parent_path());
	}
	std::cout << "Writing file " << path.string() << " ..." << std::endl;
	Timer timer;
	timer.start();
	auto u8arr = sparse_tensor_write_wxf(tensor);
	std::ofstream ofs(path, std::ios::binary);
	if (!ofs) {
		throw std::runtime_error("Cannot write file: " + path.string());
	}
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs.close();
	u8arr.clear();
	u8arr.shrink_to_fit();
	timer.stop();
	std::cout << "** Write time: " << timer.milliseconds() << " ms" << std::endl;
	print_crc32(path.string(), path);
}

int main(int argc, char* argv[]) {
	try {
		field_t F(FIELD_QQ);

		args_t args = parse_args(argc, argv);
		validate_args(args);

		// reset() asks SparseRREF to choose the thread count at runtime.
		rref_option_t opt;
		opt->method = 0;
		opt->verbose = true;
		opt->pool.reset();
		thread_pool* pool = &(opt->pool);

		std::filesystem::path base = std::filesystem::path(argv[0]).parent_path();
		if (base.empty()) {
			base = std::filesystem::current_path();
		}

		std::cout << "threads: " << pool->get_thread_count() << std::endl;
		auto local_time = std::chrono::zoned_time{
			std::chrono::current_zone(),
			std::chrono::system_clock::now()
		};
		std::cout << "Task begin at: " << std::format("{:%Y-%m-%d %H:%M:%S %z}", local_time) << std::endl;

		Timer total_timer;
		total_timer.start();

		if (args.mode == bootstrap_mode_t::project) {
			// --project runs its own recursive pipeline (no condition/first/last/output).
			// Resolve data_dir / output_dir with defaults base/"data", base/"output".
			auto data_dir = args.data_dir.empty() ? (base / "data") : resolve_path(base, args.data_dir);
			auto output_dir = args.output_dir.empty() ? (base / "output") : resolve_path(base, args.output_dir);
			run_projection_pipeline<scalar_t, index_t>(
				args.symmetry, args.target, data_dir, output_dir, F, opt);
		}
		else if (args.mode == bootstrap_mode_t::solve_symmetry) {
			// --solve-symmetry: compute invariant space of the target's projection.
			auto data_dir = args.data_dir.empty() ? (base / "data") : resolve_path(base, args.data_dir);
			auto output_dir = args.output_dir.empty() ? (base / "output") : resolve_path(base, args.output_dir);
			run_symmetry_solver<scalar_t, index_t>(
				args.symmetry, args.target, data_dir, output_dir, F, opt);
		}
		else if (args.mode == bootstrap_mode_t::solve_collinear) {
			// --solve-collinear: finite/divergent split + expansion + linear solve.
			auto target = parse_target(args.target);

			auto data_dir = args.data_dir.empty() ? (base / "data") : resolve_path(base, args.data_dir);
			auto output_dir = args.output_dir.empty() ? (base / "output") : resolve_path(base, args.output_dir);
			auto collinear_dir = output_dir / "collinear";

			// Determine target weight and target basis path
			size_t target_weight;
			std::filesystem::path target_basis_path;
			if (target.kind == target_kind_t::SEW) {
				target_weight = target.fec_weight + target.lec_weight;
				target_basis_path = collinear_dir / (target.name + "_basis.wxf");
			} else if (target.kind == target_kind_t::FEC) {
				target_weight = target.fec_weight;
				target_basis_path = collinear_dir / ("first_w" + std::to_string(target.fec_weight) + "_basis.wxf");
			} else {
				throw std::runtime_error("--solve-collinear: LEC targets not yet supported");
			}

			// Auto-detect chain base paths (lowest weight first) for projection computation
			auto chain_base_paths = detect_chain_base_paths(target, output_dir);

			// Expansion bases (highest weight first): either provided by user,
			// or auto-detected from the standard naming convention.
			std::vector<std::filesystem::path> expansion_bases = args.basis_paths;
			if (expansion_bases.empty()) {
				// Auto-detect: weights (target_weight-1) down to 2.
				// Use signed counter because size_t would wrap past 0.
				for (long w = static_cast<long>(target_weight) - 1; w >= 2; w--) {
					expansion_bases.push_back(collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf"));
				}
			}

			// Handle --rhs 0 (empty RHS) vs --rhs <file>
		// Q7: "--rhs 0" means the RHS is an empty (all-zero) tensor. We pass the
		// literal "0" as a sentinel; run_collinear_solver constructs an all-zero b
		// in-memory with shape derived from A (product of A.dims[1..] = n_constraints).
		// This avoids writing a 0-nnz tensor to disk (the WXF writer cannot serialize
		// 0-nnz tensors, and a rank-1 dim-{1} file would fail the dimension check).
		std::filesystem::path rhs_path;
		if (args.rhs == "0") {
			rhs_path = "0";  // sentinel, not a file path
		} else {
			rhs_path = resolve_path(base, args.rhs);
		}

			// Pass SEW target name for SEW-level projection naming
			std::string sew_name;
			if (target.kind == target_kind_t::SEW) {
				sew_name = target.name;
			}

			// Resolve --letter-projection: "identity" means no projection; a path
			// is resolved against the executable directory if relative.
			std::string letter_projection = args.letter_projection;
			if (letter_projection != "identity") {
				letter_projection = resolve_path(base, letter_projection).string();
			}

			run_collinear_solver<scalar_t, index_t>(
				target_basis_path, rhs_path, args.projection_type,
				expansion_bases, chain_base_paths,
				target_weight, data_dir, output_dir, F, opt, sew_name,
			letter_projection, args.solver);
		}
		else {
			std::filesystem::path condition_path = resolve_path(base, args.condition);
			std::filesystem::path first_path = resolve_path(base, args.first);
			std::filesystem::path last_path = resolve_path(base, args.last);
			std::filesystem::path output_path = resolve_path(base, args.output);

			auto dlogmat = read_tensor(condition_path, F, pool);

			// Each mode consumes dlogmat exactly once, matching the moved-tensor API
			// in bootstrap.hpp and avoiding accidental large tensor copies.
			if (args.mode == bootstrap_mode_t::extend && !args.first.empty()) {
				std::cout << "Symbol bootstrap: " << first_path.string() << " -> " << output_path.string()
						  << " @ " << condition_path.string() << std::endl;
				auto FEC = read_tensor(first_path, F, pool);
				auto output = extend_forward(std::move(dlogmat), std::move(FEC), F, opt);
				write_tensor(output_path, std::move(output));
			}
			else if (args.mode == bootstrap_mode_t::extend && !args.last.empty()) {
				std::cout << "Symbol bootstrap: " << last_path.string() << " -> " << output_path.string()
						  << " @ " << condition_path.string() << std::endl;
				auto LEC = read_tensor(last_path, F, pool);
				auto output = extend_backward(std::move(dlogmat), std::move(LEC), F, opt);
				write_tensor(output_path, std::move(output));
			}
			else if (args.mode == bootstrap_mode_t::sew) {
				std::cout << "Symbol bootstrap: " << first_path.string() << " + " << last_path.string()
						  << " -> " << output_path.string() << " @ " << condition_path.string() << std::endl;
				auto FEC = read_tensor(first_path, F, pool);
				auto LEC = read_tensor(last_path, F, pool);
				auto output = sew_first_last(std::move(dlogmat), std::move(FEC), std::move(LEC), F, opt);
				write_tensor(output_path, std::move(output));
			}
		}

		total_timer.stop();
		local_time = std::chrono::zoned_time{
			std::chrono::current_zone(),
			std::chrono::system_clock::now()
		};
		std::cout << "Task end at: " << std::format("{:%Y-%m-%d %H:%M:%S %z}", local_time) << std::endl;
		std::cout << "Done in " << total_timer.milliseconds() << " ms" << std::endl << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		print_usage(argv[0]);
		return 1;
	}
	return 0;
}
