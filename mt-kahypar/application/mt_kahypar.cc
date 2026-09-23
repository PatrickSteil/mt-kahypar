/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <iostream>
#include <chrono>
#include <algorithm>
#include <exception>
#include <vector>

#include "include/lib_generic_impls.h"
#include "include/lib_helper_functions.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/parallel/thread_management.h"
#include "mt-kahypar/partition/partitioner_facade.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/registries/registry.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/exception.h"

using namespace mt_kahypar;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

// ! Validates the initial partition and determines whether it is a k-way partition
// ! (all IDs < k) or a fragment clustering (e.g., computed by PUNCH). Returns the
// ! number of blocks required to store it in a partitioned hypergraph.
PartitionID checkInitialPartition(const std::vector<PartitionID>& initial_partition,
                                  Context& context,
                                  const mt_kahypar_hypergraph_t hypergraph) {
  PartitionID max_id = -1;
  for ( const PartitionID id : initial_partition ) {
    if ( id < 0 ) {
      throw InvalidInputException("Initial partition file contains negative block IDs!");
    }
    max_id = std::max(max_id, id);
  }
  context.partition.initial_partition_is_kway = max_id < context.partition.k;

  if ( context.partition.mode != Mode::direct || context.isNLevelPartitioning() ) {
    throw InvalidParameterException(
      "Initial partitions are only supported for multilevel presets in direct mode (e.g., --preset-type=default)!");
  }
  if ( context.partition.fixed_vertex_filename != "" ) {
    throw InvalidParameterException("Initial partitions can not be combined with fixed vertices!");
  }
  if ( !context.partition.initial_partition_is_kway ) {
    // A partitioned hypergraph with one block per fragment requires O(m * #fragments) memory
    // for hypergraphs, while partitioned graphs only store O(#fragments) block weights.
    if ( hypergraph.type != STATIC_GRAPH && hypergraph.type != DYNAMIC_GRAPH ) {
      throw InvalidParameterException("Initial partitions with more than k blocks are only supported for graphs!");
    }
    if ( context.partition.objective == Objective::steiner_tree ) {
      throw InvalidParameterException(
        "Initial partitions with more than k blocks are not supported for the steiner_tree objective!");
    }
  }
  return std::max(context.partition.k, max_id + 1);
}

int run(int argc, char* argv[]) {
  Context context(false);
  processCommandLineInput(context, argc, argv);

  if ( context.partition.preset_type == PresetType::UNDEFINED ) {
    ERR("No preset specified (--preset-type)");
  }

  // Determine instance (graph or hypergraph) and partition type
  if ( context.partition.instance_type == InstanceType::UNDEFINED ) {
    context.partition.instance_type = to_instance_type(context.partition.file_format);
  }
  context.partition.partition_type = to_partition_c_type(
    context.partition.preset_type, context.partition.instance_type);


  context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
  if (context.partition.enable_logging) {
    io::printBanner();
  }

  utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    utils::Randomize::instance().enableLocalizedParallelShuffle(
      context.shared_memory.shuffle_block_size);
  }

  if constexpr (parallel::provides_hardware_information) {
    size_t num_available_cpus = parallel::num_hardware_cpus();
    if ( num_available_cpus < context.shared_memory.num_threads ) {
      WARNING("There are currently only " << num_available_cpus << " cpus available. "
        << "Setting number of threads from " << context.shared_memory.num_threads
        << " to " << num_available_cpus);
      context.shared_memory.num_threads = num_available_cpus;
    }
  }

  // Initialize TBB task arenas on numa nodes
  parallel::initialize_tbb(context.shared_memory.num_threads);

  if constexpr (parallel::provides_hardware_information) {
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    parallel::activate_interleaved_membind_policy();
  }

  // Read Hypergraph
  utils::Timer& timer =
    utils::Utilities::instance().getTimer(context.utility_id);
  timer.start_timer("io_hypergraph", "I/O Hypergraph");
  mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
      context.partition.graph_filename, context.partition.preset_type,
      context.partition.instance_type, context.partition.file_format,
      context.preprocessing.stable_construction_of_incident_edges,
      /*remove_single_pin_hes=*/true, /*print_warnings=*/true);
  timer.stop_timer("io_hypergraph");

  // Read Target Graph
  std::unique_ptr<TargetGraph> target_graph;
  if ( context.partition.objective == Objective::steiner_tree ) {
    if ( context.mapping.target_graph_file != "" ) {
      target_graph = std::make_unique<TargetGraph>(
        io::readInputFile<ds::StaticGraph>(
          context.mapping.target_graph_file, FileFormat::Metis,
          /*stable_construnction=*/true, /*remove_single_pin_hes=*/true, /*print_warnings=*/true));
    } else {
      throw InvalidInputException("No target graph file specified (use -g <file> or --target-graph=<file>)!");
    }
  }

  if ( context.partition.fixed_vertex_filename != "" ) {
    timer.start_timer("read_fixed_vertices", "Read Fixed Vertex File");
    io::addFixedVerticesFromFile(hypergraph,
      context.partition.fixed_vertex_filename, context.partition.k);
    timer.stop_timer("read_fixed_vertices");
  }

  // Initialize Memory Pool and Algorithm/Policy Registries
  register_memory_pool(hypergraph, context);
  register_algorithms_and_policies();

  // Read Initial Partition
  mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph { nullptr, NULLPTR_PARTITION };
  if ( context.partition.initial_partition_filename != "" ) {
    timer.start_timer("read_initial_partition", "Read Initial Partition File");
    std::vector<PartitionID> initial_partition;
    io::readPartitionFile(context.partition.initial_partition_filename,
      lib::num_nodes<false>(hypergraph), initial_partition);
    timer.stop_timer("read_initial_partition");
    const PartitionID num_blocks = checkInitialPartition(initial_partition, context, hypergraph);
    partitioned_hypergraph = lib::create_partitioned_hypergraph(
      hypergraph, context, num_blocks, initial_partition.data());
  }

  // Partition Hypergraph
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  if ( context.partition.initial_partition_filename != "" ) {
    // Start from the given partition: a single V-cycle that either uses the
    // input as initial solution (k-way) or only restricts coarsening (fragments)
    context.partition.num_vcycles = std::max(context.partition.num_vcycles, UL(1));
    PartitionerFacade::improve(partitioned_hypergraph, context, target_graph.get());
  } else {
    partitioned_hypergraph = PartitionerFacade::partition(hypergraph, context, target_graph.get());
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  PartitionerFacade::printPartitioningResults(
    partitioned_hypergraph, context, elapsed_seconds);

  if ( context.partition.sp_process_output ) {
    std::cout << PartitionerFacade::serializeResultLine(
      partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if ( context.partition.csv_output ) {
    std::cout << PartitionerFacade::serializeCSV(
      partitioned_hypergraph, context, elapsed_seconds) << std::endl;
  }

  if (context.partition.write_partition_file) {
    PartitionerFacade::writePartitionFile(
      partitioned_hypergraph, context.partition.graph_partition_filename);
  }

  parallel::MemoryPool::instance().free_memory_chunks();
  parallel::terminate_tbb();

  utils::delete_hypergraph(hypergraph);
  utils::delete_partitioned_hypergraph(partitioned_hypergraph);

  return 0;
}

int main(int argc, char* argv[]) {
#ifdef NDEBUG
  try {
    return run(argc, argv);
  } catch (const InvalidInputException& e) {
    std::cerr << "\n " << e.what() << std::endl;
    return 1;
  } catch (const InvalidParameterException& e) {
    std::cerr << "\n" << e.what() << std::endl;
    return 1;
  } catch (const UnsupportedOperationException& e) {
    std::cerr << "\n" << e.what() << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
    return 1;
  }
#else
  return run(argc, argv);
#endif
}