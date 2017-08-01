#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>

#include "succinct/mapper.hpp"

#include "configuration.hpp"
#include "util.hpp"
#include "verify_collection.hpp"
#include "index_build_utils.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

#include "clustered_index_types.hpp"
#include "clustered_binary_freq_collection.hpp"
#include "partitioned_sequence.hpp"
#include "reference_selector.hpp"

using ds2i::logger;

void create_clustered_collection(ds2i::clustered_binary_freq_collection& input,
                                 ds2i::global_parameters const& params,
                                 const char* cluster_filename,
                                 const char* output_filename,
                                 bool check,
                                 uint32_t MAX_REF_SIZE,
                                 uint32_t divisor)
{
    using namespace ds2i;
    uint64_t universe = input.num_docs();
    logger() << "Processing " << universe << " documents" << std::endl;
    double tick = get_time_usecs();
    double user_tick = get_user_time_usecs();
    clustered_opt_index::builder builder(universe, params);
    progress_logger plog;

    std::ifstream file(cluster_filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    std::istream instream(&inbuf);

    std::string line;
    uint32_t c = 1;
    std::vector<sequence_t> partitions;
    std::vector<pos_t> docs_positions;
    std::vector<pos_t> freqs_positions;
    while (std::getline(instream, line))
    {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "));
        uint32_t cluster_size = std::atoi(data.front().data());

        std::cout << "cluster_size: " << cluster_size << std::endl;

        docs_positions.reserve(cluster_size);
        freqs_positions.reserve(cluster_size);

        std::for_each(++data.begin(), data.end(),
            [&docs_positions, &freqs_positions](const std::string& s)
            {
                pos_t p = std::stoull(s.data());
                docs_positions.push_back(p);
                freqs_positions.push_back(p - 2);
            });

        input.add_positions(docs_positions, freqs_positions);
        input.set_positions();

        partitions.clear();
        partitions.reserve(cluster_size);

        for (auto const& plist : input)
        {
            auto const& docs = plist.docs;
            partitions.push_back(
                std::move(
                    ds2i::partitioned_sequence<>::compute_partition(docs.begin(),
                                                                    universe,
                                                                    docs.size(),
                                                                    params)
                        )
                );
        }

        reference_selector rs(input, partitions, params, cluster_size, MAX_REF_SIZE);
        sequence_t reference = rs.select_reference(divisor);

        std::cout << "reference size: " << reference.size() << std::endl;

        builder.add_reference(reference.size(), reference.begin());
        uint32_t i = 0;
        for (auto const& plist : input)
        {
            uint64_t freqs_sum =
                std::accumulate(plist.freqs.begin(),
                                plist.freqs.end(),
                                uint64_t(0));
            builder.add_posting_list(plist.docs,
                                     plist.freqs,
                                     freqs_sum,
                                     partitions[i++]);
            plog.done_sequence(plist.docs.size());
        }
        logger() << "cluster-" << c++ << " encoded" << std::endl;

        docs_positions.clear();
        freqs_positions.clear();
    }

    file.close();
    plog.log();
    clustered_opt_index coll;
    builder.build(coll);
    double elapsed_secs = (get_time_usecs() - tick) / 1000000;
    double user_elapsed_secs = (get_user_time_usecs() - user_tick) / 1000000;
    logger() << "clustered_opt collection built in "
             << elapsed_secs << " seconds" << std::endl;

    stats_line()
        ("type", "clustered_opt_greedy")
        ("worker_threads", configuration::get().worker_threads)
        ("construction_time", elapsed_secs)
        ("construction_user_time", user_elapsed_secs)
        ;

    dump_stats(coll, "clustered_opt", plog.postings);

    if (output_filename) {
        succinct::mapper::freeze(coll, output_filename);
        if (check) {
            verify_clustered_collection(input, output_filename);
        }
    }
}

int main(int argc, const char** argv)
{
    using namespace ds2i;

    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <collection basename> <cluster filename> <MAX_REF_SIZE> <divisor> [<output filename>] [--check]"
                  << std::endl;
        std::cerr << "<divisor> represents the amount by which the number of candidates is divided to define the length of an apoch"
                  << std::endl;
        return 1;
    }

    const char* input_basename = argv[1];
    const char* cluster_filename = argv[2];
    const uint32_t MAX_REF_SIZE = std::atoi(argv[3]);
    const uint32_t divisor = std::atoi(argv[4]);

    const char* output_filename = nullptr;
    if (argc > 5) {
        output_filename = argv[5];
    }

    bool check = false;
    if (argc > 6 && std::string(argv[6]) == "--check") {
        check = true;
    }

    clustered_binary_freq_collection input(input_basename, check);
    global_parameters params;
    params.log_partition_size = configuration::get().log_partition_size;
    create_clustered_collection(
        input, params, cluster_filename,
        output_filename, check, MAX_REF_SIZE, divisor
    );

    return 0;
}
