#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>

#include "succinct/mapper.hpp"

#include "configuration.hpp"
#include "index_types.hpp"
#include "util.hpp"
#include "verify_collection.hpp"
#include "index_build_utils.hpp"
#include "binary_freq_collection.hpp"
#include "clustered_binary_freq_collection.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

using ds2i::logger;

template <typename Collection>
void dump_index_specific_stats(Collection const&, std::string const&)
{}


void dump_index_specific_stats(ds2i::uniform_index const& coll,
                               std::string const& type)
{
    ds2i::stats_line()
        ("type", type)
        ("log_partition_size", int(coll.params().log_partition_size))
        ;
}


void dump_index_specific_stats(ds2i::opt_index const& coll,
                               std::string const& type)
{
    auto const& conf = ds2i::configuration::get();

    uint64_t length_threshold = 4096;
    double long_postings = 0;
    double docs_partitions = 0;
    double freqs_partitions = 0;

    for (size_t s = 0; s < coll.size(); ++s) {
        auto const& list = coll[s];
        if (list.size() >= length_threshold) {
            long_postings += list.size();
            docs_partitions += list.docs_enum().num_partitions();
            freqs_partitions += list.freqs_enum().base().num_partitions();
        }
    }

    ds2i::stats_line()
        ("type", type)
        ("eps1", conf.eps1)
        ("eps2", conf.eps2)
        ("fix_cost", conf.fix_cost)
        ("docs_avg_part", long_postings / docs_partitions)
        ("freqs_avg_part", long_postings / freqs_partitions)
        ;
}

template <typename InputCollection, typename CollectionType>
void create_collection(InputCollection& input,
                       ds2i::global_parameters const& params,
                       const char* output_filename,
                       const char* cluster_filename,
                       bool check,
                       std::string const& seq_type,
                       uint64_t universe)
{
    using namespace ds2i;

    logger() << "Processing " << input.num_docs() << " documents" << std::endl;
    double tick = get_time_usecs();
    double user_tick = get_user_time_usecs();

    typename CollectionType::builder builder(universe, params);
    progress_logger plog;

    std::ifstream file(cluster_filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    std::istream instream(&inbuf);

    std::string line;
    uint32_t c = 1;
    std::vector<pos_t> docs_positions;
    std::vector<pos_t> freqs_positions;
    while (std::getline(instream, line))
    {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "));
        uint32_t cluster_size = std::atoi(data.front().data());
        std::cout << "cluster size: " << cluster_size << "\n";

        docs_positions.reserve(cluster_size);
        freqs_positions.reserve(cluster_size);

        std::for_each(++data.begin(), data.end(),
        [&docs_positions, &freqs_positions](std::string const& s)
        {
            pos_t p = std::stoull(s.data());
            docs_positions.push_back(p);
            freqs_positions.push_back(p - 2);
        });

        input.add_positions(docs_positions, freqs_positions);
        input.set_positions();

        for (auto const& plist : input)
        {
            uint64_t freqs_sum = std::accumulate(plist.freqs.begin(),
                                                 plist.freqs.end(),
                                                 uint64_t(0));
            builder.add_posting_list(plist.docs.size(),
                                     plist.docs.begin(),
                                     plist.freqs.begin(),
                                     freqs_sum);
            plog.done_sequence(plist.docs.size());
        }

        logger() << "cluster-" << c++ << " encoded:\n";
        builder.print_cluster_stat();

        docs_positions.clear();
        freqs_positions.clear();
    }

    plog.log();
    CollectionType coll;
    builder.build(coll);
    double elapsed_secs = (get_time_usecs() - tick) / 1000000;
    double user_elapsed_secs = (get_user_time_usecs() - user_tick) / 1000000;
    logger() << seq_type << " collection built in "
             << elapsed_secs << " seconds" << std::endl;

    stats_line()
        ("type", seq_type)
        ("worker_threads", configuration::get().worker_threads)
        ("construction_time", elapsed_secs)
        ("construction_user_time", user_elapsed_secs)
        ;

    dump_stats(coll, seq_type, plog.postings);
    dump_index_specific_stats(coll, seq_type);

    if (output_filename) {
        succinct::mapper::freeze(coll, output_filename);
        if (check) {
            verify_collection<InputCollection, CollectionType>(input, output_filename);
        }
    }
}

template <typename InputCollection, typename CollectionType>
void create_collection(InputCollection const& input,
                       ds2i::global_parameters const& params,
                       const char* output_filename, bool check,
                       std::string const& seq_type,
                       uint64_t universe)
{
    using namespace ds2i;

    logger() << "Processing " << input.num_docs() << " documents" << std::endl;
    double tick = get_time_usecs();
    double user_tick = get_user_time_usecs();

    typename CollectionType::builder builder(universe, params);
    progress_logger plog;

    for (auto const& plist: input) {
        uint64_t freqs_sum = std::accumulate(plist.freqs.begin(),
                                             plist.freqs.end(), uint64_t(0));

        builder.add_posting_list(plist.docs.size(), plist.docs.begin(),
                                 plist.freqs.begin(), freqs_sum);
        plog.done_sequence(plist.docs.size());
    }

    plog.log();
    CollectionType coll;
    builder.build(coll);
    double elapsed_secs = (get_time_usecs() - tick) / 1000000;
    double user_elapsed_secs = (get_user_time_usecs() - user_tick) / 1000000;
    logger() << seq_type << " collection built in "
             << elapsed_secs << " seconds" << std::endl;

    stats_line()
        ("type", seq_type)
        ("worker_threads", configuration::get().worker_threads)
        ("construction_time", elapsed_secs)
        ("construction_user_time", user_elapsed_secs)
        ;

    dump_stats(coll, seq_type, plog.postings);
    dump_index_specific_stats(coll, seq_type);

    if (output_filename) {
        succinct::mapper::freeze(coll, output_filename);
        if (check) {
            verify_collection<InputCollection, CollectionType>(input, output_filename);
        }
    }
}

void set_check(int argc, const char** argv, int k, bool& check) {
    if (argc > k && std::string(argv[k]) == "--check") {
        check = true;
    }
}

int main(int argc, const char** argv)
{
    using namespace ds2i;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <index type> <collection basename> <universe>\
                    [--clusters] [cluster filename] [<output filename>] [--check]"
                  << std::endl;
        return 1;
    }

    std::string type = argv[1];
    const char* input_basename = argv[2];
    const uint64_t universe = std::stoull(argv[3]);

    const char* cluster_filename = nullptr;
    const char* output_filename = nullptr;
    bool check = false;

    if (argc > 4) {
        if (std::string(argv[4]) == "--clusters") {
            cluster_filename = argv[5];
            if (argc > 6 && std::string(argv[6]) == "--check") {
                check = true;
            } else {
                output_filename = argv[6];
                if (argc > 7 && std::string(argv[7]) == "--check") {
                    check = true;
                }
            }
        } else {
            output_filename = argv[5];
            set_check(argc, argv, 5, check);
        }
    }

    global_parameters params;
    params.log_partition_size = configuration::get().log_partition_size;

    binary_freq_collection input(input_basename);
    clustered_binary_freq_collection clustered_input(input_basename, check);

    if (false) {
#define LOOP_BODY(R, DATA, T)                                                                             \
        } else if (type == BOOST_PP_STRINGIZE(T)) {                                                       \
            if (cluster_filename) {                                                                       \
                create_collection<clustered_binary_freq_collection,                                       \
                                  BOOST_PP_CAT(T, _index)>                                                \
                    (clustered_input, params, output_filename, cluster_filename, check, type, universe);  \
            } else {                                                                                      \
                create_collection<binary_freq_collection,                                                 \
                                  BOOST_PP_CAT(T, _index)>                                                \
                    (input, params, output_filename, check, type, universe);                              \
            }                                                                                             \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, DS2I_INDEX_TYPES);
#undef LOOP_BODY
    } else {
        logger() << "ERROR: Unknown type " << type << std::endl;
    }

    return 0;
}
