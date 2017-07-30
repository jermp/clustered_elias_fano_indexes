#include <fstream>
#include <iostream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

#include "succinct/mapper.hpp"
#include "binary_freq_collection.hpp"
#include "binary_collection.hpp"
#include "wand_data.hpp"
#include "util.hpp"

void fill_clusters(ds2i::clustered_binary_freq_collection& coll,
                   const char* cluster_filename)
{
    using ds2i::logger;

    std::ifstream file(cluster_filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    std::istream instream(&inbuf);

    std::string line;
    std::vector<pos_t> docs_positions;
    std::vector<pos_t> freqs_positions;
    
    while (std::getline(instream, line))
    {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "));
        uint32_t cluster_size = std::atoi(data.front().data());
        
        logger() << "cluster_size: " << cluster_size << std::endl;
        
        docs_positions.reserve(cluster_size);
        freqs_positions.reserve(cluster_size);

        std::for_each(++data.begin(), data.end(),
        [&docs_positions, &freqs_positions](std::string const& s)
        {
            pos_t p = std::stoull(s.data());
            docs_positions.push_back(p);
            freqs_positions.push_back(p - 2);
        });

        coll.add_positions(docs_positions, freqs_positions);
        docs_positions.clear();
        freqs_positions.clear();
    }
}

int main(int argc, const char** argv)
{
    using namespace ds2i;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <collection basename> <output filename> [<cluster filename>]"
                  << std::endl;
        return 1;
    }

    const char* input_basename = argv[1];
    const char* output_filename = argv[2];
    const char* cluster_filename = nullptr;

    if (argc > 3) {
        cluster_filename = argv[3];
    }

    binary_collection sizes_coll((std::string(input_basename) + ".sizes").data());
    
    if (cluster_filename) {
        clustered_binary_freq_collection coll(input_basename, true);
        fill_clusters(coll, cluster_filename);
        wand_data<> wdata(sizes_coll.begin()->begin(), coll.num_docs(), coll);
        succinct::mapper::freeze(wdata, output_filename);
    } else {
        binary_freq_collection coll(input_basename);
        wand_data<> wdata(sizes_coll.begin()->begin(), coll.num_docs(), coll);
        succinct::mapper::freeze(wdata, output_filename);
    }

    return 0;
}
