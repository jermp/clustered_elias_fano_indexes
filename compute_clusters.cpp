#include <algorithm>
#include <dirent.h>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <cmath>

#include <unordered_map>
#include <list>
#include <deque>
#include <vector>

#include "clustered_binary_freq_collection.hpp"
#include "clustered_binary_collection.hpp"
#include "partitioned_sequence.hpp"
#include "indexed_sequence.hpp"
#include "configuration.hpp"
#include "cluster.hpp"
#include "util.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/copy.hpp>

std::vector<pos_t>
read_plists_positions(const char* positions_fn,
                      uint32_t num_lists)
{
    std::ifstream file(positions_fn, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    std::istream instream(&inbuf);
    std::string line;
    std::vector<pos_t> plists_positions;
    plists_positions.reserve(num_lists);
    while (std::getline(instream, line)) {
        auto pos = stoull(line);
        if (pos) {
            plists_positions.push_back(pos);
        }
    }
    file.close();
    return plists_positions;
}

std::vector<plist_t>
preprocess_collection(const char* bin_coll_fn,
                      std::vector<pos_t>& plists_positions,
                      std::vector<double>& itfs,
                      uint64_t universe)
{
    uint32_t terms = plists_positions.size();
    std::vector<plist_t> plists;
    plists.reserve(terms);

    ds2i::clustered_binary_freq_collection in(bin_coll_fn, false);
    std::string full_name = std::string(bin_coll_fn) + ".docs";
    ds2i::clustered_binary_collection input(full_name.data());
    ds2i::global_parameters params;
    input.set_positions(plists_positions);

    int elias_fano_type =
        ds2i::indexed_sequence::index_type::elias_fano;

    uint32_t i = 0;
    for (auto const& list: input)
    {
        auto partition =
            ds2i::partitioned_sequence<>::compute_partition(list.begin(),
                                                            universe,
                                                            list.size(),
                                                            params);
        sequence_t pruned_list;
        pruned_list.reserve(list.size()); // actual size will be less than this
        for (uint32_t i = 0; i < partition.size(); ++i) {
            auto const& v = get_block_info(i, list, partition);
            uint32_t lo = v.lo, hi = v.hi, n = v.n, u = v.u;
            int block_best_type =
                ds2i::indexed_sequence::best_type(params, u, n);
            if (block_best_type == elias_fano_type) {
                for (uint32_t k = lo; k < hi; ++k) {
                    auto doc_id = list[k];
                    pruned_list.push_back(doc_id);
                    if (itfs[doc_id]) itfs[doc_id] += 1.0;
                    else itfs[doc_id] = 1.0;
                }
            }
        }

        pruned_list.shrink_to_fit();
        plists.emplace_back(plists_positions[i++],
                            std::move(pruned_list));
    }

    // sort (in ascending order) plists on length:
    // allows quicker seed_selection and defines
    // the order in which plists are added to clusters
    std::sort(plists.begin(), plists.end(),
    [&](plist_t const& x, plist_t const& y) {
        return x.second.size() < y.second.size();
    });

    double M = (double)terms;
    for (auto& x: itfs) {
        x = std::log(M / x);
    }

    ds2i::logger() << "collection preprocessed" << std::endl;

    return plists;
}

template<typename DistanceFunction,
         typename SeedSelectionFunction>
std::vector<ds2i::cluster>
kmeans(std::vector<plist_t>& plists,
       std::vector<ds2i::cluster::index_t> const& plists_indexes,
       std::vector<double> const& itfs,
       DistanceFunction d,
       SeedSelectionFunction seed_selection,
       uint32_t MAX_ITER,
       uint32_t TRIALS,
       uint32_t CUT_OFF)
{
    auto const& centroids_indexes
        = seed_selection(plists_indexes, d, plists,
                         itfs, TRIALS, CUT_OFF);

    uint32_t k = centroids_indexes.size();

    std::vector<ds2i::cluster> clusters;
    clusters.reserve(k);

    for (uint32_t i = 0; i < k; ++i) {
        clusters.emplace_back(centroids_indexes[i], plists, itfs);
    }

    bool termination = false;
    uint32_t iterations = 0;
    while (!termination)
    {
        // std::cout << "determining closest cluster" << std::endl;
        for (auto& c : clusters) c.incr_iter();
        ++iterations;

        auto const& conf = ds2i::configuration::get();
        boost::experimental::parallel::v2::task_region(*conf.executor,
        [&](ds2i::task_region_handle& trh)
        {   // assign plist to best cluster
            for (auto const& p_index : plists_indexes)
            {
                trh.run([&, p_index]
                {
                    double smallest_d
                        = std::numeric_limits<double>::max();
                    uint32_t closer_cluster_index = 0;
                    for (uint32_t i = 0; i < k; ++i)
                    {   // select best cluster
                        double cur_d = d(p_index,
                                         clusters[i],
                                         plists, itfs);
                        if (cur_d < smallest_d)
                        {
                            smallest_d = cur_d;
                            closer_cluster_index = i;
                        }
                    }
                    clusters[closer_cluster_index].add_plist_index(p_index);
                });
            }
        });

        bool empty_cluster = false;
        for (uint32_t i = 0; i < k; ++i) {
            if (clusters[i].size() < 2) {
                empty_cluster = true;
            }
        }

        if (empty_cluster) {
            break;
        }

        termination = true;
        if (iterations != MAX_ITER)
        {
            // check if we can stop:
            // sufficient to check k - 1 clusters
            for (uint32_t i = 0; i < k - 1; ++i)
            {
                auto& c = clusters[i];
                if (!c.same_as_before()) {
                    c.dump();
                    termination = false;
                }
            }
            if (!termination) {
                auto& c = clusters[k - 1];
                c.dump();
            }
        }
    }

    return clusters;
}

template<typename DistanceFunction,
         typename SeedSelectionFunction,
         typename NeedsPartitionPredicate>
std::list<ds2i::cluster>
compute_clusters(std::vector<plist_t>& plists
                , std::vector<double> const& itfs
                , uint32_t MAX_ITER
                , uint32_t MAX_REF_SIZE
                , uint32_t F
                , uint32_t TRIALS
                , uint32_t CUT_OFF
                , DistanceFunction d
                , SeedSelectionFunction seed_selection
                , NeedsPartitionPredicate needs_partition)
{
    ds2i::cluster root(plists, itfs);
    uint32_t terms = plists.size();
    for (uint32_t i = 0; i < terms; ++i)
    {   // initial cluster containing all plists' positions
        root.add_plist_index(i);
    }

    std::deque<ds2i::cluster> to_split;
    to_split.push_back(std::move(root));
    std::list<ds2i::cluster> final_clusters;

    while (to_split.size())
    {
        auto& parent = to_split.front();
        auto children = std::move(kmeans(plists
                                        , parent.plists_indexes()
                                        , itfs, d, seed_selection
                                        , MAX_ITER, TRIALS, CUT_OFF));
        to_split.pop_front(); // destroy front item
        for (auto& c : children)
        {
            if (c.size()) // cluster may be empty
            {
                bool np = needs_partition(c, itfs, MAX_REF_SIZE, F);
                c.clean();
                if (np) {
                    // ds2i::logger() << "pushing new cluster" << std::endl;
                    to_split.push_back(std::move(c));
                } else {
                    final_clusters.push_back(std::move(c));
                    // ds2i::logger() << "final cluster created" << std::endl;
                    // free plists associated to the positions
                    // in the cluster
                    auto const& plists_indexes = c.plists_indexes();
                    for (auto index: plists_indexes) {
                        std::vector<posting_t>().swap(plists[index].second);
                    }
                    // ds2i::logger() << "memory freed" << std::endl;
                }
            }
        }
    }
    return final_clusters;
}

// argv[1] --> binary collection filename
// argv[2] --> plists' positions gzipped file
// argv[3] --> universe
// argv[4] --> num of plists
// argv[5] --> threshold frequency of postings' selection
// argv[6] --> number of trials in seed_selection
// argv[7] --> cut-off threshold in seed_selection
// argv[8] --> divisor of universe to define maximum reference size
// argv[9] --> maximum number of k-means iterations
int main(int argc, char** argv)
{
    if (argc < 10) return 1;

    const char* bin_coll_fn = argv[1];
    const char* positions_fn = argv[2];
    uint64_t U = std::stoull(argv[3]);
    uint32_t num_lists = std::atoi(argv[4]);
    uint32_t F = std::atoi(argv[5]);
    uint32_t TRIALS = std::atoi(argv[6]);
    uint32_t CUT_OFF = std::atoi(argv[7]);
    uint32_t DIV = std::atoi(argv[8]);
    uint32_t MAX_ITER = std::atoi(argv[9]);

    auto plists_positions =
        read_plists_positions(positions_fn, num_lists);

    std::vector<double> itfs(U, 0.0);

    auto plists = preprocess_collection(bin_coll_fn,
                                        plists_positions,
                                        itfs, U);
    typedef std::function<double(uint32_t
                                , ds2i::cluster const&
                                , std::vector<plist_t>&
                                , std::vector<double> const&)> d_fun_t;
    // distance is defined as 1 - cosine_similarity
    d_fun_t d = [&](uint32_t plist_index
                   , ds2i::cluster const& cluster
                   , std::vector<plist_t>& plists
                   , std::vector<double> const& itfs)
    {
        double dot_product = 0.0;
        double squared_list_sum = 0.0;
        auto const& plist = plists[plist_index].second;
        auto const& centroid = cluster.centroid();
        for (auto const& i : plist) {
            double itfs_i = itfs[i];
            dot_product += centroid[i] * itfs_i;
            squared_list_sum += itfs_i * itfs_i;
        }

        double sim = dot_product /
            (cluster.centroid_norm2() *
             std::sqrt(squared_list_sum));

        return 1.0 - sim;
    };

    auto seed_selection = [&](std::vector<ds2i::cluster::index_t> const& plists_indexes
                        , d_fun_t d
                        , std::vector<plist_t>& plists
                        , std::vector<double> const& itfs
                        , uint32_t TRIALS
                        , uint32_t CUT_OFF)
    {
        // std::cout << "selecting seeds" << std::endl;
        std::vector<uint32_t> seeds; // seeds to be returned
        seeds.reserve(2);

        std::vector<ds2i::cluster::index_t> potential_seeds_indexes;
        uint32_t k = plists_indexes.size();

        if (plists_indexes.size() > CUT_OFF)
        {
            // pick as base_length the avg length
            size_t s = 0.0;
            for (auto const& i : plists_indexes) {
                s += plists[i].second.size();
            }
            double base_length = double(s) / plists_indexes.size();

            auto comp = [&](double val, uint32_t i) {
                return val < plists[i].second.size();
            };

            auto it_first = plists_indexes.begin();
            auto it_last = plists_indexes.end();
            float i = 1.0;
            uint32_t kk = 0;
            while (kk < 2)
            {
                // std::cout << "looping" << std::endl;
                double perc_less = (1.0 - i * 0.15) * base_length;
                double perc_more = (1.0 + i * 0.15) * base_length;

                it_first = std::upper_bound(plists_indexes.begin(),
                                            plists_indexes.end(),
                                            perc_less, comp);
                it_last = std::upper_bound(plists_indexes.begin(),
                                            plists_indexes.end(),
                                            perc_more, comp);
                kk = it_last == it_first ? 0 : it_last - it_first;
                k = kk;
                i += 0.05; // 5%
            }

            potential_seeds_indexes.reserve(k);
            potential_seeds_indexes.insert(potential_seeds_indexes.begin(),
                                            it_first, it_last);
        }

        auto const& candidate_seeds_indexes =
            plists_indexes.size() > CUT_OFF
            ? potential_seeds_indexes
            : plists_indexes;

        std::random_device rd;
        std::mt19937 index_dis_gen(rd());

        std::uniform_int_distribution<uint32_t> index_dis(0, k - 1);

        // First seed is drawn uniformly at random.
        uint32_t seed1_index =
            candidate_seeds_indexes[index_dis(index_dis_gen)];
        seeds.push_back(seed1_index);

        // Second seed is chosen using a random selection
        // from a non-uniform distribution as described in
        // "kmeans++: The Advantages of Careful Seeding",
        // SODA 2007, by D. Arthur and S. Vassilvitski.
        ds2i::cluster seed1_cluster(seed1_index
                                  , plists
                                  , itfs);

        std::vector<double> distribution;
        distribution.reserve(k - 1);
        double squared_sum = 0;
        for (auto const& i : candidate_seeds_indexes)
        {
            if (i != seed1_index)
            {   // generate non-uniform distribution using d^2
                double dist = d(i, seed1_cluster, plists, itfs);
                double squared_dist = dist * dist;
                squared_sum += squared_dist;
                distribution.push_back(squared_dist);
            }
        }

        std::default_random_engine prob_dis_gen;
        std::uniform_real_distribution<double> prob_dis(0.0, squared_sum);
        double best_d = 0.0;
        uint32_t best_seed_index = 0;

        for (uint32_t i = 0; i < TRIALS; ++i)
        {   // repeat body for TRIALS times:
            // mitigate randomness
            double P = prob_dis(prob_dis_gen);
            double cumulative_prob = 0.0;
            uint32_t index = 0;
            for (auto const& p : distribution)
            {
                cumulative_prob += p;
                if (cumulative_prob > P) break;
                ++index;
            }
            double cur_d = distribution[index];
            if (cur_d > best_d) // pick farthest seeds
            {
                best_d = cur_d;
                best_seed_index = index;
            }
        }
        seeds.push_back(candidate_seeds_indexes[best_seed_index]);

        return seeds;
    };

    auto needs_partition = [&](ds2i::cluster const& c
                            , std::vector<double> const& itfs
                            , uint32_t MAX_REF_SIZE
                            , uint32_t F)
    {
        // adjust F
        if (c.size() <= 3) return false;
        if (c.size() <= 10) F = 2;

        uint32_t reference_size = 0;
        auto const& centroid = c.centroid();
        uint32_t iterations = c.iterations();
        for (uint32_t i = 0; i < centroid.size(); ++i)
        {
            if (centroid[i]) {
                if (std::ceil(centroid[i] / (itfs[i] * iterations)) > F) {
                    ++reference_size;
                }
            }
        }

        // Reference size must not exceed the given threshold
        // and reference must be used for at least the given
        // percentage of its ints.
        // std::cout << "\treference_size: " << reference_size << std::endl;
        // std::cout << "\tmax_reference_size: " << MAX_REF_SIZE << std::endl;
        double redundance = double(reference_size) / c.ints();
        // std::cout << "\tredundance: " << std::setprecision(4)
        //         << redundance * 100.0 << "%" << std::endl;
        return reference_size > MAX_REF_SIZE
            || redundance > 0.3; // at most 30% of redundance
    };

    auto tick = ds2i::get_time_usecs();
    auto clusters = compute_clusters(plists, itfs
                                    , MAX_ITER, U / DIV
                                    , F, TRIALS, CUT_OFF
                                    , d, seed_selection
                                    , needs_partition); // needs_partition
    double elapsed = double(ds2i::get_time_usecs() - tick);
    ds2i::logger() << "clustering takes " << elapsed / 1000000 << " seconds" << std::endl;


    std::vector<pos_t> outliers; // group outliers in one cluster

    for (auto& c : clusters) { // print all clusters
        if (c.size() < 2) {
            auto plists_indexes = c.plists_indexes();
            std::copy(plists_indexes.begin(),
                      plists_indexes.end(),
                      std::back_inserter(outliers));
        } else {
            // print all the others
            auto plists_indexes = c.plists_indexes();
            size_t cluster_size = plists_indexes.size();
            std::cout << cluster_size << " ";
            uint32_t k = 0;
            for (auto i : plists_indexes) {
                std::cout << plists[i].first;
                if (++k != cluster_size) std::cout << " ";
            }
            std::cout << std::endl;
        }
    }

    if (outliers.size()) {
        size_t outliers_size = outliers.size();
        std::cout << outliers_size << " ";
        uint32_t k = 0;
        for (auto i : outliers) {
            std::cout << plists[i].first;
            if (++k != outliers_size) std::cout << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
