#pragma once

#include "clustered_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "indexed_sequence.hpp"
#include "global_parameters.hpp"
#include "configuration.hpp"
#include "util.hpp"

#include <vector>
#include <unordered_map>
#include <iterator>
#include <functional>
#include <stdexcept>
#include <algorithm>

namespace ds2i
{
	struct reference_selector
	{
		typedef uint32_t block_id_t;
		typedef uint16_t plist_id_t;
		typedef uint32_t posting_t;
		typedef uint32_t endpoint_t;
	    typedef uint64_t cost_t;
		typedef std::pair<posting_t, int64_t> delta_t;
		typedef std::function<cost_t(uint64_t, uint64_t)> cost_fun_t;

		typedef partitioned_sequence<> ps;

		reference_selector(clustered_binary_freq_collection const& collection,
						   std::vector<std::vector<endpoint_t>> const& partitions,	
						   global_parameters const& params,
						   uint32_t cluster_size,
						   uint32_t max_ref_size)
			: m_collection(collection)
			, m_partitions(partitions)
			, m_universe(collection.num_docs())
			, m_K(0)
			, m_F(0)
			, m_cluster_size(cluster_size)
			, m_cur_global_cost(0)
			, m_initial_global_cost(0)
			, m_total_ints(0)
			, m_max_ref_size(max_ref_size)
			, m_params(params)
			, m_plists_data(cluster_size, std::move(plist_data()))
		{
			if (!cluster_size) {
				throw std::invalid_argument("Error: cluster must not be empty");
			}

		    if (cluster_size <= 3) m_F = 1;
		    else if (cluster_size <= 10) m_F = 2;
		    else if (cluster_size > 500) m_F = 5;
		    else m_F = 3;
			
			uint32_t p = 0;
		    for (auto const& plist : m_collection)
		    {	// compute initial cost
		    	auto const& docs = plist.docs;
		    	m_total_ints += docs.size();
		    	m_initial_global_cost +=
	    			opt_cost(docs, m_partitions[p++]);
		    }
		    m_cur_global_cost = m_initial_global_cost;
		}

		double gain() const
		{
			return (int64_t(m_initial_global_cost)
				  - int64_t(m_cur_global_cost))
			      * 100.0 / int64_t(m_initial_global_cost);
		}

		std::vector<posting_t> select_reference(uint32_t divisor)
		{
			select_candidates();
			m_K = m_candidates.size() / divisor;
			m_cur_reference.reserve(m_candidates.size()); // at most
			compute_deltas();

            cost_t best_cost = cost_t(-1);
            uint32_t epochs = 1;
            std::vector<posting_t> best_reference;
            best_reference.reserve(m_candidates.size()); // at most
            uint32_t consecutive_degrading_epochs = 0;

            auto tick = get_time_usecs();
			while (m_candidates.size())
			{
				for (size_t i = 0; i < m_K && i < m_candidates.size(); ++i)
				{
					auto const& d = m_cur_deltas[i];
					if (d.second < 0) break;
					m_cur_reference.push_back(d.first);
					m_candidates.erase(d.first);

					if (m_cur_reference.size() >= m_max_ref_size) break;
				}

				++epochs;
				encode();

				if (m_cur_global_cost < best_cost)
				{
					best_cost = m_cur_global_cost;
					best_reference = m_cur_reference;
					consecutive_degrading_epochs = 0;
				}
				else ++consecutive_degrading_epochs;

				if (m_cur_reference.size() >= m_max_ref_size ||
					consecutive_degrading_epochs == 5) break;

				compute_deltas();
			}

			double elapsed = double(get_time_usecs() - tick);
			m_cur_global_cost = best_cost; // for gain()

			std::cout << "{"
					  << "\"candidates\": " << m_candidates.size()
					  << ", \"epoch_size\": " << m_K
					  << ", \"num_epochs\": " << epochs
					  << ", \"gain\": " << gain()
					  << ", \"elapsed[sec]\": " << elapsed / 1000000 // seconds
					  << "}" << std::endl;

			return best_reference;
		}

		private:
		struct plist_data {
			std::vector<posting_t> mi; // mapped intersection
			std::vector<endpoint_t> mi_p; // mapped intersection's partition
		};

		void select_candidates()
		{
			std::unordered_map<posting_t, uint32_t> freqs;
		    auto freqs_cend = freqs.cend();
		    uint32_t p = 0;
		    for (auto const& plist : m_collection)
		    {
		    	auto const& docs = plist.docs;
		    	auto const& partition = m_partitions[p++];
		        for (uint32_t i = 0; i < partition.size(); ++i)
				{
					auto const& v =
						ds2i::get_block_info(i, docs, partition);
		            uint32_t lo = v.lo, hi = v.hi, n = v.n, u = v.u;
		            int block_best_type =
		                ds2i::indexed_sequence::best_type(m_params, u, n);
		            if (block_best_type == m_elias_fano_type)
		            {
		                for (uint32_t k = lo; k < hi; ++k)
		                {
		                    auto doc_id = docs[k];
		                    if (freqs.find(doc_id) != freqs_cend) {
		                    	freqs[doc_id] += 1;
		                    } else {
		                    	freqs.emplace(doc_id, 1);
		                    }
		                }
		            }
		        }
		    }

		    auto candidates_cend = m_candidates.cend();
		    plist_id_t plist_id = 0;
		    p = 0;
		    for (auto const& plist : m_collection)
		    {
		    	auto const& docs = plist.docs;
		    	auto const& partition = m_partitions[p++];
				for (uint32_t i = 0; i < partition.size(); ++i)
				{
					auto const& v =
						ds2i::get_block_info(i, docs, partition);
		            uint32_t lo = v.lo, hi = v.hi, n = v.n, u = v.u;
		            int block_best_type =
		                ds2i::indexed_sequence::best_type(m_params, u, n);
		            if (block_best_type == m_elias_fano_type)
		            {
		                for (uint32_t k = lo; k < hi; ++k)
		                {
		                    auto candidate = docs[k];
		                    if (freqs[candidate] > m_F)
		                    {
		                    	if (m_candidates.find(candidate) != candidates_cend) {
		                    		m_candidates[candidate].emplace_back(plist_id, i);
		                    	} else {
		                    		std::vector<std::pair<plist_id_t, block_id_t>> blocks_ids;
		                    		blocks_ids.emplace_back(plist_id, i);
		                    		m_candidates.emplace(candidate, std::move(blocks_ids));
		                    	}
		                    }
		                }
		            }
		        }
		        ++plist_id;
		    }

		    for (auto const& x : freqs) {
		    	double px = double(x.second) / m_cluster_size;
		    	m_binary_entropies[x.first] = px * std::log2(1.0 / px);
		    }
		}

		void encode() {
			std::sort(m_cur_reference.begin(), m_cur_reference.end());

			m_cur_reference_partition =
				std::move(ps::compute_partition(m_cur_reference.begin(),
												m_universe,
												m_cur_reference.size(),
			    								m_params));
			m_cur_global_cost = opt_cost(m_cur_reference,
										 m_cur_reference_partition);

			uint32_t p = 0;
        	for (auto it = m_collection.begin();
        		it != m_collection.end(); ++it, ++p)
			{
				auto const& docs = it->docs;
            	task_region(*m_conf.executor, [&](task_region_handle& trh)
            	{
	                trh.run([&, p]
	                {
			            auto const& blocks
			                = clustered_sequence<>::compute_partition(m_cur_reference.begin(),
			                                    					  m_cur_reference.end(),
			                                    					  docs,
			                                    					  m_partitions[p],
			                                    					  m_params);
                		auto& mapped_intersection =
                			m_plists_data[p].mi;
                		mapped_intersection = std::move(blocks.mapped_intersection);

				        auto& mapped_intersection_partition =
				        	m_plists_data[p].mi_p;
						set_if_not_empty(mapped_intersection, mapped_intersection_partition);

		    			m_cur_global_cost +=
		    				opt_cost(mapped_intersection, mapped_intersection_partition);


	                    auto const& residual = blocks.residual;
                		if (residual.size())
                		{
			    			m_cur_global_cost +=
			    				opt_cost(residual,
			    						ps::compute_partition(residual.begin(),
			    											  m_universe,
	    									 				  residual.size(),
	    									 				  m_params));
                		}
					});
				});
			}
		}

		void set_if_not_empty(std::vector<posting_t> const& s,
							  std::vector<endpoint_t>& s_partition)
		{
	        if (!s.size()) {
            	s_partition.clear();
            } else {
            	s_partition = std::move(
    				ps::compute_partition(s.begin(), m_universe,
    									  s.size(), m_params));
    		}
		}

		template<typename RandomAccessIterator>
		cost_t opt_cost(RandomAccessIterator const& s,
						std::vector<endpoint_t> const& s_partition)
		{
			cost_t opt_cost_s = 0;
			for (uint32_t i = 0; i < s_partition.size(); ++i)
			{
				auto const& v = ds2i::get_block_info(i, s, s_partition);
            	uint32_t n = v.n, u = v.u;
				opt_cost_s += m_cost_fun(u, n);
			}
			return opt_cost_s;
		}

		void compute_deltas()
		{
			m_cur_deltas.clear();
			size_t candidates = m_candidates.size();
			m_cur_deltas.reserve(candidates);

			auto threads = std::thread::hardware_concurrency();
			size_t chunk = candidates / threads;
			std::vector<std::vector<delta_t>> regions;
			regions.reserve(threads);
			for (uint32_t i = 0; i < threads; ++i) {
				regions.emplace_back(std::vector<delta_t>());
			}

			task_region(*m_conf.executor,
				[&](task_region_handle& trh) {
					for (uint32_t i = 0; i < threads; ++i)
					{
		                trh.run([&, i] {
		                	auto& r = regions[i];
		                	auto it_begin = std::next(m_candidates.begin(), i * chunk);
		                	auto it_end = (i != threads - 1
										? std::next(m_candidates.begin(), (i + 1) * chunk)
										: m_candidates.end());
		                	r.reserve(std::distance(it_begin, it_end));
		                	for (auto it = it_begin; it != it_end; ++it) {
								auto c1 = it->first;
		                		r.emplace_back(c1, delta_x(c1));
							}
			            });
		            }
            	});

			for (uint32_t i = 0; i < threads; ++i)
			{
				std::copy(regions[i].begin(),
						  regions[i].end(),
						  std::back_inserter(m_cur_deltas));
			}

			std::sort(m_cur_deltas.begin(),
					  m_cur_deltas.end(),
					  [&](delta_t const& x, delta_t const& y)
			{
				// sort on delta cost
    			return x.second * m_binary_entropies[x.first]
    				 > y.second * m_binary_entropies[y.first];
			});
		}

		static inline
		uint32_t position(posting_t x,
						  std::vector<posting_t> const& s)
		{
			auto s_begin = s.begin();
			return std::lower_bound(s_begin, s.end(), x) - s_begin;
		}

		static inline
		uint32_t block(uint32_t x,
					   std::vector<endpoint_t> const& s)
		{
			auto s_begin = s.begin();
			uint32_t b = std::upper_bound(s_begin, s.end(), x) - s_begin;
			return b < s.size() ? b : s.size() - 1;
		}

		/* cost contribution to s of adding x to s */
		/* pos_x is an approximation because it does not take into account
		for all other postings added within an epoch. We can count them
		if we want a more precise estimation */
		delta_t delta_x_s(posting_t x,
						  std::vector<posting_t> const& s,
						  std::vector<endpoint_t> const& s_partition)
		{
			if (!s.size()) // s and s_partitions may be empty (e.g., at the beginning)
			{
				return delta_t(m_min_cost, 0);
			}

			auto pos_x = position(x, s);
			auto block_x = block(pos_x, s_partition);
			auto const& v = ds2i::get_block_info(block_x, s, s_partition);
            uint32_t n = v.n, u = v.u;
			/* also here we have an approximation: the added posting could be
			an upperbound */
			int64_t c = int64_t(m_cost_fun(u, n + 1))
            		  - int64_t(m_cost_fun(u, n));

            return delta_t(pos_x, c);
		}

		/* cost contribution to whole cluster of adding x to reference. */
		int64_t delta_x(posting_t x)
		{
			auto const& pos_cost = delta_x_s(x, m_cur_reference,
									  		 m_cur_reference_partition);
			auto m_x = pos_cost.first; // mapping of x to reference
			auto delta_x_r = pos_cost.second;

			int64_t delta_cur_residuals_sum = 0
				 , delta_cur_mapped_segments_sum = 0;
			
			auto const& blocks_x = m_candidates[x];
			auto const it = m_collection.begin();

			for (auto const& p_id_b_id : blocks_x) {

				auto plist_id = p_id_b_id.first;
				auto block_id = p_id_b_id.second;

				auto it_plist = std::next(it, plist_id);
				auto const& docs = it_plist->docs;

				auto const& partition = m_partitions[plist_id];
				auto const& v = ds2i::get_block_info(block_id, docs,
											   partition);
            	uint32_t n = v.n, u = v.u;
				/* again here, another approximation: the removed posting could be
				an upperbound */
				delta_cur_residuals_sum += int64_t(m_cost_fun(u, n))
										 - int64_t(m_cost_fun(u, n - 1));

				auto const& pd = m_plists_data[plist_id];
				delta_cur_mapped_segments_sum +=
					/* approximation: instead of 0 we should count for how many remapped
					integers are present during an epoch in mapped_segment_buffer (to be
					maintained) */
					delta_x_s(m_x, pd.mi, pd.mi_p).second;
			}
			return delta_cur_residuals_sum
					- delta_cur_mapped_segments_sum
					- delta_x_r;
		}

		clustered_binary_freq_collection const& m_collection;
		std::vector<std::vector<endpoint_t>> const& m_partitions;
		uint64_t m_universe;
		size_t m_K;
		uint32_t m_F;
		size_t m_cluster_size;

		cost_t m_cur_global_cost;
		cost_t m_initial_global_cost;
		uint32_t m_total_ints;
		uint32_t m_max_ref_size;
		
		global_parameters m_params;
		configuration const& m_conf = configuration::get();
		const int m_elias_fano_type = indexed_sequence::index_type::elias_fano;
		
		std::unordered_map<
					posting_t,
					std::vector<std::pair<plist_id_t, block_id_t>>
				> m_candidates;
		
		std::unordered_map<posting_t, double> m_binary_entropies;

		std::vector<plist_data> m_plists_data;
		std::vector<posting_t> m_cur_reference;
		std::vector<endpoint_t> m_cur_reference_partition;
		std::vector<delta_t> m_cur_deltas;

		// does not take into account fix cost
		const cost_t m_min_cost = indexed_sequence::bitsize(m_params, 1, 1);

        cost_fun_t m_cost_fun = [&](uint64_t universe, uint64_t n) {
            return indexed_sequence::bitsize(m_params, universe, n) + m_conf.fix_cost;
        };
	};
}