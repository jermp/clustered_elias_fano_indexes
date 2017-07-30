#pragma once

#include <algorithm>
#include <cmath>
#include <mutex>

namespace ds2i
{
	struct cluster {
		typedef uint16_t index_t;

		cluster(std::vector<plist_t> const& plists,
				std::vector<double> const& itfs)
			: m_total_ints(0)
			, m_centroid_norm2(0.0)
			, m_iterations(0)
			, m_plists(plists)
			, m_itfs(itfs)
		{}

		cluster(const uint32_t centroid_index,
				std::vector<plist_t> const& plists,
				std::vector<double> const& itfs)
			: m_total_ints(0)
			, m_centroid_norm2(0.0)
			, m_iterations(0)
			, m_plists(plists)
			, m_itfs(itfs)
			, m_centroid(itfs.size(), 0.0)
		{
			auto const& c = plists[centroid_index].second;
			m_total_ints = c.size();
			for (auto const& i : c) {
				double itfs_i = itfs[i];
				m_centroid[i] = itfs_i;
				m_centroid_norm2 += itfs_i * itfs_i;
			}
			m_centroid_norm2 = std::sqrt(m_centroid_norm2);
		}

		cluster(cluster&& rhs)
			: m_plists(rhs.m_plists)
			, m_itfs(rhs.m_itfs)
		{
			m_total_ints = rhs.m_total_ints;
			m_iterations = rhs.m_iterations;
			m_centroid = std::move(rhs.m_centroid);
			m_cur_plists_indexes
				= std::move(rhs.m_cur_plists_indexes);
			m_prv_plists_indexes
				= std::move(rhs.m_prv_plists_indexes);
			rhs.m_total_ints = 0;
			rhs.m_centroid_norm2 = 0.0;
			rhs.m_iterations = 0;
		}

		inline
		void incr_iter() {
			++m_iterations;
		}

		inline
		uint32_t ints() const {
			return m_total_ints;
		}

		inline
		double centroid_norm2() const {
			return m_centroid_norm2;
		}

		inline
		size_t iterations() const {
			return m_iterations;
		}

		size_t size() const {
			return m_cur_plists_indexes.size();
		}

		std::vector<double> const&
		centroid() const {
			return m_centroid;
		}

		std::vector<index_t>&
		plists_indexes() {
			return m_cur_plists_indexes;
		}

		bool same_as_before() {

			std::sort(m_prv_plists_indexes.begin(),
					  m_prv_plists_indexes.end());
			std::sort(m_cur_plists_indexes.begin(),
					  m_cur_plists_indexes.end());
			return !m_prv_plists_indexes.size()
					? false
					: m_cur_plists_indexes.begin() ==
					  m_cur_plists_indexes.end() ? false :
				  	  m_prv_plists_indexes == m_cur_plists_indexes; 
		}

		void add_plist_index(index_t plist_index)
		{
			mutex.lock();
			m_cur_plists_indexes.push_back(plist_index);
			auto const& plist = m_plists[plist_index].second;
			m_total_ints += plist.size();
			
			if (m_centroid.size()) // incremental updating of centroids
			{
				for (auto const& i : plist) {
					double& c_i = m_centroid[i];
					m_centroid_norm2 -= std::sqrt(c_i * c_i);
					c_i += m_itfs[i];
					m_centroid_norm2 += std::sqrt(c_i * c_i);
				}
			}
			mutex.unlock();
		}

		void dump() {
			m_prv_plists_indexes.swap(m_cur_plists_indexes);
			m_cur_plists_indexes.clear();
		}

		void clean() {
			// these are useless after cluster is formed
			// threfore we can free them
			std::vector<double>().swap(m_centroid);
			std::vector<index_t>().swap(m_prv_plists_indexes);
			// shrink to fit
			m_cur_plists_indexes.shrink_to_fit();
		}

	private:
		uint32_t m_total_ints;
		double m_centroid_norm2;
		uint32_t m_iterations;
		std::vector<plist_t> const& m_plists;
		std::vector<double> const& m_itfs;
		std::vector<double> m_centroid;

		std::vector<index_t> m_cur_plists_indexes;
		std::vector<index_t> m_prv_plists_indexes;
		
		std::mutex mutex; // to sych access to m_cur_plists_indexes
	};
}
