#pragma once

#include <stdexcept>
#include <iterator>
#include <stdint.h>

#include "clustered_binary_collection.hpp"

namespace ds2i {

    class clustered_binary_freq_collection {
    public:

        clustered_binary_freq_collection(const char* basename, bool check)
            : m_docs((std::string(basename) + ".docs").c_str())
            , m_freqs((std::string(basename) + ".freqs").c_str())
            , m_clusters(0)
            , m_terms(0)
            , m_check(check)
        {
            auto firstseq = *m_docs.begin();
            if (firstseq.size() != 1) {
                throw std::invalid_argument("First sequence should only contain number of documents");
            }
            m_num_docs = *firstseq.begin();
        }

        void add_positions(std::vector<pos_t> const& docs_positions,
                           std::vector<pos_t> const& freqs_positions)
        {
            ++m_clusters;
            m_terms += docs_positions.size();

            if (m_check) {
                m_docs_positions.push_back(docs_positions);
                m_freqs_positions.push_back(freqs_positions);
            } else {
                m_docs.set_positions(docs_positions);
                m_freqs.set_positions(freqs_positions);
            }
        }

        void set_positions()
        {
            if (m_check) {
                m_docs.set_positions(m_docs_positions.back());
                m_freqs.set_positions(m_freqs_positions.back());
            }
        }

        void set_positions(uint32_t i)
        {
            m_docs.set_positions(m_docs_positions[i]);
            m_freqs.set_positions(m_freqs_positions[i]);
        }

        uint64_t num_terms() const
        {
            return m_terms;
        }

        uint32_t num_clusters() const
        {
            return m_clusters;
        }

        uint64_t num_docs() const
        {
            return m_num_docs;
        }

        class iterator;

        iterator begin() const
        {
            return iterator(m_docs.begin(), m_freqs.begin());
        }

        iterator end() const
        {
            return iterator(m_docs.end(), m_freqs.end());
        }

        struct sequence {
            clustered_binary_collection::sequence docs;
            clustered_binary_collection::sequence freqs;
        };

        class iterator : public std::iterator<std::forward_iterator_tag,
                                              sequence> {
        public:
            iterator()
            {}

            value_type const& operator*() const
            {
                return m_cur_seq;
            }

            value_type const* operator->() const
            {
                return &m_cur_seq;
            }

            iterator& operator++()
            {
                m_cur_seq.docs = *++m_docs_it;
                m_cur_seq.freqs = *++m_freqs_it;
                return *this;
            }

            bool operator==(iterator const& other) const
            {
                return m_docs_it == other.m_docs_it;
            }

            bool operator!=(iterator const& other) const
            {
                return !(*this == other);
            }

        private:
            friend class clustered_binary_freq_collection;

            iterator(clustered_binary_collection::iterator docs_it,
                     clustered_binary_collection::iterator freqs_it)
                : m_docs_it(docs_it)
                , m_freqs_it(freqs_it)
            {
                m_cur_seq.docs = *m_docs_it;
                m_cur_seq.freqs = *m_freqs_it;
            }

            clustered_binary_collection::iterator m_docs_it;
            clustered_binary_collection::iterator m_freqs_it;
            sequence m_cur_seq;
        };

    private:
        clustered_binary_collection m_docs;
        clustered_binary_collection m_freqs;
        std::vector<std::vector<pos_t>> m_docs_positions;
        std::vector<std::vector<pos_t>> m_freqs_positions;
        uint64_t m_num_docs;
        uint32_t m_clusters;
        uint64_t m_terms;
        bool m_check;
    };
}
