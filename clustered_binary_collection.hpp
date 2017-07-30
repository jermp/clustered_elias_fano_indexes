#pragma once

#include <boost/iostreams/device/mapped_file.hpp>
#include <stdexcept>
#include <iterator>
#include <stdint.h>
#include <sys/mman.h>

#include "util.hpp"

namespace ds2i
{
    class clustered_binary_collection
    {
        public:
        clustered_binary_collection(const char* filename)
            : m_plists_positions(0)
            , m_terms(1)
        {
            m_file.open(filename);
            if (!m_file.is_open()) throw std::runtime_error("Error opening file");
            m_data = (sequence_iterator_t)m_file.data();
            auto ret = posix_madvise((void*)m_data, m_file.size() / sizeof(m_data[0]), POSIX_MADV_SEQUENTIAL);
            if (ret) logger() << "Error calling madvice: " << errno << std::endl;
        }

        class iterator;

        void set_positions(std::vector<pos_t> const& plists_positions)
        {
            m_plists_positions = &(plists_positions.data()[0]);
            m_terms = plists_positions.size();
        }

        iterator begin() const
        {
            return iterator(this, m_plists_positions, 0);
        }

        iterator end() const
        {
            return iterator(this, m_plists_positions, m_terms);
        }

        class sequence
        {
            public:
            sequence()
                : m_begin(nullptr)
                , m_end(nullptr)
            {}

            sequence_iterator_t begin() const
            {
                return m_begin;
            }

            sequence_iterator_t end() const
            {
                return m_end;
            }

            posting_t back() const
            {
                assert(size());
                return *(m_end - 1);
            }

            posting_t operator[](uint32_t i) const
            {
                assert(i < size());
                return *(m_begin + i);
            }

            size_t size() const
            {
                return m_end - m_begin;
            }

            private:
            friend class clustered_binary_collection::iterator;

            sequence(sequence_iterator_t begin,
                     sequence_iterator_t end)
                : m_begin(begin)
                , m_end(end)
            {}

            sequence_iterator_t m_begin;
            sequence_iterator_t m_end;
        };

        class iterator : public std::iterator<std::forward_iterator_tag,
                                               sequence>
        {
            public:
            iterator()
                : m_collection(nullptr)
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
                m_pos++;
                read();
                return *this;
            }

            bool operator==(iterator const& other) const
            {
                assert(m_collection == other.m_collection);
                return m_pos == other.m_pos;
            }

            bool operator!=(iterator const& other) const
            {
                return !(*this == other);
            }

            private:
            friend class clustered_binary_collection;

            iterator(clustered_binary_collection const* coll,
                     pos_t const* plists_positions,
                     size_t pos)
                : m_collection(coll)
                , m_pos(pos)
                , m_plists_positions(plists_positions)
            {
                read();
            }

            void read()
            {
                if (m_pos == m_collection->m_terms) return;
                
                size_t n = 0;
                pos_t pos = !m_plists_positions ? // if no position is given, read first singleton plist
                            0 :
                            m_plists_positions[m_pos];
                while (!(n = m_collection->m_data[pos++])); // skip empty seqs
                sequence_iterator_t begin = &m_collection->m_data[pos];
                sequence_iterator_t end = begin + n;
                m_cur_seq = sequence(begin, end);
            }

            clustered_binary_collection const* m_collection;
            size_t m_pos;
            sequence m_cur_seq;
            pos_t const* m_plists_positions;
        };

        private:
        boost::iostreams::mapped_file_source m_file;
        sequence_iterator_t m_data;
        pos_t const* m_plists_positions;
        uint32_t m_terms;
    };
}
