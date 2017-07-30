#pragma once

#include "global_parameters.hpp"
#include "indexed_sequence.hpp"
#include "integer_codes.hpp"
#include "util.hpp"
#include "optimal_partition.hpp"
#include "partitioned_sequence.hpp"
#include "clustered_binary_collection.hpp"

namespace ds2i
{
    template <typename BaseSequence = partitioned_sequence<>>
    struct clustered_sequence
    {
        typedef BaseSequence base_sequence_type;
        typedef typename base_sequence_type::enumerator base_sequence_enumerator;

        template <typename ReferenceIterator>
        static sequence_t
        write(succinct::bit_vector_builder& bvb,
              uint64_t reference_id,
              ReferenceIterator reference_begin,
              uint64_t reference_size,
              ds2i::clustered_binary_collection::sequence const& docs,
              uint64_t universe,
              sequence_t const& partition,
              global_parameters const& params)
        {
            auto const& blocks
                = compute_partition(reference_begin,
                                    reference_begin + reference_size,
                                    docs, partition, params);

            write_gamma(bvb, reference_id);

            auto const& mapped_intersection = blocks.mapped_intersection;
            size_t mi_size = mapped_intersection.size();
            uint64_t log_mi_size = ds2i::ceil_log2(mi_size);

            write_gamma(bvb, log_mi_size);
            bvb.append_bits(mi_size, log_mi_size + 1);

            sequence_t freqs_positions;
            freqs_positions.reserve(docs.size());
            auto const& intersection = blocks.intersection;
            compute_mapping(docs.begin(),
                            intersection.begin(),
                            intersection.end(),
                            std::back_inserter(freqs_positions));
            auto const& residual = blocks.residual;
            compute_mapping(docs.begin(),
                            residual.begin(),
                            residual.end(),
                            std::back_inserter(freqs_positions));

            if (!mi_size) {
                base_sequence_type::write(bvb, docs.begin(),
                                          universe,
                                          docs.size(), params);
            } else {
                succinct::bit_vector_builder bvb_mi;
                base_sequence_type::write(bvb_mi,
                                          mapped_intersection.begin(),
                                          reference_size,
                                          mapped_intersection.size(),
                                          params);
                bvb.append_bits(bvb_mi.size(), 32);
                bvb.append(bvb_mi);
                if (residual.size()) {
                    base_sequence_type::write(bvb, residual.begin(),
                                              universe, residual.size(),
                                              params);
                }
            }

            return freqs_positions;
        }

        template <typename ReferenceIterator>
        static void write_reference(succinct::bit_vector_builder& bvb,
                                    ReferenceIterator begin,
                                    uint64_t universe,
                                    size_t n,
                                    global_parameters const& params)
        {
            base_sequence_type::write(bvb, begin, universe, n, params);
        }

        class enumerator
        {
            public:

            typedef std::pair<uint64_t, uint64_t> value_type; // (position, value)

            enumerator()
            {}

            enumerator(succinct::bit_vector const& bv_docs,
                       uint64_t offset_docs,
                       succinct::bit_vector const& bv_refs,
                       uint64_t offset_ref,
                       uint64_t universe, uint64_t n,
                       global_parameters const& params)
                : m_last_returned_val(0)
                , m_last_res_val(0)
                , m_last_mi_val(0)
                , m_last_returned_pos(0)
                , m_last_res_pos(0)
                , m_last_mi_pos(0)
                , m_size(n)
                , m_mi_size(0)
                , m_res_size(0)
                , m_params(params)
            {
                succinct::bit_vector::enumerator docs_it(bv_docs, offset_docs);
                uint64_t log_mi_size = read_gamma(docs_it);
                uint64_t residual_offset = 0;
                m_mi_size = docs_it.take(log_mi_size + 1);

                if (DS2I_LIKELY(m_mi_size))
                {
                    succinct::bit_vector::enumerator ref_it(bv_refs, offset_ref);
                    uint64_t log_reference_size = read_gamma_nonzero(ref_it);
                    uint64_t reference_size = 1;

                    if (log_reference_size > 1) {
                        reference_size = ref_it.take(log_reference_size + 1);
                    }

                    m_ref_enum = base_sequence_enumerator(bv_refs, ref_it.position(),
                                                          universe, reference_size,
                                                          m_params);

                    residual_offset = docs_it.take(32);
                    m_mi_enum = base_sequence_enumerator(bv_docs, docs_it.position(),
                                                         reference_size, m_mi_size,
                                                         m_params);
                }

                m_res_size = n - m_mi_size;
                if (DS2I_LIKELY(m_res_size)) {
                    m_res_enum = base_sequence_enumerator(bv_docs,
                                                          docs_it.position() + residual_offset,
                                                          universe, m_res_size, m_params);
                }
            }

            value_type DS2I_ALWAYSINLINE next_geq(uint64_t lower_bound)
            {
                if (DS2I_LIKELY(m_res_size && m_mi_size)) {

                    if (DS2I_LIKELY(lower_bound > m_last_mi_val)) {
                        mapped_intersection_next_geq(lower_bound);
                        if (DS2I_LIKELY(lower_bound > m_last_res_val)) {
                            residual_next_geq(lower_bound);
                        }
                        return value();
                    } else if (DS2I_LIKELY(lower_bound > m_last_res_val)) {
                        residual_next_geq(lower_bound);
                        return value();
                    } else {
                        return value_type(m_last_returned_pos,
                                          m_last_returned_val);
                    }
                }
                else if (DS2I_UNLIKELY(m_res_size == 0)) {

                    auto ref_pv = m_ref_enum.next_geq(lower_bound);
                    auto mi_pv = m_mi_enum.next_geq(ref_pv.first);
                    if (mi_pv.second == ref_pv.first) {
                        return value_type(mi_pv.first, ref_pv.second);
                    } else {
                        return value_type(mi_pv.first,
                            m_ref_enum.move(mi_pv.second).second);
                    }

                } else {
                    auto val = m_res_enum.next_geq(lower_bound);
                    return value_type(val.first, val.second);
                }
            }

            value_type DS2I_ALWAYSINLINE next()
            {
                if (DS2I_LIKELY(m_res_size && m_mi_size))
                {
                    if (DS2I_LIKELY(m_last_returned_val == m_last_res_val)) {
                        residual_next();
                    } else {
                        mapped_intersection_next();
                    }
                    return value();
                } else if (DS2I_UNLIKELY(m_res_size == 0)) {
                    auto val = m_mi_enum.next();
                    return value_type(val.first,
                        m_ref_enum.move(val.second).second);
                } else {
                    auto val = m_res_enum.next();
                    return value_type(val.first, val.second);
                }
            }

            value_type DS2I_ALWAYSINLINE init()
            {
                if (DS2I_LIKELY(m_res_size && m_mi_size)) {
                    auto mi_pv = m_mi_enum.move(0);
                    m_last_res_val = m_res_enum.move(0).second;
                    m_last_mi_val = m_ref_enum.move(mi_pv.second).second;
                    return value();
                } else if (DS2I_UNLIKELY(m_res_size == 0)) {
                    auto val = m_mi_enum.move(0);
                    return value_type(0, m_ref_enum.move(val.second).second);
                } else {
                    return value_type(0, m_res_enum.move(0).second);
                }
            }

            uint64_t size() const {
                return m_size;
            }

            private:

            value_type DS2I_ALWAYSINLINE value()
            {
                if (m_last_res_val < m_last_mi_val) {
                    m_last_returned_val = m_last_res_val;
                    m_last_returned_pos = m_mi_size + m_last_res_pos;
                } else {
                    m_last_returned_val = m_last_mi_val;
                    m_last_returned_pos = m_last_mi_pos;
                }
                return value_type(m_last_returned_pos,
                                  m_last_returned_val);
            }

            void DS2I_ALWAYSINLINE mapped_intersection_next_geq(uint64_t lower_bound)
            {
                auto ref_pv = m_ref_enum.next_geq(lower_bound);
                auto mi_pv = m_mi_enum.next_geq(ref_pv.first);
                m_last_mi_pos = mi_pv.first;
                if (mi_pv.second == ref_pv.first) {
                    m_last_mi_val = ref_pv.second;
                } else {
                    m_last_mi_val = m_ref_enum.move(mi_pv.second).second;
                }
            }

            void DS2I_ALWAYSINLINE residual_next_geq(uint64_t lower_bound)
            {
                auto val = m_res_enum.next_geq(lower_bound);
                m_last_res_pos = val.first;
                m_last_res_val = val.second;
            }

            void DS2I_ALWAYSINLINE mapped_intersection_next()
            {
                ++m_last_mi_pos;
                m_last_mi_val = m_ref_enum.move(m_mi_enum.next().second).second;
            }

            void DS2I_ALWAYSINLINE residual_next()
            {
                ++m_last_res_pos;
                m_last_res_val = m_res_enum.next().second;
            }

            uint64_t m_last_returned_val;
            uint64_t m_last_res_val;
            uint64_t m_last_mi_val;

            uint64_t m_last_returned_pos;
            uint64_t m_last_res_pos;
            uint64_t m_last_mi_pos;

            uint64_t m_size;
            uint64_t m_mi_size;
            uint64_t m_res_size;

            global_parameters m_params;
            base_sequence_enumerator m_ref_enum;
            base_sequence_enumerator m_mi_enum;
            base_sequence_enumerator m_res_enum;
        };

        struct blocks {
            sequence_t intersection;
            sequence_t mapped_intersection;
            sequence_t residual;
        };

        template <typename ReferenceIterator>
        static blocks
        compute_partition(ReferenceIterator reference_begin,
                          ReferenceIterator reference_end,
                          ds2i::clustered_binary_collection::sequence const& docs,
                          sequence_t const& partition,
                          global_parameters const& params)
        {
            sequence_t to_be_intersected;
            to_be_intersected.reserve(docs.size());

            int elias_fano_type =
                ds2i::indexed_sequence::index_type::elias_fano;

            for (uint32_t i = 0; i < partition.size(); ++i)
            {
                auto const& v = ds2i::get_block_info(i, docs, partition);
                uint32_t lo = v.lo, hi = v.hi, n = v.n, u = v.u;
                int block_best_type
                    = ds2i::indexed_sequence::best_type(params, u, n);
                if (block_best_type == elias_fano_type) {
                    for (uint32_t k = lo; k < hi; ++k) {
                        to_be_intersected.push_back(docs[k]);
                    }
                }
            }

            sequence_t intersection;
            intersection.reserve(to_be_intersected.size());
            std::set_intersection(reference_begin,
                                  reference_end,
                                  to_be_intersected.begin(),
                                  to_be_intersected.end(),
                                  std::back_inserter(intersection));

            sequence_t mapped_intersection;
            mapped_intersection.reserve(intersection.size());
            compute_mapping(reference_begin,
                            intersection.begin(),
                            intersection.end(),
                            std::back_inserter(mapped_intersection));

            sequence_t residual;
            residual.reserve(docs.size() - intersection.size());
            std::set_difference(docs.begin(),
                                docs.end(),
                                intersection.begin(),
                                intersection.end(),
                                std::back_inserter(residual));

            return  {
                std::move(intersection),
                std::move(mapped_intersection),
                std::move(residual)
            };
        }

        template <typename ReferenceIterator,
                  typename ListIterator,
                  typename OutputIterator>
        static OutputIterator
        compute_mapping(ReferenceIterator reference_begin,
                        ListIterator list_begin,
                        ListIterator list_end,
                        OutputIterator result)
        {
            uint32_t pos = 0;
            while (list_begin != list_end) {
                auto v = *list_begin++;
                while (*reference_begin != v) {
                    pos++;
                    reference_begin++;
                }
                *result++ = pos;
            }
            return result;
        }
    };
}
