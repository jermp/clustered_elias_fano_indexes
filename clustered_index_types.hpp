#pragma once

#include "clustered_freq_index.hpp"
#include "positive_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "strict_sequence.hpp"

namespace ds2i {
    typedef clustered_freq_index<
                positive_sequence<partitioned_sequence<strict_sequence>>
            > clustered_opt_index;
}
