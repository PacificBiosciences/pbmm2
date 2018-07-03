#pragma once

#include <minimap.h>

namespace PacBio {
namespace minimap2 {

struct ThreadBuffer
{
    ThreadBuffer() { tbuf_ = mm_tbuf_init(); }

    ~ThreadBuffer() { mm_tbuf_destroy(tbuf_); }

    mm_tbuf_t* tbuf_;
};
}  // namespace minimap2
}  // namespace PacBio