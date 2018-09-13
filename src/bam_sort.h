// Author: Armin Töpfer

#pragma once

extern "C"
{
    int bam_sort(const char *inputName, const char *outputName, int numThreads, int memory,
                 int *numFiles, int *numBlocks);
}
