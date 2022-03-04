#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

    int bam_sort(const char* inputName, const char* outputName, const char* tmpDir, bool useTmpDir,
                 int numThreads, int merge_threads, size_t memory, int* numFiles, int* numBlocks);

#ifdef __cplusplus
}
#endif
