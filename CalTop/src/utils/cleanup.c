void cleanup_filter_files(void)
{
    const char *files[] = {
        "drow.bin",
        "dnnz.bin",
        "dcol.bin",
        "dsum.bin",
        "dval.bin"
    };

    int nfiles = sizeof(files) / sizeof(files[0]);

    for (int i = 0; i < nfiles; i++)
    {
        remove(files[i]);
    }
}