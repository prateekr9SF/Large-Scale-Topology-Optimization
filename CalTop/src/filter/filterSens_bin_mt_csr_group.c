/* CalTop: grouped CSR transpose sensitivity filtering.
 * Computes H^T D^-1 g_k = H D^-1 g_k for symmetric H, k=0..nsens-1.
 * CSR binaries: row_ptr.bin (int64, ne+1), col_idx.bin (int64, nnz),
 * dval.bin (double, nnz), dsum.bin (double, ne); all indices zero-based.
 * Independent pthread row ranges; reads CSR values from disk in chunks.
 */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>

#ifndef CSR_GROUP_CHUNK
#define CSR_GROUP_CHUNK 65536U
#endif

typedef struct {
    int id, begin, end, ne, nsens, failed;
    const int64_t *rowptr;
    double * const *scaled;
    double * const *out;
} CSRGroupArgs;

static void cg_fail(const char *s) { fprintf(stderr,"CSR grouped sensitivity: %s\n",s); exit(EXIT_FAILURE); }
static int64_t cg_file_items(const char *name, size_t bytes) {
    struct stat st;
    if (stat(name,&st) != 0) { perror(name); exit(EXIT_FAILURE); }
    if (st.st_size < 0 || (uint64_t)st.st_size % bytes) cg_fail("invalid file length");
    return (int64_t)((uint64_t)st.st_size / bytes);
}
static void cg_read(const char *name, void *p, size_t size, size_t count) {
    FILE *f=fopen(name,"rb");
    if (!f) { perror(name); exit(EXIT_FAILURE); }
    if (fread(p,size,count,f)!=count) cg_fail("short metadata read");
    if (fclose(f)) cg_fail("metadata close error");
}
static void *cg_worker(void *opaque) {
    CSRGroupArgs *a=(CSRGroupArgs*)opaque;
    FILE *fc=NULL,*fw=NULL;
    int64_t *cols=NULL;
    double *weights=NULL;
    double *sums=NULL;
    fc=fopen("col_idx.bin","rb");
    fw=fopen("dval.bin","rb");
    cols=malloc(CSR_GROUP_CHUNK*sizeof(*cols));
    weights=malloc(CSR_GROUP_CHUNK*sizeof(*weights));
    sums=malloc((size_t)a->nsens*sizeof(*sums));
    if (!fc||!fw||!cols||!weights||!sums) { a->failed=1; goto done; }
    int64_t first=a->rowptr[a->begin];
    if (first > INT64_MAX/(int64_t)sizeof(int64_t)) { a->failed=1; goto done; }
    off_t off=(off_t)(first*(int64_t)sizeof(int64_t));
    if ((int64_t)off != first*(int64_t)sizeof(int64_t) ||
        fseeko(fc,off,SEEK_SET) || fseeko(fw,off,SEEK_SET)) { a->failed=1; goto done; }
    int64_t cursor=first, last=a->rowptr[a->end];
    size_t buffered=0,pos=0;
    for (int i=a->begin;i<a->end;++i) {
        for (int k=0;k<a->nsens;++k) sums[k]=0.0;
        int64_t row_end=a->rowptr[i+1];
        while (cursor<row_end) {
            if (pos==buffered) {
                int64_t remaining=last-cursor;
                size_t take=remaining<(int64_t)CSR_GROUP_CHUNK?(size_t)remaining:CSR_GROUP_CHUNK;
                if (!take || fread(cols,sizeof(*cols),take,fc)!=take ||
                    fread(weights,sizeof(*weights),take,fw)!=take) { a->failed=1; goto done; }
                pos=0; buffered=take;
            }
            size_t take=buffered-pos;
            if ((int64_t)take>row_end-cursor) take=(size_t)(row_end-cursor);
            for (size_t n=0;n<take;++n) {
                int64_t j=cols[pos+n];
                double w=weights[pos+n];
                if (j<0 || j>=a->ne || !isfinite(w)) { a->failed=1; goto done; }
                for (int k=0;k<a->nsens;++k) sums[k]+=w*a->scaled[k][j];
            }
            pos+=take;
            cursor+=(int64_t)take;
        }
        for (int k=0;k<a->nsens;++k) a->out[k][i]=sums[k];
    }
done:
    if (fc && fclose(fc)) a->failed=1;
    if (fw && fclose(fw)) a->failed=1;
    free(cols);free(weights);free(sums);
    return NULL;
}

void filterSensitivity_bin_buffered_mts_multi_csr(
    const double * const *SensInArr, double * const *SensOutArr,
    int nsens, int ne, long long nnz_total)
{
    (void)nnz_total; /* row_ptr[ne] is authoritative; COO nnz may differ. */
    if (!SensInArr||!SensOutArr||nsens<=0||ne<=0) cg_fail("invalid arguments");
    if (cg_file_items("row_ptr.bin",sizeof(int64_t))!=(int64_t)ne+1 ||
        cg_file_items("dsum.bin",sizeof(double))!=ne) cg_fail("metadata dimensions mismatch");
    int64_t *rowptr=malloc(((size_t)ne+1)*sizeof(*rowptr));
    double *dsum=malloc((size_t)ne*sizeof(*dsum));
    double **scaled=calloc((size_t)nsens,sizeof(*scaled));
    if (!rowptr||!dsum||!scaled) cg_fail("allocation failed");
    cg_read("row_ptr.bin",rowptr,sizeof(*rowptr),(size_t)ne+1);
    cg_read("dsum.bin",dsum,sizeof(*dsum),(size_t)ne);
    if (rowptr[0]!=0) cg_fail("row_ptr[0] != 0");
    for (int i=0;i<ne;++i) {
        if (rowptr[i+1]<=rowptr[i]) cg_fail("empty/invalid CSR row");
        if (!(dsum[i]>0.0)||!isfinite(dsum[i])) cg_fail("invalid dsum");
    }
    int64_t nnz=rowptr[ne];
    if (cg_file_items("col_idx.bin",sizeof(int64_t))!=nnz ||
        cg_file_items("dval.bin",sizeof(double))!=nnz) cg_fail("CSR nnz mismatch");
    for (int k=0;k<nsens;++k) {
        if (!SensInArr[k]||!SensOutArr[k]) cg_fail("null sensitivity array");
        for (int j=0;j<nsens;++j)
            if ((const void*)SensOutArr[k]==(const void*)SensInArr[j])
                cg_fail("input/output aliasing is not supported");
        scaled[k]=malloc((size_t)ne*sizeof(double));
        if (!scaled[k]) cg_fail("scaled allocation failed");
        for (int i=0;i<ne;++i) {
            if (!isfinite(SensInArr[k][i])) cg_fail("nonfinite sensitivity input");
            scaled[k][i]=SensInArr[k][i]/dsum[i];
        }
    }
    int nt=4;
    const char *env=getenv("OMP_NUM_THREADS");
    if (env && *env) {
        char *end=NULL; long n=strtol(env,&end,10);
        if (end!=env && *end=='\0' && n>0 && n<=INT_MAX) nt=(int)n;
    }
    if (nt>ne) nt=ne;
    pthread_t *threads=calloc((size_t)nt,sizeof(*threads));
    CSRGroupArgs *args=calloc((size_t)nt,sizeof(*args));
    int *bounds=malloc(((size_t)nt+1)*sizeof(*bounds));
    if (!threads||!args||!bounds) cg_fail("thread allocation failed");
    bounds[0]=0;bounds[nt]=ne;
    for (int t=1;t<nt;++t) {
        int64_t target=(nnz/nt)*t+(nnz%nt)*t/nt;
        int lo=bounds[t-1],hi=ne;
        while (lo<hi) { int mid=lo+(hi-lo)/2; if (rowptr[mid]<target) lo=mid+1; else hi=mid; }
        bounds[t]=lo;
    }
    int started=0;
    for (int t=0;t<nt;++t) {
        args[t]=(CSRGroupArgs){.id=t,.begin=bounds[t],.end=bounds[t+1],
            .ne=ne,.nsens=nsens,.rowptr=rowptr,.scaled=scaled,.out=SensOutArr,.failed=0};
        int rc=pthread_create(&threads[t],NULL,cg_worker,&args[t]);
        if (rc) {
            fprintf(stderr,"pthread_create: %s\n",strerror(rc));
            for (int j=0;j<started;++j) pthread_join(threads[j],NULL);
            cg_fail("pthread_create failed");
        }
        ++started;
    }
    int failed=0;
    for (int t=0;t<started;++t) {
        int rc=pthread_join(threads[t],NULL);
        if (rc) { fprintf(stderr,"pthread_join: %s\n",strerror(rc));failed=1; }
        failed|=args[t].failed;
    }
    for (int k=0;k<nsens;++k) free(scaled[k]);
    free(scaled);free(rowptr);free(dsum);free(threads);free(args);free(bounds);
    if (failed) cg_fail("CSR worker I/O or data failure");
}

/* Drop-in 3-vector CG wrapper. */
void filterSensitivity_bin_buffered_mts3_csr(
    const double *dCGx, const double *dCGy, const double *dCGz,
    double *dCGxFiltered, double *dCGyFiltered, double *dCGzFiltered,
    int ne, long long nnz_total)
{
    const double *in[3]={dCGx,dCGy,dCGz};
    double *out[3]={dCGxFiltered,dCGyFiltered,dCGzFiltered};
    filterSensitivity_bin_buffered_mts_multi_csr(in,out,3,ne,nnz_total);
}
