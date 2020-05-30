/* assert */
#include <assert.h>
/* INFINITY */
#include <math.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, calloc, free */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <time.h>
#include <mpi.h>
#include <memory.h>

#define ROWMJR(R, C, NR, NC) (R*NC+C)
#define COLMJR(R, C, NR, NC) (C*NR+R)
/* define access directions for matrices */
#define a(R, C) a[ROWMJR(R,C,ln,n)]
#define b(R, C) b[ROWMJR(R,C,nn,n)]

#define MAIN_PROCESS 0
#define SEND_NUM_TAG 0
#define DISPLS_TAG 1
#define ELEMNTS_TAG 2
#define WEIGHT_TAG 3
#define COUNTS_TAG 4

static void cal(int ** displs, int ** elements, int processes, int size){
    int local, reminder;
    if (size >= processes){
        local = size / processes;
        reminder = size % processes;

        //calculate numbers to send to each node.
        *elements = malloc(processes * sizeof(**elements));
        for (int i = 0; i < processes; i++) {
            *(*elements + i) = local;
        }
        *(*elements + processes - 1) += reminder;
        
    } else {//If there are more rows than available nodes, separate into almost equal pieces.
        local = 1;
        *elements = malloc(processes * sizeof(**elements));
        int i = 0;
        for (; i < size; i++) {
            *(*elements + i) = local;
        }
        for (; i < processes; i++) {
            *(*elements + i) = 0;
        }
    }
    //calculate offset that points from send buffer(vals).
    *displs = malloc(processes * sizeof(**displs));
    for (int i = 0; i < processes; i++) {
        *(*displs + i) = 0;
        for (int j = 0; j < i; j++) {
            *(*displs + i) += *(*elements + j);
        }
    }
}
static void
load(
    char const *const filename,
    int *const np,
    float **const ap, 
    int processes,
    int ** displs, 
    int ** elements,
    int rank
) {
    int n;
    float *a = NULL;
    if (rank == MAIN_PROCESS) {//Main node read the file and send to other nodes piece by piece.
        int i, j, k, ret;
        FILE *fp = NULL;


        /* open the file */
        fp = fopen(filename, "r");
        assert(fp);

        /* get the number of nodes in the graph */
        ret = fscanf(fp, "%d", &n);
        assert(1 == ret);

        //Calculate how many rows each node will hold. And their offsets(displs).
        cal(displs, elements, processes, n);

        /* allocate memory for local values */
        a = malloc(n * *(*elements) * sizeof(*a));
        for (j = 0; j < *(*elements) * n; ++j) {
            ret = fscanf(fp, "%f", &a[j]);
            assert(1 == ret);
        }
        *ap = a;

        for (i = 1; i < processes; ++i) {
            a = malloc(n * *(*elements + i) * sizeof(*a));
            //Send just collected info to that node
            MPI_Send(&n, 1, MPI_INTEGER, i, SEND_NUM_TAG, MPI_COMM_WORLD);

            MPI_Send(*displs, processes, MPI_INTEGER, i, DISPLS_TAG, MPI_COMM_WORLD);
            MPI_Send(*elements, processes, MPI_INTEGER, i, ELEMNTS_TAG, MPI_COMM_WORLD);
            //Read file
            for (j = 0; j < *(*elements + i) * n; ++j) {
                ret = fscanf(fp, "%f", &a[j]);
                assert(1 == ret);
            }
            //printf("%d. %d info read\n", i, j);//TODO debug use only
            MPI_Send(&j, 1, MPI_INTEGER, i, COUNTS_TAG, MPI_COMM_WORLD);
            MPI_Send(a, j, MPI_FLOAT, i, WEIGHT_TAG, MPI_COMM_WORLD);
            free(a);//Free memory after send.
        }

        /* close file */
        ret = fclose(fp);
        assert(!ret);
    } else {//All nodes except MAIN NODE will receive piece(rows) of graph data.
        int count;
        MPI_Recv(&n, 1, MPI_INTEGER, MAIN_PROCESS, SEND_NUM_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *displs = malloc(processes * sizeof(**displs));
        MPI_Recv(*displs, processes, MPI_INTEGER, MAIN_PROCESS, DISPLS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *elements = malloc(processes * sizeof(**elements));
        MPI_Recv(*elements, processes, MPI_INTEGER, MAIN_PROCESS, ELEMNTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&count, 1, MPI_INTEGER, MAIN_PROCESS, COUNTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        a = malloc(count * sizeof(*a));
        MPI_Recv(a, count, MPI_FLOAT, MAIN_PROCESS, WEIGHT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *ap = a;
    }
    /* record output values */

    *np = n;
    MPI_Barrier(MPI_COMM_WORLD);
}

static void
dijkstra(
    int const s,
    int const n,
    float const *const a,
    float **const lp, 
    int rank, 
    int * displs, 
    int * elements,
    int processes
) {
    int i, j, k, sourceNode = 0;
    struct float_int {
        float distance;
        int u;
    } min;
    char *m = NULL;
    float *l = NULL;
    float * ll= NULL;


    m = calloc(n, sizeof(*m));
    assert(m);

    l = malloc(n * sizeof(*l));
    assert(l);

    ll = malloc(n * sizeof(*l));
    assert(ll);

    for (i = 0; i < processes; i++){
        if (s < displs[i]){
            sourceNode = i - 1;
            break;
        }
    }
    
    //The initial result distance is what currently the distance between source to each vertex.
    if (rank == sourceNode) {//Copy and prepare for broad cast.
        for (i = 0; i < n; ++i) {
            l[i] = a[i + n * (s - displs[sourceNode])];
            
        }
    }
    MPI_Bcast(l, n, MPI_FLOAT, sourceNode, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   
    m[s] = 1;
    min.u = -1; /* avoid compiler warning */

    //Iterate through each.
    for (i = 1; i < n; ++i) {
        min.distance = INFINITY;

        for (j = 0; j < n; ++j) {
            if (!m[j] && l[j] < min.distance) {
                min.distance = l[j];
                min.u = j;
            }
            ll[j] = l[j];
        }
       
        m[min.u] = 1;
        for (j = 0; j < elements[rank]; j++){
            if (m[j + displs[rank]]){
                continue;
            }
            if (a(j, min.u) + min.distance < ll[j + displs[rank]]){
                ll[j + displs[rank]] = a(j, min.u) + min.distance;
            }
        }

        MPI_Allreduce(ll, l, n, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(m);

    *lp = l;
}

static void
print_time(double const seconds) {
    printf("Operation Time: %0.04fs\n", seconds);
}

static void
print_numbers(
    char const *const filename,
    int const n,
    float const *const numbers
){
    int i;
    FILE *fout;

    /* open file */
    if (NULL == (fout = fopen(filename, "w"))) {
        fprintf(stderr, "error opening '%s'\n", filename);
        abort();
    }

    /* write numbers to fout */
    for (i = 0; i < n; ++i) {
        fprintf(fout, "%10.4f\n", numbers[i]);
    }

    fclose(fout);
}

int
main(int argc, char **argv) {
    int n, npes, myrank;
    double ts, te;
    float *a = NULL, *result = NULL;
    int * displs = NULL, * elements = NULL;

    if (argc < 4) {
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <source> <output_file>.\n");
        return EXIT_FAILURE;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    load(argv[1], &n, &a, npes, &displs, &elements, myrank);

    MPI_Barrier(MPI_COMM_WORLD);
    ts = MPI_Wtime();
    dijkstra(atoi(argv[2]), n, a, &result, myrank, displs, elements, npes);
    te = MPI_Wtime();

    print_time((te - ts) / CLOCKS_PER_SEC);
    if (myrank == MAIN_PROCESS) {
        print_time(te - ts);
        print_numbers(argv[3], n, result);
    }
    free(a);
    free(result);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
