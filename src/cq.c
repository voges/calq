#include "cqcodec.h"
#include "cqconfig.h"
#include "cqlib.h"
#include "str.h"
#include <getopt.h>
#include <signal.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

typedef enum {
    CQ_MODE_COMPRESS,
    CQ_MODE_DECOMPRESS,
    CQ_MODE_INFO
} cq_mode_t;

// global variables
bool cq_stats = false;
cq_mode_t cq_mode = CQ_MODE_COMPRESS;

// local options/flags from getopt
static const char *opt_fname_in = NULL;
static const char *opt_fname_out = NULL;
static size_t opt_blocksz = 0;
static bool opt_force = false;

static void print_version(void)
{
    printf("calq %d.%d.%d\n", CQ_VERSION_MAJOR, CQ_VERSION_MINOR, CQ_VERSION_PATCH);
}

static void print_copyright(void)
{
    printf("Copyright (c) 2016\n");
    printf("Leibniz Universitaet Hannover, Institut fuer "
           "Informationsverarbeitung (TNT)\n");
    printf("Contact: Jan Voges <voges@tnt.uni-hannover.de>\n");
}

static void print_help(void)
{
    printf("\n");
    printf("Usage:\n");
    printf("  Compress  : calq [-o FILE] [-b SIZE] [-f] file.sam\n");
    printf("  Decompress: calq -d [-o FILE] [-f] file.cq\n");
    printf("  Info      : calq -i file.cq\n");
    printf("\n");
    printf("Options:\n");
    printf("  -b  --blocksz=SIZE Specify block SIZE\n");
    printf("  -d  --decompress   Decompress\n");
    printf("  -f, --force        Force overwriting of output file(s)\n");
    printf("  -h, --help         Print this help\n");
    printf("  -i, --info         Print information about CQ file\n");
    printf("  -o, --output=FILE  Specify output FILE\n");
    printf("  -v, --version      Display program version\n");
    printf("\n");
}

static void parse_options(int argc, char *argv[])
{
    int opt;

    static struct option long_options[] = {
        { "blocksz",    required_argument, NULL, 'b'},
        { "decompress", no_argument,       NULL, 'd'},
        { "force",      no_argument,       NULL, 'f'},
        { "help",       no_argument,       NULL, 'h'},
        { "info",       no_argument,       NULL, 'i'},
        { "output",     required_argument, NULL, 'o'},
        { "version",    no_argument,       NULL, 'v'},
        { NULL,         0,                 NULL,  0 }
    };

    const char *short_options = "b:dfhio:v";

    do {
        int opt_idx = 0;
        opt = getopt_long(argc, argv, short_options, long_options, &opt_idx);
        switch (opt) {
        case -1:
            break;
        case 'b':
            if (atoi(optarg) <= 0) cq_error("Block size must be positive\n");
            else opt_blocksz = (size_t)atoi(optarg);
            break;
        case 'd':
            if (cq_mode == CQ_MODE_INFO)
                cq_error("Cannot decompress and get info at once\n");
            else
                cq_mode = CQ_MODE_DECOMPRESS;
            break;
        case 'f':
            opt_force = true;
            break;
        case 'h':
            print_version();
            print_copyright();
            print_help();
            exit(EXIT_SUCCESS);
            break;
        case 'i':
            if (cq_mode == CQ_MODE_DECOMPRESS)
                cq_error("Cannot decompress and get info at once\n");
            else
                cq_mode = CQ_MODE_INFO;
            break;
        case 'o':
            opt_fname_out = optarg;
            break;
        case 'v':
            print_version();
            print_copyright();
            exit(EXIT_SUCCESS);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    } while (opt != -1);

    // the input file must be the one remaining command line argument
    if (argc - optind > 1)
        cq_error("Only one input file allowed\n");
    else if (argc - optind < 1)
        cq_error("Input file missing\n");
    else
        opt_fname_in = argv[optind];
}

static const char * fname_extension(const char *path)
{
    const char *dot = strrchr(path, '.');
    if (!dot || dot == path) { return ""; }
    return (dot + 1);
}

static void handle_signal(int sig)
{
    signal(sig, SIG_IGN); // ignore the signal
    cq_log("Catched signal: %d\n", sig);
    signal(sig, SIG_DFL); // invoke default signal action
    raise(sig);
}

int main(int argc, char *argv[])
{
    str_t *fname_in = str_new();
    str_t *fname_out = str_new();

    // register custom signal handler(s)
    signal(SIGHUP,  handle_signal);
    signal(SIGQUIT, handle_signal);
    signal(SIGABRT, handle_signal);
    signal(SIGPIPE, handle_signal);
    signal(SIGTERM, handle_signal);
    signal(SIGXCPU, handle_signal);
    signal(SIGXFSZ, handle_signal);

    // parse command line options and check them for sanity
    parse_options(argc, argv);

    if (cq_mode == CQ_MODE_COMPRESS) {
        // all possible options are legal in compression mode
        if (!opt_blocksz) {
            cq_log("Using default block size 10,000\n");
            opt_blocksz = 10000; // Default value
        }
    } else if (cq_mode == CQ_MODE_DECOMPRESS){
        // option -b is illegal in decompression mode
        if (opt_blocksz)
            cq_error("Illegal option(s) detected\n");
    } else { // CQ_MODE_INFO
        // options -bf are illegal in info mode
        if (opt_blocksz || opt_force)
            cq_error("Illegal option(s) detected\n");
    }

    // check if input file is accessible
    str_copy_cstr(fname_in, opt_fname_in);
    if (access(fname_in->s, F_OK | R_OK))
        cq_error("Cannot access input file: %s\n", fname_in->s);

    if (cq_mode == CQ_MODE_COMPRESS) {
        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "sam"))
            cq_error("Input file extension must be 'sam'\n");

        // create correct output file name
        if (opt_fname_out == NULL) {
            str_copy_str(fname_out, fname_in);
            str_append_cstr(fname_out, ".cq");
        } else {
            str_copy_cstr(fname_out, opt_fname_out);
        }

        // check if output file is accessible
        if (!access(fname_out->s, F_OK | W_OK) && opt_force == false) {
            cq_log("Output file already exists: %s\n", fname_out->s);
            cq_log("Do you want to overwrite %s? ", fname_out->s);
            if (yesno()) ; // proceed
            else exit(EXIT_SUCCESS);
        }

        // invoke compressor
        FILE *fp_in = cq_fopen(fname_in->s, "r");
        FILE *fp_out = cq_fopen(fname_out->s, "wb");
        cq_log("Compressing: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, fp_out, opt_blocksz);
        cqcodec_encode(cqcodec);
        cqcodec_free(cqcodec);
        cq_log("Finished: %s\n", fname_out->s);
        cq_fclose(fp_in);
        cq_fclose(fp_out);
    } else if (cq_mode == CQ_MODE_DECOMPRESS) {
        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "cq"))
            cq_error("Input file extension must be 'cq'\n");

        // create correct output file name
        if (opt_fname_out == NULL) {
            str_copy_str(fname_out, fname_in);
            str_trunc(fname_out, 3); // strip '.cq'
            if (strcmp(fname_extension(fname_out->s), "sam")) 
                str_append_cstr(fname_out, ".sam");
        } else {
            str_copy_cstr(fname_out, opt_fname_out);
        }

        // check if output file is accessible
        if (!access(fname_out->s, F_OK | W_OK) && opt_force == false) {
            cq_log("Output file already exists: %s\n", fname_out->s);
            cq_log("Do you want to overwrite %s? ", fname_out->s);
            if (yesno()) ; // proceed
            else exit(EXIT_SUCCESS);
        }

        // invoke decompressor
        FILE *fp_in = cq_fopen(fname_in->s, "rb");
        FILE *fp_out = cq_fopen(fname_out->s, "w");
        cq_log("Decompressing: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, fp_out, 0);
        cqcodec_decode(cqcodec);
        cqcodec_free(cqcodec);
        cq_log("Finished: %s\n", fname_out->s);
        cq_fclose(fp_in);
        cq_fclose(fp_out);
    } else { // CQ_MODE_INFO
        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "cq"))
            cq_error("Input file extension must be 'cq'\n");
        
        // invoke information tool
        FILE *fp_in = cq_fopen(fname_in->s, "rb");
        cq_log("Reading information: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, NULL, 0);
        cqcodec_info(cqcodec);
        cqcodec_free(cqcodec);
        cq_fclose(fp_in);
    }

    str_free(fname_in);
    str_free(fname_out);

    return EXIT_SUCCESS;
}

