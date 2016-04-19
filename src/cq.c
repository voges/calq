#include "cqcodec.h"
#include "cqconfig.h"
#include "cqlib.h"
#include "misc/str.h"
#include <getopt.h>
#include <signal.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

cq_opts_t cq_opts;

static void print_copyright(void)
{
    printf("Copyright (c) 2015-%d\n", CQ_BUILD_YEAR);
    printf("Leibniz Universitaet Hannover, Institut fuer "
           "Informationsverarbeitung (TNT)\n");
    printf("Contact: Jan Voges <voges@tnt.uni-hannover.de>\n");
}

static void print_version(void)
{
    printf("calq %d.%d.%d\n", CQ_VERSION_MAJOR, CQ_VERSION_MINOR, CQ_VERSION_PATCH);
    printf("Build time: %s\n", CQ_UTCTIMESTAMP);
    printf("Git revision: %s\n", CQ_GITREVISION_LONG);
    printf("\n");
    print_copyright();
}

static void print_help(void)
{
    print_version();
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

static int parse_options(int argc, char *argv[])
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
            if (atoi(optarg) <= 0) {
                cq_err("Block size must be positive\n");
                return CQ_FAILURE;
            } else {
                cq_opts.blocksz = (size_t)atoi(optarg);
            }
            break;
        case 'd':
            if (cq_opts.mode == CQ_OPTS_MODE_INFO) {
                cq_err("Cannot decompress and get info at once\n");
                return CQ_FAILURE;
            } else {
                cq_opts.mode = CQ_OPTS_MODE_DECOMPRESS;
            }
            break;
        case 'f':
            cq_opts.force = true;
            break;
        case 'h':
            print_help();
            exit(EXIT_SUCCESS);
            break;
        case 'i':
            if (cq_opts.mode == CQ_OPTS_MODE_DECOMPRESS) {
                cq_err("Cannot decompress and get info at once\n");
                return CQ_FAILURE;
            } else {
                cq_opts.mode = CQ_OPTS_MODE_INFO;
            }
            break;
        case 'o':
            cq_opts.fname_out = optarg;
            break;
        case 'v':
            print_version();
            exit(EXIT_SUCCESS);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    } while (opt != -1);

    // the input file must be the one remaining command line argument
    if (argc - optind > 1) {
        cq_err("Only one input file allowed\n");
        return CQ_FAILURE;
    } else if (argc - optind < 1) {
        cq_err("Input file missing\n");
        return CQ_FAILURE;
    } else {
        cq_opts.fname_in = argv[optind];
    }

    // sanity checks
    if (cq_opts.mode == CQ_OPTS_MODE_COMPRESS) {
        // all possible options are legal in compression mode
        if (!cq_opts.blocksz) {
            cq_out("Using default block size 10,000\n");
            cq_opts.blocksz = 10000; // default value
        }
    } else if (cq_opts.mode == CQ_OPTS_MODE_DECOMPRESS){
        // option -b is illegal in decompression mode
        if (cq_opts.blocksz) {
            cq_err("Illegal option(s) detected\n");
            exit(EXIT_FAILURE);
        }
    } else { // CQ_OPTS_MODE_INFO
        // options -bf are illegal in info mode
        if (cq_opts.blocksz || cq_opts.force) {
            cq_err("Illegal option(s) detected\n");
            exit(EXIT_FAILURE);
        }
    }

    return CQ_SUCCESS;
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
    cq_out("Catched signal: %d\n", sig);
    signal(sig, SIG_DFL); // invoke default signal action
    raise(sig);
}

int main(int argc, char *argv[])
{
    // init CLI options
    cq_opts.blocksz = 0;
    cq_opts.fname_in = NULL;
    cq_opts.fname_out = NULL;
    cq_opts.force = false;
    cq_opts.mode = CQ_OPTS_MODE_COMPRESS;

    // register custom signal handler(s)
    signal(SIGHUP,  handle_signal);
    signal(SIGQUIT, handle_signal);
    signal(SIGABRT, handle_signal);
    signal(SIGPIPE, handle_signal);
    signal(SIGTERM, handle_signal);
    signal(SIGXCPU, handle_signal);
    signal(SIGXFSZ, handle_signal);

    // parse command line options
    if (CQ_SUCCESS != parse_options(argc, argv)) {
        cq_err("Failed to parse options\n");
        exit(EXIT_FAILURE);
    }

    str_t *fname_in = str_new();
    str_t *fname_out = str_new();

    if (cq_opts.mode == CQ_OPTS_MODE_COMPRESS) {
        str_copy_cstr(fname_in, cq_opts.fname_in);

        // check if input file is accessible
        if (access(fname_in->s, F_OK | R_OK)) {
            cq_err("Cannot access input file: %s\n", cq_opts.fname_in);
            goto cleanup;
        }

        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "sam")) {
            cq_err("Input file extension must be 'sam'\n");
            goto cleanup;
        }

        // create correct output file name
        if (cq_opts.fname_out == NULL) {
            str_copy_str(fname_out, fname_in);
            str_append_cstr(fname_out, ".cq");
        } else {
            str_copy_cstr(fname_out, cq_opts.fname_out);
        }

        // check if output file is already there and if the user wants to
        // overwrite it in this case
        if (!access(fname_out->s, F_OK | W_OK) && cq_opts.force == false) {
            cq_out("Output file already exists: %s\n", fname_out->s);
            cq_out("Do you want to overwrite %s? ", fname_out->s);
            if (!yesno()) goto cleanup;;
        }

        // invoke compressor
        FILE *fp_in = cq_fopen(fname_in->s, "r");
        FILE *fp_out = cq_fopen(fname_out->s, "wb");
        cq_out("Compressing: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, fp_out, cq_opts.blocksz);
        if (CQ_SUCCESS != cqcodec_encode(cqcodec)) {
            cq_err("Encoding failed\n");
        } else {
            cq_out("Finished: %s\n", fname_out->s);
        }
        cqcodec_delete(cqcodec);
        cq_fclose(fp_in);
        cq_fclose(fp_out);
    } else if (cq_opts.mode == CQ_OPTS_MODE_DECOMPRESS) {
        str_copy_cstr(fname_in, cq_opts.fname_in);

        // check if input file is accessible
        if (access(fname_in->s, F_OK | R_OK)) {
            cq_err("Cannot access input file: %s\n", cq_opts.fname_in);
            goto cleanup;
        }

        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "cq")) {
            cq_err("Input file extension must be 'cq'\n");
            goto cleanup;
        }

        // create correct output file name
        if (cq_opts.fname_out == NULL) {
            str_copy_str(fname_out, fname_in);
            str_append_cstr(fname_out, ".sam");
            //str_trunc(fname_out, 3); // strip '.cq'
            //if (strcmp(fname_extension(fname_out->s), "sam")) 
            //    str_append_cstr(fname_out, ".sam");
        } else {
            str_copy_cstr(fname_out, cq_opts.fname_out);
        }

        // check if output file is accessible
        if (!access(fname_out->s, F_OK | W_OK) && cq_opts.force == false) {
            cq_out("Output file already exists: %s\n", fname_out->s);
            cq_out("Do you want to overwrite %s? ", fname_out->s);
            if (!yesno()) goto cleanup;
        }

        // invoke decompressor
        FILE *fp_in = cq_fopen(fname_in->s, "rb");
        FILE *fp_out = cq_fopen(fname_out->s, "w");
        cq_out("Decompressing: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, fp_out, 0);
        if (CQ_SUCCESS != cqcodec_decode(cqcodec)) {
            cq_err("Decoding failed\n");
        } else {
            cq_out("Finished: %s\n", fname_out->s);
        }
        cqcodec_delete(cqcodec);
        cq_fclose(fp_in);
        cq_fclose(fp_out);
    } else { // CQ_OPTS_MODE_INFO
        str_copy_cstr(fname_in, cq_opts.fname_in);

        // check if input file is accessible
        if (access(fname_in->s, F_OK | R_OK)) {
            cq_err("Cannot access input file: %s\n", cq_opts.fname_in);
            goto cleanup;
        }

        // check correct file name extension of input file
        if (strcmp(fname_extension(fname_in->s), "cq")) {
            cq_err("Input file extension must be 'cq'\n");
            goto cleanup;
        }

        // invoke information tool
        FILE *fp_in = cq_fopen(fname_in->s, "rb");
        cq_out("Reading information: %s\n", fname_in->s);
        cqcodec_t *cqcodec = cqcodec_new(fp_in, NULL, 0);
        if (CQ_SUCCESS != cqcodec_info(cqcodec)) {
            cq_err("Info failed\n");
        } else {
            cq_out("Finished: %s\n", fname_in->s);
        }
        cqcodec_delete(cqcodec);
        cq_fclose(fp_in);
    }

cleanup:
    str_free(fname_in);
    str_free(fname_out);
    return EXIT_SUCCESS;
}

