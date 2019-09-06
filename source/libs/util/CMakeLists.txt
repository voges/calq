project("util")

set(source_files
        fasta-file-reader.cc
        file-reader.cc
        # file-writer.cc
        file.cc
        helpers.cc
        sam-file-reader.cc
        sam-record.cc
)

add_library(util ${source_files})
