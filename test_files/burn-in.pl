#!/usr/bin/perl

use strict;
use warnings;

my $tsc = "tsc";
my $dir = ".";
opendir(DIR, $dir) or die $!;

while (my $file = readdir(DIR)) {
    next if ($file =~ m/^\./);
    next if !($file =~ m/.*\.sam$/); # only SAM files
    $file = $dir."/".$file;
    print "Testing tsc with $file ...";

    # Compress and decompress
    my $file_com = $file.".tsc";
    my $file_dec = $file_com.".sam";
    my $com = `$tsc $file -o $file_com -f 1>/dev/null`;
    my $dec = `$tsc -d  $file_com -o $file_dec -f 1>/dev/null`;

    # Check correctness of decoded file.
    my $ret = `diff $file $file_dec`;
    if ($ret eq "") {
        print " passed\n";
    } else {
        print " NOT passed: $ret\n";
        cleanup();
        exit -1;
    }

    # Cleanup
    unlink $file_com or warn "Could not unlink $file_com: $!";
    unlink $file_dec or warn "Could not unlink $file_dec: $!";
}

closedir(DIR);
exit 0;

