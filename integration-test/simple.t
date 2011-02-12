#!/usr/bin/env perl

use File::Basename qw/dirname/;
use POSIX qw/WEXITSTATUS/;
use Test::More tests => 2;

my $dir = dirname(__FILE__);
my $exe = "$dir/../bin/bincoeff";

my $output = `$exe 10 5`;
chomp $output;
ok(!WEXITSTATUS($?), 'command returned success');

is($output, "log(10 choose 5): 5.529429\n10 choose 5: 252", 'output is as expected');

done_testing();
